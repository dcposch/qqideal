"""Polynomial rings and polynomials, over Q or a prime field.

Every call into python-flint is made from this module and nowhere else, so an
API break upstream is one file to fix. Rings carry the graded reverse
lexicographic order, because that is the only order msolve prints.

Parsing does not use ``eval``. A source string is first handed to
:func:`msolveio.emit_system`, which rejects anything msolve would mis-parse
(parentheses, ``x/2``, unknown identifiers, repeated monomials), and only then
re-tokenized into exponent vectors.

Prime fields are backed by ``fmpq_mpoly`` too: coefficients live in a rational
polynomial, but every :class:`Poly` on a prime-field ring is normalized so each
coefficient is the canonical representative in ``0 .. p-1``. Arithmetic is
therefore genuinely arithmetic in ``F_p`` -- ``p*x`` is zero, not a nonzero
rational polynomial -- and emission to msolve is a straight print.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable, Iterator, Sequence

from flint import Ordering, fmpq, fmpq_mpoly, fmpq_mpoly_ctx
from msolveio import emit_system

from .errors import MsolveInputError, RingMismatch

__all__ = ["Ring", "Poly"]

_TOKEN_RE = re.compile(r"[A-Za-z][A-Za-z0-9]*|[0-9]+|[+\-*/^]|\s+")

Scalar = int | Fraction


@dataclass(frozen=True, init=False, repr=False)
class Ring:
    """A polynomial ring over Q (``characteristic=0``) or over ``F_p``.

    ``Ring("x", "y")`` is ``QQ[x, y]`` with the grevlex order. Variable names
    match ``[A-Za-z][A-Za-z0-9]*``: no underscores, because msolve's parser does
    not accept them.

    :param names: the variable names, in order.
    :param characteristic: ``0`` for Q, or a prime below ``2**31``.
    :raises MsolveInputError: for an illegal name or a non-prime characteristic.
    """

    names: tuple[str, ...]
    characteristic: int

    def __init__(self, *names: str, characteristic: int = 0) -> None:
        name_tuple = tuple(names)
        # emit_system owns the rules for what msolve will accept as a variable
        # name and as a characteristic; ask it rather than re-deriving them.
        emit_system(["0"], variables=name_tuple, characteristic=characteristic)
        object.__setattr__(self, "names", name_tuple)
        object.__setattr__(self, "characteristic", characteristic)
        object.__setattr__(
            self, "_ctx", fmpq_mpoly_ctx.get(name_tuple, Ordering.degrevlex)
        )

    # -- construction ----------------------------------------------------

    def __call__(self, source: str) -> Poly:
        """Parse ``source`` into a :class:`Poly` of this ring.

        The grammar is msolve's: an expanded sum of monomials, no parentheses,
        with rational coefficients written as a leading ``1/2*x``. Rational
        coefficients are accepted on prime-field rings too and are reduced
        immediately; a denominator divisible by the characteristic is an error.

        :raises MsolveInputError: if ``source`` is not a legal polynomial.
        """
        if not isinstance(source, str):
            raise MsolveInputError(
                f"expected a polynomial string, got {type(source).__name__}"
            )
        # Validate against msolve's own grammar first. characteristic=0 here on
        # purpose: '1/2*x' is a legal way to write a prime-field element, and
        # the reduction happens below.
        emit_system([source], variables=self.names, characteristic=0)
        return self._wrap(self._ctx.from_dict(_parse(source, self.names)))

    def constant(self, value: Scalar) -> Poly:
        """The constant polynomial ``value``."""
        return self._wrap(self._ctx.from_dict({(0,) * self.nvars: _to_fmpq(value)}))

    def gens(self) -> tuple[Poly, ...]:
        """The variables, as polynomials."""
        return tuple(self._wrap(g) for g in self._ctx.gens())

    def gen(self, index: int) -> Poly:
        """The ``index``-th variable, as a polynomial."""
        return self._wrap(self._ctx.gen(index))

    @property
    def nvars(self) -> int:
        """How many variables this ring has."""
        return len(self.names)

    def extend(self, *extra: str) -> Ring:
        """This ring with ``extra`` variables appended, same characteristic."""
        return Ring(*self.names, *extra, characteristic=self.characteristic)

    def fresh_variable(self, prefix: str = "u") -> str:
        """A variable name not already used here: ``u``, else ``u1``, ``u2``..."""
        if prefix not in self.names:
            return prefix
        index = 1
        while f"{prefix}{index}" in self.names:
            index += 1
        return f"{prefix}{index}"

    def embed(self, poly: Poly) -> Poly:
        """Re-express ``poly`` in this ring, matching variables by name.

        Every variable of ``poly.ring`` must be a variable of this ring, and the
        characteristics must agree.

        :raises RingMismatch: if the embedding does not exist.
        """
        source = poly.ring
        if source.characteristic != self.characteristic:
            raise RingMismatch(
                f"cannot embed a polynomial of {source} into {self}: "
                f"different characteristic"
            )
        try:
            positions = [self.names.index(name) for name in source.names]
        except ValueError:
            missing = [n for n in source.names if n not in self.names]
            raise RingMismatch(
                f"cannot embed a polynomial of {source} into {self}: "
                f"{self} has no variable(s) {', '.join(missing)}"
            ) from None
        terms: dict[tuple[int, ...], fmpq] = {}
        for exponents, coeff in poly.terms():
            key = [0] * self.nvars
            for position, exponent in zip(positions, exponents):
                key[position] = exponent
            terms[tuple(key)] = coeff
        return self._wrap(self._ctx.from_dict(terms))

    # -- elimination-free helpers ----------------------------------------

    def resultant(self, f: Poly, g: Poly, var: str) -> Poly:
        """The resultant of ``f`` and ``g`` with respect to ``var``.

        The result is free of ``var`` -- it lies in the subring generated by the
        remaining variables -- but is returned as a :class:`Poly` of this same
        ring, so it composes with everything else here. It is a constant
        polynomial when ``var`` was the only variable.

        :raises NotImplementedError: on a prime-field ring. Reducing a resultant
            computed over Q is not the resultant over ``F_p`` when the leading
            coefficient in ``var`` vanishes mod ``p``, and v0.1 will not guess
            which case it is in.
        """
        self._require_characteristic_zero("resultant")
        raw = self._check(f)._raw.resultant(self._check(g)._raw, self._check_var(var))
        return self._wrap(raw)

    def discriminant(self, f: Poly, var: str) -> Poly:
        """The discriminant of ``f`` with respect to ``var``.

        :raises NotImplementedError: on a prime-field ring, for the reason given
            in :meth:`resultant`.
        """
        self._require_characteristic_zero("discriminant")
        return self._wrap(self._check(f)._raw.discriminant(self._check_var(var)))

    # -- internals -------------------------------------------------------

    def _wrap(self, raw: fmpq_mpoly) -> Poly:
        if self.characteristic:
            raw = _reduce_mod_p(raw, self.characteristic, self._ctx, self.nvars)
        return Poly(self, raw)

    def _check(self, poly: Poly) -> Poly:
        if not isinstance(poly, Poly):
            raise RingMismatch(f"expected a Poly, got {type(poly).__name__}")
        if poly.ring != self:
            raise RingMismatch(f"polynomial belongs to {poly.ring}, not {self}")
        return poly

    def _check_var(self, var: str) -> str:
        if var not in self.names:
            raise MsolveInputError(f"{var!r} is not a variable of {self}")
        return var

    def _require_characteristic_zero(self, operation: str) -> None:
        if self.characteristic:
            raise NotImplementedError(
                f"{operation} is implemented over Q only; this ring has "
                f"characteristic {self.characteristic}"
            )

    def __str__(self) -> str:
        field = "QQ" if self.characteristic == 0 else f"GF({self.characteristic})"
        return f"{field}[{', '.join(self.names)}]"

    def __repr__(self) -> str:
        args = ", ".join(repr(name) for name in self.names)
        if self.characteristic:
            return f"Ring({args}, characteristic={self.characteristic})"
        return f"Ring({args})"


class Poly:
    """A polynomial of a :class:`Ring`. Immutable.

    Build these with ``ring("x^2 - y")`` or from :meth:`Ring.gens`; the
    two-argument constructor is internal.
    """

    __slots__ = ("_ring", "_raw")

    def __init__(self, ring: Ring, raw: fmpq_mpoly) -> None:
        self._ring = ring
        self._raw = raw

    @property
    def ring(self) -> Ring:
        """The ring this polynomial belongs to."""
        return self._ring

    def degree(self) -> int:
        """Total degree. The zero polynomial has degree ``-1``."""
        return int(self._raw.total_degree())

    def is_zero(self) -> bool:
        """Whether this is the zero polynomial."""
        return bool(self._raw.is_zero())

    def is_one(self) -> bool:
        """Whether this is the constant polynomial ``1``."""
        return bool(self._raw.is_one())

    def is_constant(self) -> bool:
        """Whether this polynomial involves no variables."""
        return bool(self._raw.is_constant())

    def terms(self) -> tuple[tuple[tuple[int, ...], fmpq], ...]:
        """The terms as ``(exponent vector, coefficient)``, grevlex-descending."""
        return tuple(
            (tuple(int(e) for e in monom), coeff)
            for monom, coeff in zip(self._raw.monoms(), self._raw.coeffs())
        )

    def leading_monomial(self) -> tuple[int, ...]:
        """The grevlex-leading monomial, as an exponent vector.

        :raises ValueError: for the zero polynomial, which has no leading term.
        """
        monoms = self._raw.monoms()
        if not monoms:
            raise ValueError("the zero polynomial has no leading monomial")
        return tuple(int(e) for e in monoms[0])

    def to_msolve(self) -> str:
        """Render as msolve input text: an expanded sum of monomials."""
        return _msolve_text(self.terms(), self._ring.names)

    # -- arithmetic ------------------------------------------------------

    def __add__(self, other: Poly | Scalar) -> Poly:
        return self._ring._wrap(self._raw + self._coerce(other)._raw)

    __radd__ = __add__

    def __sub__(self, other: Poly | Scalar) -> Poly:
        return self._ring._wrap(self._raw - self._coerce(other)._raw)

    def __rsub__(self, other: Poly | Scalar) -> Poly:
        return self._ring._wrap(self._coerce(other)._raw - self._raw)

    def __mul__(self, other: Poly | Scalar) -> Poly:
        return self._ring._wrap(self._raw * self._coerce(other)._raw)

    __rmul__ = __mul__

    def __neg__(self) -> Poly:
        return self._ring._wrap(-self._raw)

    def __pos__(self) -> Poly:
        return self

    def __pow__(self, exponent: int) -> Poly:
        if not isinstance(exponent, int) or isinstance(exponent, bool):
            raise TypeError(f"exponent must be an int, got {type(exponent).__name__}")
        if exponent < 0:
            raise NotImplementedError(
                "negative powers are not polynomials; qqideal has no fraction field"
            )
        return self._ring._wrap(self._raw**exponent)

    def _coerce(self, other: Poly | Scalar) -> Poly:
        if isinstance(other, Poly):
            return self._ring._check(other)
        if isinstance(other, (int, Fraction)) and not isinstance(other, bool):
            return self._ring.constant(other)
        raise RingMismatch(
            f"cannot combine a Poly with {type(other).__name__}; use a Poly of "
            f"the same ring, an int, or a Fraction"
        )

    # -- identity --------------------------------------------------------

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Poly):
            return NotImplemented
        return self._ring == other._ring and self._raw == other._raw

    def __hash__(self) -> int:
        return hash((self._ring, tuple((m, str(c)) for m, c in self.terms())))

    def __str__(self) -> str:
        return str(self._raw)

    def __repr__(self) -> str:
        return f"Poly({str(self._raw)!r}, {self._ring!r})"


def _to_fmpq(value: Scalar) -> fmpq:
    if isinstance(value, bool):
        raise MsolveInputError("bool is not a coefficient")
    if isinstance(value, int):
        return fmpq(value)
    if isinstance(value, Fraction):
        return fmpq(value.numerator, value.denominator)
    raise MsolveInputError(
        f"coefficient must be an int or a Fraction, got {type(value).__name__}"
    )


def _reduce_mod_p(
    raw: fmpq_mpoly, characteristic: int, ctx: fmpq_mpoly_ctx, nvars: int
) -> fmpq_mpoly:
    """Put every coefficient in ``0 .. p-1``, dropping the ones that vanish."""
    reduced: dict[tuple[int, ...], fmpq] = {}
    for monom, coeff in zip(raw.monoms(), raw.coeffs()):
        numerator = int(coeff.numer()) % characteristic
        denominator = int(coeff.denom()) % characteristic
        if denominator == 0:
            raise MsolveInputError(
                f"coefficient {coeff} has a denominator divisible by the "
                f"characteristic {characteristic}, so it is not a field element"
            )
        value = numerator * pow(denominator, -1, characteristic) % characteristic
        if value:
            reduced[tuple(int(e) for e in monom)] = fmpq(value)
    return ctx.from_dict(reduced)


def _tokenize(source: str) -> list[str]:
    """Split into identifiers, integers, and operators. Whitespace is dropped."""
    tokens: list[str] = []
    position = 0
    for match in _TOKEN_RE.finditer(source):
        if match.start() != position:
            raise MsolveInputError(
                f"unexpected character {source[position]!r} in {source!r}"
            )
        position = match.end()
        token = match.group()
        if not token.isspace():
            tokens.append(token)
    if position != len(source):
        raise MsolveInputError(f"unexpected character {source[position]!r} in {source!r}")
    return tokens


def _parse(source: str, names: Sequence[str]) -> dict[tuple[int, ...], fmpq]:
    """Parse an expanded monomial sum into ``{exponent vector: coefficient}``.

    Assumes :func:`msolveio.emit_system` has already accepted ``source``; the
    checks here are a second line of defence, not the first.
    """
    tokens = _tokenize(source)
    if not tokens:
        raise MsolveInputError(f"empty polynomial: {source!r}")
    index = {name: position for position, name in enumerate(names)}
    nvars = len(names)

    terms: dict[tuple[int, ...], fmpq] = {}
    position = 0
    count = len(tokens)
    while position < count:
        sign = 1
        if tokens[position] in ("+", "-"):
            sign = -1 if tokens[position] == "-" else 1
            position += 1
        coefficient = fmpq(1)
        exponents = [0] * nvars
        leading_factor = True
        while True:
            if position >= count:
                raise MsolveInputError(f"polynomial ends mid-term: {source!r}")
            token = tokens[position]
            if token.isdigit():
                if not leading_factor:
                    raise MsolveInputError(
                        f"numeric coefficient {token} must come first in its "
                        f"term: {source!r}"
                    )
                position, coefficient = _parse_coefficient(tokens, position, source)
            elif token[0].isalpha():
                if token not in index:
                    raise MsolveInputError(
                        f"unknown identifier {token!r} in {source!r}"
                    )
                position += 1
                exponent = 1
                if position < count and tokens[position] == "^":
                    position += 1
                    if position >= count or not tokens[position].isdigit():
                        raise MsolveInputError(
                            f"exponent of {token!r} must be a positive integer: "
                            f"{source!r}"
                        )
                    exponent = int(tokens[position])
                    position += 1
                exponents[index[token]] += exponent
            else:
                raise MsolveInputError(f"unexpected {token!r} in {source!r}")
            leading_factor = False
            if position >= count or tokens[position] in ("+", "-"):
                break
            if tokens[position] != "*":
                raise MsolveInputError(
                    f"expected '*' between factors, got {tokens[position]!r} in "
                    f"{source!r}"
                )
            position += 1
        key = tuple(exponents)
        terms[key] = terms.get(key, fmpq(0)) + sign * coefficient
    return {key: value for key, value in terms.items() if value != 0}


def _parse_coefficient(
    tokens: list[str], position: int, source: str
) -> tuple[int, fmpq]:
    numerator = int(tokens[position])
    position += 1
    denominator = 1
    if position < len(tokens) and tokens[position] == "/":
        position += 1
        if position >= len(tokens) or not tokens[position].isdigit():
            raise MsolveInputError(
                f"expected an integer denominator after '/' in {source!r}"
            )
        denominator = int(tokens[position])
        if denominator == 0:
            raise MsolveInputError(f"zero denominator in {source!r}")
        position += 1
    return position, fmpq(numerator, denominator)


def _msolve_text(
    terms: Iterable[tuple[tuple[int, ...], fmpq]], names: Sequence[str]
) -> str:
    """Render terms as msolve input text. The zero polynomial prints as ``0``."""
    pieces: list[str] = []
    for exponents, coefficient in terms:
        numerator = int(coefficient.numer())
        denominator = int(coefficient.denom())
        sign = "-" if numerator < 0 else "+"
        numerator = abs(numerator)
        factors = [
            name if exponent == 1 else f"{name}^{exponent}"
            for name, exponent in zip(names, exponents)
            if exponent
        ]
        if factors and numerator == 1 and denominator == 1:
            body = "*".join(factors)
        else:
            head = str(numerator) if denominator == 1 else f"{numerator}/{denominator}"
            body = "*".join([head, *factors])
        pieces.append(sign + body if pieces else ("-" + body if sign == "-" else body))
    return "".join(pieces) if pieces else "0"
