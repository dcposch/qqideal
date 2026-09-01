"""Krull dimension and 0-dimensional degree, read off a leading ideal.

Both quantities are combinatorial once you have the leading monomials of a
Groebner basis: ``k[x]/I`` and ``k[x]/LT(I)`` have the same Hilbert function, so
they have the same dimension, and in the 0-dimensional case the same vector
space dimension over ``k``. Nothing here talks to msolve; the input is just
exponent vectors, which makes it cheap to test.

The convention for the unit ideal is ``dimension == -1``: ``LT(I)`` contains the
constant monomial, ``k[x]/I`` is the zero ring, and the variety is empty.
"""

from __future__ import annotations

from itertools import product
from typing import Iterable, Sequence

__all__ = ["dimension", "degree", "dim_and_degree"]


def dimension(leading: Iterable[Sequence[int]], nvars: int) -> int:
    """Krull dimension of ``k[x1..xn]/I`` from the leading monomials of ``I``.

    The dimension is the largest number of variables ``S`` such that ``LT(I)``
    contains no monomial supported entirely inside ``S``.

    :param leading: exponent vectors of the leading monomials, each of length
        ``nvars``. An empty collection means the zero ideal.
    :param nvars: the number of variables.
    :returns: the dimension, or ``-1`` for the unit ideal.
    """
    supports = _supports(leading, nvars)
    best = -1

    def search(start: int, chosen: int, size: int) -> None:
        nonlocal best
        # Independence is downward closed, so a dependent set prunes the whole
        # subtree below it.
        if any(support & ~chosen == 0 for support in supports):
            return
        if size > best:
            best = size
        for index in range(start, nvars):
            if size + (nvars - index) <= best:
                break
            search(index + 1, chosen | (1 << index), size + 1)

    search(0, 0, 0)
    return best


def degree(leading: Iterable[Sequence[int]], nvars: int) -> int:
    """Number of standard monomials, i.e. ``dim_k k[x1..xn]/I``.

    This is the affine degree of a 0-dimensional ideal, counted with
    multiplicity.

    :raises ValueError: unless the ideal is 0-dimensional. In positive dimension
        there are infinitely many standard monomials and no such number.
    """
    monomials = _normalize(leading, nvars)
    bounds: list[int] = []
    for index in range(nvars):
        pure = [
            monomial[index]
            for monomial in monomials
            if all(e == 0 for position, e in enumerate(monomial) if position != index)
        ]
        if not pure:
            raise ValueError(
                f"the ideal is not 0-dimensional: the leading ideal contains no "
                f"pure power of variable {index}"
            )
        bounds.append(min(pure))
    if any(bound == 0 for bound in bounds):
        # A pure power with exponent 0 is the constant monomial: the unit ideal.
        return 0

    count = 0
    for candidate in product(*(range(bound) for bound in bounds)):
        if not any(_divides(monomial, candidate) for monomial in monomials):
            count += 1
    return count


def dim_and_degree(
    leading: Iterable[Sequence[int]], nvars: int
) -> tuple[int, int | None]:
    """``(dimension, degree)``, with ``degree`` only when the dimension is 0."""
    monomials = _normalize(leading, nvars)
    dim = dimension(monomials, nvars)
    if dim != 0:
        return dim, None
    return dim, degree(monomials, nvars)


def _normalize(
    leading: Iterable[Sequence[int]], nvars: int
) -> tuple[tuple[int, ...], ...]:
    monomials: list[tuple[int, ...]] = []
    for monomial in leading:
        exponents = tuple(int(e) for e in monomial)
        if len(exponents) != nvars:
            raise ValueError(
                f"leading monomial {exponents} has {len(exponents)} exponents, "
                f"expected {nvars}"
            )
        if any(e < 0 for e in exponents):
            raise ValueError(f"leading monomial {exponents} has a negative exponent")
        monomials.append(exponents)
    return tuple(monomials)


def _supports(leading: Iterable[Sequence[int]], nvars: int) -> tuple[int, ...]:
    """One bitmask per leading monomial, marking which variables it uses."""
    masks = set()
    for monomial in _normalize(leading, nvars):
        mask = 0
        for index, exponent in enumerate(monomial):
            if exponent:
                mask |= 1 << index
        masks.add(mask)
    return tuple(masks)


def _divides(monomial: Sequence[int], candidate: Sequence[int]) -> bool:
    return all(a <= b for a, b in zip(monomial, candidate))
