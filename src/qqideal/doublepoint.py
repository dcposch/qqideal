"""The double-point ideal of a polynomial parametrization.

For a plane curve parametrized by ``t -> (p(t), q(t))``, two distinct parameters
land on the same point exactly when ``p(s) = p(t)`` and ``q(s) = q(t)`` with
``s != t``. That is the ideal

    I_DP = (p(s) - p(t), q(s) - q(t)) : (s - t)^oo

in ``QQ[s, t]``: saturating at ``s - t`` throws away the diagonal, which solves
the equations trivially and says nothing.

This module builds that ideal and stops there. It is not a test for whether a
parametrization is an embedding: an injective map can still fail to be one, and
the cusp ``t -> (t^2, t^3)`` is the standard example -- it has no double point
at all, and fails on the derivative instead. Deciding embedding needs the rest
of the postcheck (length, immersivity, distinct tangents), which is not here.
"""

from __future__ import annotations

from .errors import MsolveInputError, RingMismatch
from .ideal import Ideal
from .ring import Poly, Ring

__all__ = ["double_point_ideal"]


def double_point_ideal(
    p: "Poly | str",
    q: "Poly | str",
    *,
    ring: Ring | None = None,
    names: tuple[str, str] = ("s", "t"),
) -> Ideal:
    """Build ``(p(s)-p(t), q(s)-q(t)) : (s-t)^oo`` for the parametrization
    ``(p, q)``.

    :param p: first coordinate, a univariate polynomial or a string.
    :param q: second coordinate, likewise.
    :param ring: the univariate source ring, required if ``p`` or ``q`` is a
        string.
    :param names: the two parameter names to use.
    :returns: an :class:`~qqideal.Ideal`. As with
        :meth:`~qqideal.Ideal.saturate`, the saturation is carried as a
        Rabinowitsch ideal, so the ideal lives in ``QQ[s, t, u]`` with
        ``u*(s-t) - 1`` among its generators, and its variety is the set of
        ordered pairs ``s != t`` with the same image. Emptiness, dimension and
        degree are the ones of the double-point locus.
    :raises MsolveInputError: if the source ring is not univariate.
    """
    source = _source_ring(p, q, ring)
    if source.nvars != 1:
        raise MsolveInputError(
            f"a parametrization is univariate; got {source} with "
            f"{source.nvars} variables"
        )
    if len(names) != 2 or names[0] == names[1]:
        raise MsolveInputError(f"names must be two distinct variables, got {names!r}")

    poly_p = _coerce(source, p)
    poly_q = _coerce(source, q)

    plane = Ring(*names, characteristic=source.characteristic)
    s, t = plane.gens()
    equations = [
        _substitute(poly_p, plane, 0) - _substitute(poly_p, plane, 1),
        _substitute(poly_q, plane, 0) - _substitute(poly_q, plane, 1),
    ]
    return Ideal(equations, ring=plane).saturate(s - t)


def _source_ring(p: "Poly | str", q: "Poly | str", ring: Ring | None) -> Ring:
    rings = {value.ring for value in (p, q) if isinstance(value, Poly)}
    if ring is not None:
        rings.add(ring)
    if not rings:
        raise MsolveInputError("ring= is required when p and q are strings")
    if len(rings) > 1:
        raise RingMismatch(
            f"p and q must live in one ring, got {sorted(str(r) for r in rings)}"
        )
    return rings.pop()


def _coerce(ring: Ring, value: "Poly | str") -> Poly:
    return ring(value) if isinstance(value, str) else ring._check(value)


def _substitute(poly: Poly, plane: Ring, position: int) -> Poly:
    """Send the source variable to variable ``position`` of the plane ring."""
    terms = {}
    for exponents, coeff in poly.terms():
        key = [0, 0]
        key[position] = exponents[0]
        terms[tuple(key)] = coeff
    return plane._wrap(plane._ctx.from_dict(terms))
