"""Witness extraction. Tests that shell out to msolve skip when it is missing.

Every point asserted here is a known solution of a hand-checkable system, so
these tests pin the whole chain: msolve's -P conventions (through msolveio),
the factorization of the eliminating polynomial, the number-field coordinate
arithmetic, the chart un-permutation, and the substitution proof.
"""

from __future__ import annotations

import shutil
from fractions import Fraction

import pytest

from qqideal import (
    AlgebraicPoint,
    Certainty,
    Ideal,
    RationalPoint,
    Ring,
    WitnessKind,
    Witnesses,
    witness_points,
)
from qqideal.ring import coordinate_mod, factor_univariate

needs_msolve = pytest.mark.skipif(
    shutil.which("msolve") is None, reason="msolve binary not on PATH"
)

F = Fraction


# -- the flint facade, no msolve needed --------------------------------------


def test_factor_univariate_is_monic_and_sorted() -> None:
    # 2*(t^2 - 3)*(2t + 1) = 4t^3 + 2t^2 - 12t - 6
    factors = factor_univariate((-6, -12, 2, 4))
    assert factors == (
        ((F(1, 2), F(1)), 1),
        ((F(-3), F(0), F(1)), 1),
    )


def test_factor_univariate_reports_multiplicity() -> None:
    # (t - 1)^2
    assert factor_univariate((1, -2, 1)) == (((F(-1), F(1)), 2),)


def test_coordinate_mod_inverts_wprime() -> None:
    # w = t^2 - 3, w' = 2t, v = -4t, scale 1: coordinate is 4t/(2t) = 2.
    assert coordinate_mod((0, -4), 1, (0, 2), (-3, 0, 1)) == (F(2), F(0))
    # x = 1/2 as printed by msolve: v = -3, scale 2 against w = 3t - 1.
    assert coordinate_mod((-3,), 2, (3,), (F(-1, 3), F(1))) == (F(1, 2),)


def test_coordinate_mod_rejects_repeated_factor() -> None:
    # g = t - 1 divides w' of w = (t-1)^2, so w' is not invertible mod g.
    with pytest.raises(ValueError):
        coordinate_mod((1,), 1, (-2, 2), (-1, 1))


# -- rational and number-field points ----------------------------------------


@needs_msolve
def test_single_rational_point() -> None:
    R = Ring("x", "y")
    w = Ideal(["2*x-1", "3*y-1"], ring=R).witness_points()
    assert w.kind is WitnessKind.POINTS
    assert w.points == (RationalPoint(coordinates=(F(1, 2), F(1, 3))),)
    assert w.variables == ("x", "y")
    assert w.quotient_degree == 1
    assert w.completeness is Certainty.MODULAR
    assert w.run is not None and w.run.msolve_version.startswith("0.10.")
    assert len(w.run.input_sha256) == 64


@needs_msolve
def test_quadratic_field_point() -> None:
    R = Ring("x", "y")
    w = Ideal(["x-2", "y^2-3"], ring=R).witness_points()
    assert w.kind is WitnessKind.POINTS
    (point,) = w.points
    assert isinstance(point, AlgebraicPoint)
    assert point.min_poly_ascending == (F(-3), F(0), F(1))  # t^2 - 3
    assert point.degree == 2
    # x = 2, y = t, each padded to the field degree.
    assert point.coordinates == ((F(2), F(0)), (F(0), F(1)))
    assert w.quotient_degree == 2


@needs_msolve
def test_mixed_rational_and_algebraic_orbits() -> None:
    R = Ring("x", "y")
    # y = x, x^4 - 4x^2 + 3 = (x^2-1)(x^2-3): points (+-1, +-1), (+-sqrt3, +-sqrt3).
    w = Ideal(["y-x", "x^4-4*x^2+3"], ring=R).witness_points()
    assert w.kind is WitnessKind.POINTS
    assert w.quotient_degree == 4
    rational = [p for p in w.points if isinstance(p, RationalPoint)]
    algebraic = [p for p in w.points if isinstance(p, AlgebraicPoint)]
    assert [p.coordinates for p in rational] == [(F(-1), F(-1)), (F(1), F(1))]
    (orbit,) = algebraic
    assert orbit.min_poly_ascending == (F(-3), F(0), F(1))
    assert orbit.coordinates == ((F(0), F(1)), (F(0), F(1)))  # x = t, y = t


@needs_msolve
def test_points_survive_msolve_variable_reorder() -> None:
    R = Ring("x", "y")
    # msolve reorders to make x the separating variable; coordinates must
    # still come back in (x, y) order.
    w = Ideal(["y", "x^2-2"], ring=R).witness_points()
    (point,) = w.points
    assert isinstance(point, AlgebraicPoint)
    assert point.min_poly_ascending == (F(-2), F(0), F(1))
    assert point.coordinates == ((F(0), F(1)), (F(0), F(0)))  # x = t, y = 0


@needs_msolve
def test_points_survive_added_linear_form() -> None:
    R = Ring("x", "y")
    # {(0,0), (0,1), (1,1)}: neither variable separates, msolve adds one.
    w = Ideal(["x^2-x", "y^2-y", "x*y-x"], ring=R).witness_points()
    assert w.kind is WitnessKind.POINTS
    assert [p.coordinates for p in w.points] == [
        (F(0), F(0)),
        (F(0), F(1)),
        (F(1), F(1)),
    ]


@needs_msolve
def test_non_radical_ideal_yields_distinct_points() -> None:
    R = Ring("x", "y")
    # x^2, y-1: one point of multiplicity two.
    w = Ideal(["x^2", "y-1"], ring=R).witness_points()
    assert w.quotient_degree == 2
    assert w.points == (RationalPoint(coordinates=(F(0), F(1))),)


# -- non-point outcomes ------------------------------------------------------


@needs_msolve
def test_empty_locus() -> None:
    R = Ring("x", "y")
    w = Ideal(["x", "y", "x+1"], ring=R).witness_points()
    assert w.kind is WitnessKind.EMPTY
    assert w.points is None
    assert w.completeness is Certainty.MODULAR


@needs_msolve
def test_positive_dimensional_locus() -> None:
    R = Ring("x", "y")
    w = Ideal(["x-1"], ring=R).witness_points()
    assert w.kind is WitnessKind.POSITIVE_DIMENSIONAL
    assert w.points is None


def test_zero_ideal_short_circuits() -> None:
    R = Ring("x", "y")
    w = Ideal([R.constant(0)], ring=R).witness_points()
    assert w.kind is WitnessKind.POSITIVE_DIMENSIONAL
    assert w.run is None
    assert w.completeness is Certainty.PROVEN


def test_prime_field_raises() -> None:
    R = Ring("x", characteristic=65537)
    with pytest.raises(NotImplementedError):
        Ideal(["x-1"], ring=R).witness_points()


def test_witnesses_has_no_truth_value() -> None:
    w = Witnesses(
        kind=WitnessKind.EMPTY,
        points=None,
        variables=("x",),
        quotient_degree=None,
        completeness=Certainty.MODULAR,
    )
    with pytest.raises(TypeError):
        bool(w)


# -- opens= ------------------------------------------------------------------


@needs_msolve
def test_opens_projects_slack_away() -> None:
    R = Ring("x", "y")
    # V(x^2-x, y) = {(0,0), (1,0)}; away from V(x) only (1,0) survives.
    w = witness_points(["x^2-x", "y"], ring=R, opens=["x"])
    assert w.kind is WitnessKind.POINTS
    assert w.variables == ("x", "y")
    assert w.points == (RationalPoint(coordinates=(F(1), F(0))),)


@needs_msolve
def test_opens_with_algebraic_point() -> None:
    R = Ring("x", "y")
    # V(y-x, x^3-3x) = {(0,0)} and (+-sqrt3, +-sqrt3); drop the origin.
    w = witness_points(["y-x", "x^3-3*x"], ring=R, opens=["x"])
    assert w.kind is WitnessKind.POINTS
    (orbit,) = w.points
    assert isinstance(orbit, AlgebraicPoint)
    assert orbit.min_poly_ascending[-1] == F(1)  # monic
    assert orbit.degree == 2
    # y = x on the orbit, and the orbit is not the origin. Which element of
    # Q(sqrt 3) presents x depends on msolve's separating choice for the
    # saturated system, so assert the relations rather than the presentation;
    # membership itself was substitution-proven during extraction.
    x_coord, y_coord = orbit.coordinates
    assert x_coord == y_coord
    assert any(c != 0 for c in x_coord)


@needs_msolve
def test_opens_can_empty_the_locus() -> None:
    R = Ring("x", "y")
    w = witness_points(["x", "y"], ring=R, opens=["x"])
    assert w.kind is WitnessKind.EMPTY
    assert w.variables == ("x", "y")


# -- determinism -------------------------------------------------------------


@needs_msolve
def test_extraction_is_deterministic() -> None:
    R = Ring("x", "y")
    gens = ["y-x", "x^4-4*x^2+3"]
    first = Ideal(gens, ring=R).witness_points()
    second = Ideal(gens, ring=R).witness_points()
    assert first.points == second.points
    assert first.run.output_sha256 == second.run.output_sha256
