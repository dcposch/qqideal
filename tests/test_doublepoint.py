"""Double-point ideals of plane parametrizations.

Each test asks one question: does `t -> (p(t), q(t))` send two *distinct*
parameters to the same point? That is exactly what the saturated double-point
ideal records, and nothing more -- see the note on the cusp below.
"""

from __future__ import annotations

import shutil

import pytest

from qqideal import Ideal, Kind, MsolveInputError, Ring, RingMismatch, double_point_ideal

needs_msolve = pytest.mark.skipif(
    shutil.which("msolve") is None, reason="msolve binary not on PATH"
)

T = Ring("t")


def test_construction_is_symbolic() -> None:
    ideal = double_point_ideal("t", "t^2", ring=T)
    assert isinstance(ideal, Ideal)
    assert ideal.ring.names == ("s", "t", "u")
    # p(s) - p(t) and q(s) - q(t), then the slack relation for s - t.
    assert ideal.gens == (
        ideal.ring("s-t"),
        ideal.ring("s^2-t^2"),
        ideal.ring("s*u-t*u-1"),
    )


def test_polys_may_replace_strings() -> None:
    assert double_point_ideal(T("t^2"), T("t^3")).gens == double_point_ideal(
        "t^2", "t^3", ring=T
    ).gens


def test_source_ring_must_be_univariate() -> None:
    R = Ring("x", "y")
    with pytest.raises(MsolveInputError):
        double_point_ideal("x", "y", ring=R)


def test_p_and_q_must_share_a_ring() -> None:
    with pytest.raises(RingMismatch):
        double_point_ideal(T("t"), Ring("s")("s"))


def test_a_ring_is_required_for_strings() -> None:
    with pytest.raises(MsolveInputError):
        double_point_ideal("t", "t^2")


@needs_msolve
def test_smooth_embedding_has_no_double_point() -> None:
    # t -> (t, t^2) is injective, so after removing the diagonal nothing is left.
    assert double_point_ideal("t", "t^2", ring=T).contains_one().kind is Kind.EMPTY


@needs_msolve
def test_the_cusp_has_no_double_point_either() -> None:
    # t -> (t^2, t^3) is injective: s^2 = t^2 forces s = -t, and s^3 = t^3 then
    # forces t = 0, i.e. s = t, which the saturation has already removed. The
    # cusp is not an embedding, but a double-point ideal is not what shows that
    # -- the derivative vanishing at 0 is. This function builds the ideal and
    # makes no claim about embeddings.
    assert double_point_ideal("t^2", "t^3", ring=T).contains_one().kind is Kind.EMPTY


@needs_msolve
def test_the_node_has_a_double_point() -> None:
    # t -> (t^2-1, t^3-t) sends 1 and -1 to the origin: two ordered pairs.
    verdict = double_point_ideal("t^2-1", "t^3-t", ring=T).verdict()
    assert verdict.kind is Kind.NONEMPTY
    assert (verdict.dim, verdict.degree) == (0, 2)


@needs_msolve
def test_even_exponents_give_a_curve_of_double_points() -> None:
    # t -> (t^8, t^6) factors through t -> t^2, so every t != 0 shares its image
    # with -t. The double-point locus is the punctured line s + t = 0.
    verdict = double_point_ideal("t^8", "t^6", ring=T).verdict()
    assert verdict.kind is Kind.NONEMPTY
    assert verdict.dim == 1


@needs_msolve
def test_the_parameter_names_are_configurable() -> None:
    ideal = double_point_ideal("t^2-1", "t^3-t", ring=T, names=("a", "b"))
    assert ideal.ring.names == ("a", "b", "u")
    assert ideal.verdict().degree == 2
