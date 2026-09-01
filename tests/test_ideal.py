"""Ideals. Tests that shell out to msolve skip when the binary is missing."""

from __future__ import annotations

import shutil

import pytest

from qqideal import Certainty, Ideal, Kind, MsolveInputError, Ring, RingMismatch
from qqideal.dimdeg import degree, dim_and_degree, dimension
from qqideal.errors import MsolveDied, MsolveTimeout

needs_msolve = pytest.mark.skipif(
    shutil.which("msolve") is None, reason="msolve binary not on PATH"
)


def cyclic(ring: Ring) -> list[str]:
    """The cyclic-n system: small to write down, expensive to solve."""
    names = ring.names
    n = len(names)
    system = [
        "+".join("*".join(names[(i + j) % n] for j in range(d)) for i in range(n))
        for d in range(1, n)
    ]
    return [*system, "*".join(names) + "-1"]


# -- construction ------------------------------------------------------------


def test_generators_may_be_strings_or_polys() -> None:
    R = Ring("x", "y")
    assert Ideal(["x*y-1", R("y^2-x")], ring=R).gens == (R("x*y-1"), R("y^2-x"))


def test_ring_is_inferred_from_a_poly() -> None:
    R = Ring("x", "y")
    assert Ideal([R("x")]).ring == R


def test_ring_is_required_for_strings() -> None:
    with pytest.raises(MsolveInputError):
        Ideal(["x"])


def test_generators_must_share_the_ring() -> None:
    R, S = Ring("x", "y"), Ring("x")
    with pytest.raises(RingMismatch):
        Ideal([R("x"), S("x")], ring=R)


def test_a_bare_string_is_not_a_generator_list() -> None:
    with pytest.raises(MsolveInputError):
        Ideal("x*y-1", ring=Ring("x", "y"))  # type: ignore[arg-type]


def test_to_msolve_is_a_canonical_file() -> None:
    R = Ring("x", "y")
    assert Ideal(["x^2-y", "x*y-1"], ring=R).to_msolve() == "x,y\n0\nx^2-y,\nx*y-1\n"
    F = Ring("x", characteristic=65521)
    assert Ideal(["x^2-1"], ring=F).to_msolve() == "x\n65521\nx^2+65520\n"


# -- dimension and degree, straight from exponent vectors --------------------


def test_dimension_is_combinatorial() -> None:
    assert dimension([(2, 0)], 2) == 1  # V(x^2) is a line
    assert dimension([(1, 1)], 2) == 1  # V(x*y) is two lines
    assert dimension([(2, 0), (0, 1)], 2) == 0
    assert dimension([], 3) == 3  # the zero ideal
    assert dimension([(0, 0)], 2) == -1  # the unit ideal


def test_degree_counts_standard_monomials() -> None:
    assert degree([(2,)], 1) == 2
    assert degree([(2, 0), (0, 1)], 2) == 2
    assert degree([(3, 0, 0), (0, 1, 0), (0, 0, 1)], 3) == 3
    assert degree([(2, 0), (1, 1), (0, 2)], 2) == 3
    assert dim_and_degree([(2, 0)], 2) == (1, None)


def test_degree_needs_dimension_zero() -> None:
    with pytest.raises(ValueError):
        degree([(2, 0)], 2)


def test_leading_monomials_must_match_the_variable_count() -> None:
    with pytest.raises(ValueError):
        dimension([(1, 0, 0)], 2)


# -- verdicts, without msolve ------------------------------------------------


def test_zero_ideal_needs_no_solver() -> None:
    R = Ring("x", "y")
    verdict = Ideal(["0"], ring=R).verdict()
    assert verdict.kind is Kind.NONEMPTY
    assert (verdict.dim, verdict.degree) == (2, None)
    assert verdict.certainty is Certainty.PROVEN
    assert Ideal([], ring=R).groebner() == ()
    assert Ideal([], ring=R).leading_monomials() == ()


def test_timeout_becomes_a_verdict(monkeypatch: pytest.MonkeyPatch) -> None:
    def boom(*args: object, **kwargs: object) -> None:
        raise MsolveTimeout("nope", timeout=0.5)

    monkeypatch.setattr("qqideal.ideal.run_groebner", boom)
    verdict = Ideal(["x"], ring=Ring("x")).verdict(timeout=0.5)
    assert verdict.kind is Kind.TIMEOUT
    assert (verdict.dim, verdict.degree) == (None, None)
    assert verdict.detail == "nope"


def test_a_dead_solver_becomes_a_verdict(monkeypatch: pytest.MonkeyPatch) -> None:
    def boom(*args: object, **kwargs: object) -> None:
        raise MsolveDied("segfault")

    monkeypatch.setattr("qqideal.ideal.run_groebner", boom)
    assert Ideal(["x"], ring=Ring("x")).verdict().kind is Kind.ERROR


def test_radical_is_not_implemented() -> None:
    with pytest.raises(NotImplementedError):
        Ideal(["x^2"], ring=Ring("x")).radical()


def test_saturate_adjoins_one_slack_variable() -> None:
    R = Ring("x", "y")
    saturated = Ideal(["x*y"], ring=R).saturate("x")
    assert saturated.ring.names == ("x", "y", "u")
    assert saturated.gens[-1] == saturated.ring("x*u-1")
    assert Ideal(["x*y"], ring=R).colon("x").gens == saturated.gens
    # A second saturation does not collide with the first slack variable.
    assert saturated.saturate("y").ring.names == ("x", "y", "u", "u1")


# -- verdicts, with msolve ---------------------------------------------------


@needs_msolve
def test_unit_ideal_is_empty() -> None:
    verdict = Ideal(["x-1", "x"], ring=Ring("x", "y")).contains_one()
    assert verdict.kind is Kind.EMPTY
    assert (verdict.dim, verdict.degree) == (None, None)
    # msolve 0.10.1 returns after its first modular prime for the unit ideal.
    assert verdict.certainty is Certainty.MODULAR


@needs_msolve
def test_zero_dimensional_ideal_has_a_degree() -> None:
    verdict = Ideal(["x^2-1"], ring=Ring("x")).verdict()
    assert verdict.kind is Kind.NONEMPTY
    assert (verdict.dim, verdict.degree) == (0, 2)
    assert verdict.certainty is Certainty.PROVEN


@needs_msolve
def test_positive_dimensional_ideal_has_no_degree() -> None:
    verdict = Ideal(["x"], ring=Ring("x", "y")).verdict()
    assert verdict.kind is Kind.NONEMPTY
    assert (verdict.dim, verdict.degree) == (1, None)


@needs_msolve
def test_degree_counts_multiplicity() -> None:
    verdict = Ideal(["x^2", "y"], ring=Ring("x", "y")).verdict()
    assert (verdict.dim, verdict.degree) == (0, 2)


@needs_msolve
def test_groebner_basis_comes_back_as_polynomials() -> None:
    R = Ring("x", "y")
    basis = Ideal(["2*x^2-1/3", "x*y-1"], ring=R).groebner()
    assert set(basis) == {R("6*x-y"), R("y^2-6")}
    assert Ideal(["x-1", "x"], ring=R).groebner() == (R("1"),)


@needs_msolve
def test_leading_monomials() -> None:
    R = Ring("x", "y")
    assert Ideal(["x*y-1", "y^2-x"], ring=R).leading_monomials() == (
        (0, 2),
        (1, 1),
        (2, 0),
    )
    assert Ideal(["x-1", "x"], ring=R).leading_monomials() == ((0, 0),)


@needs_msolve
def test_prime_field_verdicts_are_proven() -> None:
    F = Ring("x", characteristic=65521)
    verdict = Ideal(["x^2-1"], ring=F).verdict()
    assert (verdict.kind, verdict.dim, verdict.degree) == (Kind.NONEMPTY, 0, 2)
    assert verdict.certainty is Certainty.PROVEN
    empty = Ideal(["x-1", "x"], ring=F).contains_one()
    assert empty.kind is Kind.EMPTY
    assert empty.certainty is Certainty.PROVEN


@needs_msolve
def test_prime_field_groebner_basis() -> None:
    F = Ring("x", characteristic=65521)
    assert Ideal(["x^2-1"], ring=F).groebner() == (F("x^2-1"),)


@needs_msolve
def test_saturation_removes_a_component() -> None:
    R = Ring("x", "y")
    # V(x*y, x) is the line x = 0; away from x = 0 there is nothing left.
    assert Ideal(["x*y", "x"], ring=R).saturate("x").contains_one().kind is Kind.EMPTY
    # V(x*y) away from x = 0 is the line y = 0, which is still there.
    verdict = Ideal(["x*y"], ring=R).saturate("x").verdict()
    assert (verdict.kind, verdict.dim) == (Kind.NONEMPTY, 1)


@needs_msolve
def test_saturating_at_zero_empties_the_ideal() -> None:
    R = Ring("x")
    assert Ideal(["x"], ring=R).saturate("0").contains_one().kind is Kind.EMPTY


@needs_msolve
def test_radical_membership() -> None:
    R = Ring("x", "y")
    # x^2 = 0 forces x = 0, so x is in the radical: nothing survives x != 0.
    assert Ideal(["x^2"], ring=R).radical_member("x").kind is Kind.EMPTY
    # y does not vanish on V(x), so it is not in the radical.
    assert Ideal(["x"], ring=R).radical_member("y").kind is Kind.NONEMPTY


@needs_msolve
def test_msolve_timeout_is_reported_not_raised() -> None:
    ring = Ring(*[f"x{i}" for i in range(6)])
    verdict = Ideal(cyclic(ring), ring=ring).verdict(timeout=0.001)
    assert verdict.kind is Kind.TIMEOUT
    assert verdict.detail is not None


@needs_msolve
def test_groebner_raises_where_verdict_reports() -> None:
    ring = Ring(*[f"x{i}" for i in range(6)])
    with pytest.raises(MsolveTimeout):
        Ideal(cyclic(ring), ring=ring).groebner(timeout=0.001)
