"""Verdicts: the shape of an answer, and the refusal to be a boolean."""

from __future__ import annotations

import dataclasses
import shutil

import pytest

from qqideal import Certainty, Ideal, Kind, Ring, Verdict, ideal_verdict

needs_msolve = pytest.mark.skipif(
    shutil.which("msolve") is None, reason="msolve binary not on PATH"
)

EMPTY = Verdict(
    kind=Kind.EMPTY, dim=None, degree=None, certainty=Certainty.MODULAR
)


def test_verdict_has_no_truth_value() -> None:
    with pytest.raises(TypeError):
        bool(EMPTY)
    with pytest.raises(TypeError):
        if EMPTY:  # noqa: SIM103 - the point of the test
            pass
    with pytest.raises(TypeError):
        not EMPTY


def test_str_always_names_kind_and_certainty() -> None:
    for verdict in (
        EMPTY,
        Verdict(Kind.NONEMPTY, 0, 4, Certainty.PROVEN),
        Verdict(Kind.TIMEOUT, None, None, Certainty.MODULAR, detail="killed"),
    ):
        text = str(verdict)
        assert f"kind={verdict.kind.value}" in text
        assert f"certainty={verdict.certainty.value}" in text
    assert "dim=0" in str(Verdict(Kind.NONEMPTY, 0, 4, Certainty.PROVEN))
    assert "degree=4" in str(Verdict(Kind.NONEMPTY, 0, 4, Certainty.PROVEN))


def test_verdict_is_frozen() -> None:
    with pytest.raises(dataclasses.FrozenInstanceError):
        EMPTY.kind = Kind.NONEMPTY  # type: ignore[misc]


@needs_msolve
def test_using_a_verdict_as_a_condition_is_a_type_error() -> None:
    with pytest.raises(TypeError):
        if ideal_verdict(["x-1", "x"], ring=Ring("x")):
            pass


@needs_msolve
def test_ideal_verdict_accepts_strings_polys_and_ideals() -> None:
    R = Ring("x")
    expected = (Kind.NONEMPTY, 0, 2)
    for gens in (["x^2-1"], [R("x^2-1")], Ideal(["x^2-1"], ring=R)):
        verdict = ideal_verdict(gens, ring=R)
        assert (verdict.kind, verdict.dim, verdict.degree) == expected


@needs_msolve
def test_ideal_verdict_reports_the_msolve_version() -> None:
    verdict = ideal_verdict(["x^2-1"], ring=Ring("x"))
    assert verdict.msolve_version is not None
    assert verdict.msolve_version.startswith("0.10.")


@needs_msolve
def test_opens_saturate_before_the_test() -> None:
    R = Ring("x", "y")
    # V(x*y, x) is the line x = 0, so requiring x != 0 empties it.
    assert ideal_verdict(["x*y", "x"], ring=R, opens=["x"]).kind is Kind.EMPTY
    # Without the open set it is a line.
    assert ideal_verdict(["x*y", "x"], ring=R).dim == 1
    # V(x*y) with x != 0 is still the line y = 0.
    assert ideal_verdict(["x*y"], ring=R, opens=["x"]).dim == 1


@needs_msolve
def test_several_opens_multiply() -> None:
    R = Ring("x", "y")
    # V(x*y) minus both axes is empty; every point of it has x = 0 or y = 0.
    verdict = ideal_verdict(["x*y"], ring=R, opens=["x", R("y")])
    assert verdict.kind is Kind.EMPTY
