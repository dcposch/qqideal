"""Ring and Poly: parsing, arithmetic, resultants. None of this needs msolve."""

from __future__ import annotations

from fractions import Fraction

import pytest

from qqideal import MsolveInputError, Poly, Ring, RingMismatch


def test_parse_expanded_monomial_sum() -> None:
    R = Ring("x", "y")
    f = R("x^2-y")
    assert f.degree() == 2
    assert f.terms() == (((2, 0), 1), ((0, 1), -1))
    assert f == R("x^2 - y")


def test_parse_leading_rational_coefficient() -> None:
    R = Ring("x")
    assert R("1/2*x") * 2 == R("x")


@pytest.mark.parametrize(
    "source",
    ["(x+1)", "x/2", "x_1", "z", "x^", "2x", "x + x", "x**2", ""],
)
def test_rejected_sources(source: str) -> None:
    R = Ring("x", "y")
    with pytest.raises(MsolveInputError):
        R(source)


def test_parse_rejects_non_strings() -> None:
    R = Ring("x")
    with pytest.raises(MsolveInputError):
        R(1)  # type: ignore[arg-type]


@pytest.mark.parametrize("names", [("x_1",), ("x", "x"), ("2x",), ("x y",), ()])
def test_bad_variable_names(names: tuple[str, ...]) -> None:
    with pytest.raises(MsolveInputError):
        Ring(*names)


@pytest.mark.parametrize("characteristic", [4, 1, -7, 2**31])
def test_bad_characteristic(characteristic: int) -> None:
    with pytest.raises(MsolveInputError):
        Ring("x", characteristic=characteristic)


def test_arithmetic() -> None:
    R = Ring("x", "y")
    f, g = R("x^2-y"), R("1/2*x")
    assert f + g == R("x^2 + 1/2*x - y")
    assert f - g == R("x^2 - 1/2*x - y")
    assert f * g == R("1/2*x^3 - 1/2*x*y")
    assert -f == R("y - x^2")
    assert f - f == R("0")
    assert (g * 2) ** 2 == R("x^2")
    assert 1 + R("x") == R("x+1")
    assert 1 - R("x") == R("1-x")
    assert R("x") * Fraction(1, 2) == R("1/2*x")


def test_zero_polynomial() -> None:
    R = Ring("x")
    zero = R("0")
    assert zero.is_zero()
    assert zero.degree() == -1
    with pytest.raises(ValueError):
        zero.leading_monomial()


def test_leading_monomial_is_grevlex() -> None:
    R = Ring("x", "y")
    # In grevlex, x^3 > x^2*y > y^3.
    assert R("y^3 + x^2*y + x^3").leading_monomial() == (3, 0)


def test_resultant() -> None:
    R = Ring("x", "y")
    resultant = R.resultant(R("x-y"), R("x-2"), "x")
    assert resultant in (R("y-2"), R("2-y"))


def test_discriminant() -> None:
    R = Ring("x", "y")
    assert R.discriminant(R("x^2-y"), "x") == R("4*y")


def test_resultant_needs_a_declared_variable() -> None:
    R = Ring("x", "y")
    with pytest.raises(MsolveInputError):
        R.resultant(R("x"), R("y"), "z")


def test_resultant_over_prime_field_is_not_implemented() -> None:
    F = Ring("x", "y", characteristic=65521)
    with pytest.raises(NotImplementedError):
        F.resultant(F("x-y"), F("x-2"), "x")


def test_ring_mismatch() -> None:
    R = Ring("x", "y")
    S = Ring("x")
    with pytest.raises(RingMismatch):
        R("x") + S("x")
    with pytest.raises(RingMismatch):
        R("x") + "x"  # type: ignore[operator]


def test_prime_field_coefficients_are_reduced() -> None:
    F = Ring("x", characteristic=65521)
    assert F("65520*x") == -F("x")
    assert (F("x") * 65521).is_zero()
    # 1/2 mod 65521 is 32761, and 2 * 32761 == 65522 == 1.
    assert F("1/2*x") == F("32761*x")
    assert F("1/2*x") * 2 == F("x")


def test_prime_field_rejects_denominator_divisible_by_p() -> None:
    F = Ring("x", characteristic=7)
    with pytest.raises(MsolveInputError):
        F("1/7*x")


def test_to_msolve_round_trips() -> None:
    R = Ring("x", "y")
    for source in ["x^2-y", "1/2*x*y - 3", "-x^3 + 2*x*y^2 - 1/3", "0"]:
        assert R(R(source).to_msolve()) == R(source)


def test_extend_embed_and_fresh_variable() -> None:
    R = Ring("x", "y")
    assert R.fresh_variable("u") == "u"
    assert R.extend("u").fresh_variable("u") == "u1"
    big = R.extend("u")
    assert big.embed(R("x^2-y")) == big("x^2-y")
    with pytest.raises(RingMismatch):
        R.embed(big("u"))
    with pytest.raises(RingMismatch):
        Ring("x", characteristic=7).embed(R("x"))


def test_equality_and_hashing() -> None:
    R = Ring("x", "y")
    assert Ring("x", "y") == R
    assert Ring("y", "x") != R
    assert Ring("x", "y", characteristic=7) != R
    assert hash(R("x^2-y")) == hash(R("x^2 - y"))
    assert len({R("x"), R("x"), R("y")}) == 2
    assert isinstance(R("x"), Poly)


def test_repr_and_str() -> None:
    R = Ring("x", "y")
    assert str(R) == "QQ[x, y]"
    assert str(Ring("x", characteristic=7)) == "GF(7)[x]"
    assert repr(Ring("x", characteristic=7)) == "Ring('x', characteristic=7)"
    assert str(R("x^2-y")) == "x^2 - y"
