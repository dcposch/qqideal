# qqideal

Exact ideals over QQ: python-flint for arithmetic, msolve for Gröbner bases, verdicts that
refuse to guess. Emptiness, dimension, and 0-dimensional degree from grevlex leading ideals;
exact witness points — rational and number-field — from msolve's rational parametrization.
No solver-mode I/O.

## Install

```
pip install qqideal
```

That pulls [msolveio](https://pypi.org/project/msolveio/) and python-flint from PyPI. You also
need a system `msolve` 0.10.x binary on `PATH`.

## Usage

```python
from qqideal import Ideal, Kind, Ring, ideal_verdict

R = Ring("x", "y")
I = Ideal(["x^2-1", "y-x"], ring=R)

verdict = ideal_verdict(I)
print(verdict.kind)       # Kind.NONEMPTY
print(verdict.dim)        # 0
print(verdict.degree)     # 2
print(verdict.certainty)  # Certainty.PROVEN

if verdict.kind is Kind.NONEMPTY:   # never `if verdict:` -- that raises
    print(I.groebner())             # (Poly('x - y', ...), Poly('y^2 - 1', ...))
```

`Verdict.__bool__` raises `TypeError`. `TIMEOUT` and `ERROR` are not answers, and a
truthiness test would silently fold them into one of the two that are.

`opens=` saturates before the test, so you can ask about the complement of a hypersurface:

```python
ideal_verdict(["x*y", "x"], ring=R, opens=["x"]).kind   # Kind.EMPTY
Ideal(["x^2"], ring=R).radical_member("x").kind         # Kind.EMPTY: x is in the radical
```

## Certainty

A unit ideal over Q from msolve `-g` is `Certainty.MODULAR`, not `PROVEN`: msolve 0.10.1
returns after its first modular prime and still prints characteristic 0. A nonempty result
over Q uses a lifted `-g 2` basis and is `Certainty.PROVEN`, as is any result over a prime
field.

## Saturation

`I.saturate(f)` is Rabinowitsch: it returns `I + (u*f - 1)` in the ring extended by one slack
variable. Its variety is `V(I) \ V(f)`, so emptiness, dimension, and degree are those of
`I : f^∞` -- but its generators are not that ideal written back in `R`, which would need an
elimination order msolve's Gröbner mode does not offer. `colon` is an alias, and in v0.1 the
colon is the saturation.

## Witness points

`Ideal.witness_points()` runs msolve's rational parametrization (`-P`, through msolveio's
strict parser) and factors the eliminating polynomial over Q. Each irreducible factor is one
Galois orbit of solutions:

```python
from fractions import Fraction
from qqideal import Ideal, Ring, WitnessKind, witness_points

R = Ring("x", "y")
w = Ideal(["x-2", "y^2-3"], ring=R).witness_points()

w.kind                # WitnessKind.POINTS
(point,) = w.points   # one orbit: (2, ±√3)
point.min_poly_ascending   # (Fraction(-3), Fraction(0), Fraction(1))   t² - 3
point.coordinates          # ((Fraction(2), Fraction(0)), (Fraction(0), Fraction(1)))
                           # x = 2, y = t, as vectors in the basis 1, t
```

Degree-1 factors come back as `RationalPoint`s with plain `Fraction` coordinates. Everything
is exact; nothing is a float.

Two claims come out of this with different strengths, and they are graded separately.
**Membership is proven here**: every returned point is substituted back into every generator
with exact arithmetic before it is returned, and a nonzero value raises instead of
returning. **Completeness is msolve's claim**: that these are *all* the solutions rests on a
rational lift msolve 0.10.1 does not certify, so `Witnesses.completeness` is
`Certainty.MODULAR`. `Witnesses` has no truth value, like `Verdict`; branch on `.kind`,
where `TIMEOUT` and `ERROR` stay distinct from `EMPTY` and `POSITIVE_DIMENSIONAL`.

`opens=` works the way it does for `ideal_verdict` and hands back points in the original
chart, with the Rabinowitsch slack coordinate projected away:

```python
witness_points(["x^2-x", "y"], ring=R, opens=["x"]).points
# (RationalPoint(coordinates=(Fraction(1), Fraction(0))),)   — (0,0) fell to the open
```

msolve parametrizes the radical, so `quotient_degree` counts solutions with multiplicity
while the points are the distinct ones. The full msolveio `ParamResult` — argv, versions,
input/output SHA-256 — rides along as `Witnesses.run` for custody.

## Not in v0.2

Primary decomposition, positive-dimensional radicals, radical computation of any kind
(membership only), witness points over prime fields, real-root boxes, Macaulay2. These raise
`NotImplementedError` rather than returning an approximation.

## License

MIT © 2026 DC Posch — <https://github.com/dcposch/qqideal>
