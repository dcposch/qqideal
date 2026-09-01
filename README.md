# qqideal

Exact ideals over QQ: python-flint for arithmetic, msolve for Gröbner bases, verdicts that
refuse to guess. Emptiness, dimension, and 0-dimensional degree from grevlex leading ideals.
No solver-mode I/O.

## Install

```
pip install git+https://github.com/dcposch/qqideal.git
```

That pulls [msolveio](https://github.com/dcposch/msolveio) from GitHub and python-flint from
PyPI. You also need a system `msolve` 0.10.x binary on `PATH`.

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
over Q, and any result over a prime field, is `Certainty.PROVEN`.

## Saturation

`I.saturate(f)` is Rabinowitsch: it returns `I + (u*f - 1)` in the ring extended by one slack
variable. Its variety is `V(I) \ V(f)`, so emptiness, dimension, and degree are those of
`I : f^∞` -- but its generators are not that ideal written back in `R`, which would need an
elimination order msolve's Gröbner mode does not offer. `colon` is an alias, and in v0.1 the
colon is the saturation.

## Not in v0.1

Primary decomposition, positive-dimensional radicals, radical computation of any kind
(membership only), solver mode / `-P`, Macaulay2. These raise `NotImplementedError` rather
than returning an approximation.

## License

MIT © 2026 DC Posch — <https://github.com/dcposch/qqideal>
