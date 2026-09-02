"""Witness points: exact solutions of zero-dimensional ideals over Q.

Built on msolve's rational parametrization (``-P``, via msolveio 0.2's strict
parser). The eliminating polynomial ``w`` is factored over Q; each irreducible
factor is one Galois orbit of solutions, returned as a :class:`RationalPoint`
(degree 1) or an :class:`AlgebraicPoint` whose coordinates live in the number
field ``Q[t]/(g)``. Nothing is a float and nothing is approximated.

Two claims with two different strengths come out of this, and they are graded
separately:

* **Membership is proven here.** Every returned point is substituted back into
  every generator with exact arithmetic before it is returned; a nonzero value
  raises instead of returning. A point you receive satisfies the system -- that
  is a theorem, not msolve's word.
* **Completeness is msolve's claim.** That these are *all* the solutions rests
  on msolve's parametrization, whose rational lift msolve 0.10.1 does not
  certify; :attr:`Witnesses.completeness` is therefore
  :attr:`~qqideal.Certainty.MODULAR`, never
  :attr:`~qqideal.Certainty.PROVEN`.

msolve parametrizes the radical, so multiplicity is not recoverable:
:attr:`Witnesses.quotient_degree` counts solutions with multiplicity while the
points are the distinct ones, and the two disagree exactly when the ideal is
not radical.
"""

from __future__ import annotations

import dataclasses
import enum
from dataclasses import dataclass
from fractions import Fraction
from typing import TYPE_CHECKING, Sequence

from msolveio import (
    EmptySolutionSet,
    ParamResult,
    PositiveDimensional,
    RationalParametrization,
    run_param,
)
from msolveio.errors import MsolveDied, MsolveOutputError, MsolveTimeout

from .ring import Ascending, coordinate_mod, factor_univariate, vanishes_at_algebraic
from .verdict import Certainty

if TYPE_CHECKING:  # pragma: no cover - import cycle broken at runtime
    from .ideal import Ideal
    from .ring import Poly, Ring

__all__ = [
    "WitnessKind",
    "RationalPoint",
    "AlgebraicPoint",
    "Witnesses",
    "witness_points",
]


class WitnessKind(enum.Enum):
    """What witness extraction established."""

    #: The locus is zero-dimensional and its points were extracted.
    POINTS = "points"
    #: The unit ideal: no solutions in an algebraic closure.
    EMPTY = "empty"
    #: Infinitely many solutions; there is no finite witness list.
    POSITIVE_DIMENSIONAL = "positive-dimensional"
    #: msolve exceeded the timeout. Nothing was established.
    TIMEOUT = "timeout"
    #: msolve died or emitted bytes we refuse to interpret. Nothing was
    #: established.
    ERROR = "error"


@dataclass(frozen=True)
class RationalPoint:
    """A solution with rational coordinates, in the variable order of
    :attr:`Witnesses.variables`."""

    coordinates: tuple[Fraction, ...]


@dataclass(frozen=True)
class AlgebraicPoint:
    """One Galois orbit of solutions, with coordinates in ``Q[t]/(g)``.

    ``g`` is monic and irreducible over Q of degree at least 2, and each of the
    ``deg g`` embeddings of ``t`` into an algebraic closure yields one solution
    of the system. The presentation field is exactly the field the coordinates
    generate: msolve's parameter is a separating element, so ``Q(t)`` equals
    ``Q`` adjoined the coordinates -- no larger.

    :param min_poly_ascending: ascending coefficients of ``g``.
    :param coordinates: one entry per variable of
        :attr:`Witnesses.variables`; each is the ascending coefficient vector
        of the coordinate in the basis ``1, t, .., t^(deg g - 1)``, padded to
        length ``deg g``.
    """

    min_poly_ascending: Ascending
    coordinates: tuple[Ascending, ...]

    @property
    def degree(self) -> int:
        """The degree of the presentation field over Q."""
        return len(self.min_poly_ascending) - 1


@dataclass(frozen=True)
class Witnesses:
    """The result of witness extraction. Not a boolean and not a list.

    :param kind: what was established; see :class:`WitnessKind`.
    :param points: the witness points, deterministically sorted, when ``kind``
        is :attr:`WitnessKind.POINTS`; ``None`` otherwise -- including for
        :attr:`WitnessKind.EMPTY`, so a non-answer can never be iterated as an
        empty solution set.
    :param variables: the coordinate order every point uses.
    :param quotient_degree: the number of solutions counted with multiplicity,
        when the locus is zero-dimensional. The number of returned points can
        be smaller: msolve parametrizes the radical.
    :param completeness: how firmly "these are all the solutions" is
        established. Membership of each returned point is proven by exact
        substitution regardless; this field grades only the completeness
        claim, and is :attr:`~qqideal.Certainty.MODULAR` whenever msolve
        produced the parametrization.
    :param run: the full msolveio :class:`~msolveio.ParamResult`, custody
        fields included, when msolve ran to completion.
    :param detail: why msolve timed out or died, for the non-answers.
    """

    kind: WitnessKind
    points: tuple[RationalPoint | AlgebraicPoint, ...] | None
    variables: tuple[str, ...]
    quotient_degree: int | None
    completeness: Certainty
    run: ParamResult | None = None
    detail: str | None = None

    def __bool__(self) -> bool:
        raise TypeError(
            "Witnesses has no truth value; test .kind explicitly, e.g. "
            "`if witnesses.kind is WitnessKind.POINTS:`. TIMEOUT and ERROR "
            "are not answers and must not collapse into one."
        )

    def __str__(self) -> str:
        parts = [f"kind={self.kind.value}", f"completeness={self.completeness.value}"]
        if self.points is not None:
            parts.append(f"points={len(self.points)}")
        if self.quotient_degree is not None:
            parts.append(f"quotient_degree={self.quotient_degree}")
        if self.detail:
            parts.append(f"detail={self.detail}")
        return f"Witnesses({', '.join(parts)})"


def witness_points(
    gens: "Ideal | Sequence[Poly | str]",
    *,
    ring: "Ring | None" = None,
    opens: "Sequence[Poly | str]" = (),
    timeout: float = 60,
) -> Witnesses:
    """Extract exact witness points, optionally away from hypersurfaces.

    :param gens: an :class:`~qqideal.Ideal`, or generators as
        :class:`~qqideal.Poly` or strings together with ``ring=``.
    :param ring: the ring, required when ``gens`` contains strings.
    :param opens: polynomials that must not vanish. The ideal is saturated at
        their product first, so the witnesses are points of
        ``V(I) \\ V(f1*..*fk)``; the slack coordinate of the Rabinowitsch
        trick is projected away and the returned points use the original
        variables.
    :param timeout: wall-clock limit for msolve, in seconds.
    :returns: a :class:`Witnesses`. Never raises for a solver failure; a
        msolve timeout or crash comes back as :attr:`WitnessKind.TIMEOUT` or
        :attr:`WitnessKind.ERROR`. Bad input still raises, and so does a
        point that fails its substitution check -- that is corrupt solver
        output, not a result.
    """
    from .ideal import Ideal

    ideal = gens if isinstance(gens, Ideal) else Ideal(gens, ring=ring)
    original = ideal.ring
    if opens:
        product = original.constant(1)
        for open_set in opens:
            product = product * ideal._coerce(open_set)
        ideal = ideal.saturate(product)

    witnesses = extract(ideal, timeout=timeout)
    if ideal.ring is original:
        return witnesses

    # Project the Rabinowitsch slack coordinate away: the saturated variety is
    # the graph of V(I) \ V(f), so dropping the slack is a bijection onto it.
    keep = tuple(ideal.ring.names.index(name) for name in original.names)
    points = witnesses.points
    if points is not None:
        points = tuple(
            dataclasses.replace(
                point,
                coordinates=tuple(point.coordinates[i] for i in keep),
            )
            for point in points
        )
        points = tuple(sorted(points, key=_point_sort_key))
    return dataclasses.replace(witnesses, points=points, variables=original.names)


def extract(ideal: "Ideal", *, timeout: float = 60) -> Witnesses:
    """Extract witness points of ``ideal`` in its own ring, no saturation.

    This is what :meth:`qqideal.Ideal.witness_points` calls; use the module
    function :func:`witness_points` for ``opens=`` support.

    :raises NotImplementedError: on a prime-field ring. Witness extraction is
        a characteristic-0 activity; msolve's characteristic-p ``-P`` output
        is a different grammar msolveio refuses to parse.
    """
    ring = ideal.ring
    if ring.characteristic:
        raise NotImplementedError(
            f"witness extraction is implemented over Q only; this ring has "
            f"characteristic {ring.characteristic}"
        )
    if ideal.is_zero_ideal:
        return Witnesses(
            kind=WitnessKind.POSITIVE_DIMENSIONAL,
            points=None,
            variables=ring.names,
            quotient_degree=None,
            completeness=Certainty.PROVEN,
            detail="zero ideal; no msolve call needed",
        )

    try:
        result = run_param(ideal.to_msolve(), timeout=timeout)
    except MsolveTimeout as exc:
        return Witnesses(
            kind=WitnessKind.TIMEOUT,
            points=None,
            variables=ring.names,
            quotient_degree=None,
            completeness=Certainty.MODULAR,
            detail=str(exc),
        )
    except (MsolveDied, MsolveOutputError) as exc:
        return Witnesses(
            kind=WitnessKind.ERROR,
            points=None,
            variables=ring.names,
            quotient_degree=None,
            completeness=Certainty.MODULAR,
            detail=str(exc),
        )

    output = result.output
    if isinstance(output, EmptySolutionSet):
        return Witnesses(
            kind=WitnessKind.EMPTY,
            points=None,
            variables=ring.names,
            quotient_degree=None,
            completeness=Certainty.MODULAR,
            run=result,
        )
    if isinstance(output, PositiveDimensional):
        return Witnesses(
            kind=WitnessKind.POSITIVE_DIMENSIONAL,
            points=None,
            variables=ring.names,
            quotient_degree=None,
            completeness=Certainty.MODULAR,
            run=result,
        )

    assert isinstance(output, RationalParametrization)
    assert result.chart is not None
    points = _points_from(output, result.chart, ideal)
    return Witnesses(
        kind=WitnessKind.POINTS,
        points=points,
        variables=ring.names,
        quotient_degree=output.quotient_degree,
        completeness=Certainty.MODULAR,
        run=result,
    )


def _points_from(
    parametrization: RationalParametrization,
    chart,
    ideal: "Ideal",
) -> tuple[RationalPoint | AlgebraicPoint, ...]:
    if len(parametrization.w_ascending) - 1 > parametrization.quotient_degree:
        raise MsolveOutputError(
            "msolve printed an eliminating polynomial of larger degree than "
            "the quotient ring; refusing to interpret it"
        )
    generators = [g for g in ideal.gens if not g.is_zero()]

    points: list[RationalPoint | AlgebraicPoint] = []
    for min_poly, multiplicity in factor_univariate(parametrization.w_ascending):
        if multiplicity != 1:
            raise MsolveOutputError(
                "the eliminating polynomial is not squarefree; msolve 0.10.x "
                "parametrizes the radical, so this is corrupt or drifted output"
            )
        try:
            coordinates = tuple(
                coordinate_mod(
                    numerator.v_ascending,
                    numerator.denominator_scale,
                    parametrization.wprime_ascending,
                    min_poly,
                )
                for numerator in chart.numerators_input
            )
        except ValueError as exc:
            raise MsolveOutputError(str(exc)) from None

        for generator in generators:
            if not vanishes_at_algebraic(generator.terms(), coordinates, min_poly):
                raise MsolveOutputError(
                    f"witness verification failed: a point of the "
                    f"parametrization does not satisfy the generator "
                    f"{generator}. msolve's parametrization is wrong or was "
                    f"misread; nothing from this run can be trusted."
                )

        if len(min_poly) - 1 == 1:
            points.append(
                RationalPoint(
                    coordinates=tuple(vector[0] for vector in coordinates)
                )
            )
        else:
            points.append(
                AlgebraicPoint(
                    min_poly_ascending=min_poly, coordinates=coordinates
                )
            )
    return tuple(sorted(points, key=_point_sort_key))


def _point_sort_key(point: RationalPoint | AlgebraicPoint):
    if isinstance(point, RationalPoint):
        return (1, (), point.coordinates)
    return (point.degree, point.min_poly_ascending, point.coordinates)
