"""Verdicts: what was established, and how firmly.

A :class:`Verdict` is deliberately not a boolean. ``EMPTY`` and ``NONEMPTY`` are
answers; ``TIMEOUT`` and ``ERROR`` are non-answers, and a truthiness test would
quietly fold the non-answers into one of the answers -- in whichever direction
the surrounding ``if`` happens to lean. So :meth:`Verdict.__bool__` raises, and
callers branch on :attr:`Verdict.kind`.
"""

from __future__ import annotations

import enum
from dataclasses import dataclass
from typing import TYPE_CHECKING, Sequence

if TYPE_CHECKING:  # pragma: no cover - import cycle broken at runtime
    from .ideal import Ideal
    from .ring import Poly, Ring

__all__ = ["Kind", "Certainty", "Verdict", "ideal_verdict"]


class Kind(enum.Enum):
    """What the verdict says about the variety of the ideal."""

    #: The ideal is the unit ideal: no solutions in an algebraic closure.
    EMPTY = "empty"
    #: The ideal is proper: the variety is nonempty over an algebraic closure.
    NONEMPTY = "nonempty"
    #: msolve exceeded the timeout. Nothing was established.
    TIMEOUT = "timeout"
    #: msolve died or emitted bytes we refuse to interpret. Nothing was
    #: established.
    ERROR = "error"


class Certainty(enum.Enum):
    """How firmly the verdict was established."""

    #: True over the stated field.
    PROVEN = "proven"
    #: True for one modular prime msolve chose, and not lifted to Q.
    MODULAR = "modular"


@dataclass(frozen=True)
class Verdict:
    """The result of asking msolve about an ideal.

    :param kind: emptiness, or the absence of an answer.
    :param dim: Krull dimension of the quotient, or ``None`` when the ideal is
        the unit ideal or nothing was established.
    :param degree: number of solutions with multiplicity. Only set when
        ``kind`` is :attr:`Kind.NONEMPTY` and ``dim`` is ``0``.
    :param certainty: see :class:`Certainty`. When ``kind`` is
        :attr:`Kind.TIMEOUT` or :attr:`Kind.ERROR` nothing was established at
        all; the field is then set to the weaker value,
        :attr:`Certainty.MODULAR`, so no caller can read a stronger claim out of
        it than was made.
    :param msolve_version: the binary that produced this, when one ran.
    :param detail: why msolve timed out or died, for the non-answers.
    """

    kind: Kind
    dim: int | None
    degree: int | None
    certainty: Certainty
    msolve_version: str | None = None
    detail: str | None = None

    def __bool__(self) -> bool:
        raise TypeError(
            "Verdict has no truth value; test .kind explicitly, e.g. "
            "`if verdict.kind is Kind.EMPTY:`. TIMEOUT and ERROR are not "
            "answers and must not collapse into one."
        )

    def __str__(self) -> str:
        parts = [f"kind={self.kind.value}", f"certainty={self.certainty.value}"]
        if self.dim is not None:
            parts.append(f"dim={self.dim}")
        if self.degree is not None:
            parts.append(f"degree={self.degree}")
        if self.detail:
            parts.append(f"detail={self.detail}")
        return f"Verdict({', '.join(parts)})"


def ideal_verdict(
    gens: "Ideal | Sequence[Poly | str]",
    *,
    ring: "Ring | None" = None,
    opens: "Sequence[Poly | str]" = (),
    timeout: float = 60,
) -> Verdict:
    """Decide emptiness -- and dimension and degree -- for an ideal.

    :param gens: an :class:`~qqideal.Ideal`, or generators as
        :class:`~qqideal.Poly` or strings together with ``ring=``.
    :param ring: the ring, required when ``gens`` contains strings.
    :param opens: polynomials that must not vanish. The ideal is saturated at
        their product before the test, so the verdict is about
        ``V(I) \\ V(f1*..*fk)`` -- "does the system have a solution away from
        these hypersurfaces".
    :param timeout: wall-clock limit for msolve, in seconds.
    :returns: a :class:`Verdict`. Never raises for a solver failure; a msolve
        timeout or crash comes back as :attr:`Kind.TIMEOUT` or
        :attr:`Kind.ERROR`. Bad input still raises.
    """
    from .ideal import Ideal

    ideal = gens if isinstance(gens, Ideal) else Ideal(gens, ring=ring)
    if opens:
        product = ideal.ring.constant(1)
        for open_set in opens:
            product = product * ideal._coerce(open_set)
        ideal = ideal.saturate(product)
    return ideal.verdict(timeout=timeout)
