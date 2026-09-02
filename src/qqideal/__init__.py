"""Exact ideals over Q, with python-flint for arithmetic and msolve for
Groebner bases.

qqideal answers three questions about an ideal -- is its variety empty, what is
its dimension, and (when the dimension is zero) what is its degree -- and
refuses the rest. It is not a computer algebra system: anything outside that
slate raises :class:`NotImplementedError` rather than returning something
approximate.

msolve is called in Groebner mode (``-g``) and, for witness extraction over Q,
rational-parametrization mode (``-P``) through msolveio's mode-required
parsers. Solver mode is never invoked, so no output whose meaning depends on
the mode is ever guessed at.
"""

from __future__ import annotations

from .doublepoint import double_point_ideal
from .errors import MsolveInputError, QQIdealError, RingMismatch
from .ideal import Ideal
from .ring import Poly, Ring
from .verdict import Certainty, Kind, Verdict, ideal_verdict
from .witness import (
    AlgebraicPoint,
    RationalPoint,
    Witnesses,
    WitnessKind,
    witness_points,
)

__version__ = "0.2.0"

__all__ = [
    "Ring",
    "Poly",
    "Ideal",
    "Kind",
    "Certainty",
    "Verdict",
    "ideal_verdict",
    "double_point_ideal",
    "WitnessKind",
    "RationalPoint",
    "AlgebraicPoint",
    "Witnesses",
    "witness_points",
    "QQIdealError",
    "RingMismatch",
    "MsolveInputError",
]
