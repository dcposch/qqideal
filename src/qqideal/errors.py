"""Exceptions raised by qqideal.

Two rules decide which exception a failure gets:

* A caller mistake -- an unparseable polynomial, a ring mismatch, a request
  outside the v0.1 slate -- raises. Callers should not have to inspect a
  :class:`~qqideal.Verdict` to discover they wrote ``x/2``.
* A solver failure -- msolve timing out, dying, or emitting bytes we refuse to
  interpret -- becomes a :class:`~qqideal.Verdict` with
  :attr:`~qqideal.Kind.TIMEOUT` or :attr:`~qqideal.Kind.ERROR`.

Input errors are reported as msolveio's :class:`~msolveio.MsolveInputError`
rather than a new type, so that a bad polynomial raises the same exception
whether qqideal or msolveio caught it.
"""

from __future__ import annotations

from msolveio import (
    MsolveAmbiguous,
    MsolveDied,
    MsolveError,
    MsolveInputError,
    MsolveOutputError,
    MsolveTimeout,
    MsolveVersionUnsupported,
)

__all__ = [
    "QQIdealError",
    "RingMismatch",
    "MsolveError",
    "MsolveInputError",
    "MsolveOutputError",
    "MsolveAmbiguous",
    "MsolveTimeout",
    "MsolveDied",
    "MsolveVersionUnsupported",
]


class QQIdealError(Exception):
    """Base class for errors qqideal raises on its own behalf."""


class RingMismatch(QQIdealError, TypeError):
    """Two polynomials from different rings were combined.

    qqideal never coerces across rings implicitly: ``QQ[x]`` and ``QQ[x,y]``
    are different rings, and so are ``QQ[x,y]`` and ``F_p[x,y]``. Use
    :meth:`qqideal.Ring.embed` to move a polynomial deliberately.
    """
