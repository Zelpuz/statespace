``statespace_kinds``
====================

Kind parameters, constants and status codes shared by the library.

Kinds and constants
-------------------

.. code-block:: fortran

   integer, parameter :: dp = real64
   real(dp), parameter :: log2pi = 1.8378770664093454835606594728112_dp

``dp``
   The real kind of all computations.
``log2pi``
   :math:`\log 2\pi`, for the Gaussian log density.

Status codes
------------

Routines that can fail set an ``info`` argument to one of these. The C
interface returns them, and the Python package raises
:class:`~ssfortran.StateSpaceError` with the code.

=========================  =====  ==============================================
name                       value  meaning
=========================  =====  ==============================================
``SS_OK``                  0      success
``SS_ERR_DIM``             1      inconsistent array dimensions or arguments
``SS_ERR_NOT_PD``          2      a matrix is not positive definite
``SS_ERR_INIT``            3      missing or invalid initialization
``SS_ERR_UNSUPPORTED``     4      the model is outside the routine's scope
``SS_ERR_SINGULAR``        5      singular linear system
``SS_ERR_NOT_STATIONARY``  6      the transition matrix has an eigenvalue on
                                  or outside the unit circle
``SS_ERR_NOT_CONVERGED``   7      an iteration did not converge
=========================  =====  ==============================================

``SS_ERR_NOT_PD`` from a likelihood evaluation usually means the parameters
are invalid, for example a negative variance; the optimizer treats it as a
failed step and backs off.
