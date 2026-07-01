!======================================================================
! Authors: Vladimir Meshkov, Asen Pashov, Andrey Stolyarov
! Corresponding author: Vladimir Meshkov
! Contact e-mail: MeshkovVV@my.msu.ru
! Version: 1.0
!
! SPDX-License-Identifier: MIT
! Copyright (c) 2026 Vladimir Meshkov, Asen Pashov, Andrey Stolyarov
! qsolver is distributed under the MIT License; see the LICENSE file
! distributed with this program for the full license text.
!
!  This file consists of the module qsolver_mod. The module exposes a single
!  public subroutine, qsolver, supported by private subroutines and
!  functions that implement the numerical method.
!
!  qsolver -- a program for the accurate solution of coupled radial
!             Schrodinger equations for bound and resonant states
!
!  This Fortran 2003 program solves the coupled-channel (CC) radial
!  Schrodinger system on the semi-infinite interval r in [0, +infinity)
!  (or, optionally, Dirichlet boundary-value problems on a general interval
!  [a, b], including a = -infinity or b = +infinity).
!  The CC equations may include first-derivative couplings collected in Bmn(r)
!  and potential/coupling matrix elements Vmn(r). Both real- and complex-valued
!  functions are supported.
!
!----------------------------------------------------------------------
!  Problem description
!----------------------------------------------------------------------
!  The coupled-channel radial Schrodinger equation solved by qsolver can be written as
!
!      -(hbar^2/(2*mu)) d^2/dr^2  psi(r)
!      + (hbar^2/(2*mu)) [ B(r) d/dr psi(r) - lambda d/dr ( B_dag(r) psi(r) ) ]
!      + V(r) psi(r)  =  E G(r) psi(r),
!
!  where psi(r) is an M-component channel vector, V(r) and B(r) are MxM matrices,
!  G(r) is a diagonal matrix in the channel basis, B_dag(r) is the conjugate
!  transpose of B(r), and lambda = 0 or 1.
!  The B-term can be handled in both Hermitian form (lambda = 1)
!  and non-Hermitian form (lambda = 0).
!  If G is not provided, it is assumed that G = I.
!
!  The subroutine currently supports three types of boundary conditions:
!     1) BC_DIRICHLET at either endpoint: psi_m(r_b) = 0.
!     2) BC_REGULAR at either endpoint: no zero value is imposed and the
!        solution is treated as finite/regular at the corresponding endpoint.
!     3) For resonant states, the boundary condition in open
!        channels is an outgoing wave exp(i*k_m*r) for r -> +infinity.
!  The left and right endpoint conditions are selected independently. The
!  endpoints rmin and rmax define the interval; omitted endpoints represent
!  -infinity and +infinity, respectively.
!
!  The outgoing-wave boundary condition leads to complex
!  eigenvalues of the Hamiltonian,
!     E = Er - i*Gamma/2,
!  where Er is the resonance energy and Gamma is the resonance width.
!
!----------------------------------------------------------------------
!  Numerical method (short)
!----------------------------------------------------------------------
!  1) Exterior complex scaling (ECS):
!        r(y) = y                          , y <= rc
!             = rc + (y-rc) * exp(i*phi)   , y > rc
!     which converts the outgoing-wave boundary condition into a zero boundary
!     condition at infinity along the complex path (for a suitable phi).
!
!  2) Mapping to a finite domain:
!     A smooth monotone mapping x(y) is used to transform an infinite radial
!     domain to a finite interval. With an exact analytical transformation,
!     the coupled-channel problem is rewritten on a finite domain and
!     solved by polynomial collocation.
!
!  3) Spectral collocation:
!     The discretization exhibits exponential (spectral) convergence for
!     analytic solutions. The Legendre grid type is selected automatically
!     from the endpoint boundary conditions (Gauss-Lobatto, Gauss-Radau, or
!     Gauss-Legendre).
!     If ECS is applied, the interval is split at xc = x(rc)
!     to handle the derivative discontinuity introduced by ECS; Hermite
!     interpolation is used at the matching point, preserving spectral
!     convergence.
!
!  4) Eigenvalue solution:
!     The collocation discretization yields a (generally complex) matrix
!     eigenvalue problem, which is solved with ARPACK (a library for finding
!     selected eigenvalues) or LAPACK (a library for standard dense-matrix
!     linear algebra).
!
!----------------------------------------------------------------------
!  Additional features
!----------------------------------------------------------------------
!  * Automatic extrapolation to the complex plane:
!    For resonant calculations, if the user supplies Vmn(r) (and Bmn(r))
!    only on the real axis, qsolver extrapolates these functions into the
!    complex plane using barycentric interpolation on Gauss-Legendre nodes
!    with an optimally truncated polynomial degree. This makes it unnecessary to provide an
!    explicit analytic continuation of real-axis potentials/couplings.
!
!  * Spurious ("ghost") resonance filtering:
!    Spurious eigenvalues are detected by examining the eigenfunction in a
!    Legendre spectral basis: for physical solutions the Legendre coefficients
!    decay rapidly, whereas spurious solutions typically do not. A simple
!    decay-ratio threshold is used to discard non-physical eigenpairs.
!
!  * Energy derivatives:
!    First derivatives of (possibly complex) eigenvalues with respect to
!    Hamiltonian parameters can be computed. For complex-scaled problems, a
!    non-Hermitian analogue of the Hellmann-Feynman theorem is used; for
!    real (Hermitian) problems, the standard theorem applies. The user provides
!    derivative subroutines for dVmn/dp, dBmn/dp, and dGm/dp.
!
!----------------------------------------------------------------------
!  External libraries
!----------------------------------------------------------------------
!  The program uses the following ARPACK subroutines:
!      dsaupd, dseupd, dnaupd, dneupd, znaupd, zneupd
!
!  The program uses the following BLAS subroutines:
!      dcopy, zcopy
!
!  The program uses the following LAPACK routines:
!      dlamch, dsyevr, zheevr, dgeev, zgeev,
!      dgetrf, zgetrf, dgetrs, zgetrs, dsytrf,
!      zhetrf, zsytrf, dsytrs, zhetrs, zsytrs
!
!  All explicit interfaces use the qsolver working kind dp.
!  This assumes that dp corresponds to the double-precision kind expected
!  by the linked BLAS/LAPACK/ARPACK libraries.
!
!----------------------------------------------------------------------
!  Important
!----------------------------------------------------------------------
!  Most arguments are OPTIONAL, so it is not necessary to know all features
!  in order to perform simple tasks.
!  qsolver SHOULD BE CALLED WITH KEYWORD ARGUMENTS. Arguments may be supplied
!  in any order.
!
!======================================================================
!  Example: simplest program to calculate 10 eigenvalues of an infinite box
!======================================================================
!  use qsolver_mod
!  .....
!  call qsolver(rmin = 0.0d0, & ! the problem is solved on [0, pi]
!       &       rmax = pi,    & !
!       &       mu   = 0.5d0, & ! effective mass (see description of mu)
!       &       N    = 40,    & ! number of collocation grid points
!       &       nE   = 10,    & ! number of eigenvalues to be found
!       &       E    = Ev     ) ! Ev is the output array with energies
!
!======================================================================
!  See file qtest.f90 for test examples.
!  The distributed source files qsolver.f90, qpotentials.f90, and qtest.f90
!  have been validated to conform to the Fortran 2003 standard.
!======================================================================
!  Arguments of qsolver
!======================================================================
!  Name   Intent  Description
!
!  rmin   in   real(dp), optional
!              Left boundary. If rmin is not present, it is assumed that
!              rmin = -infinity (problem on an unbounded interval).
!
!  bc_left in  integer, optional, default BC_DIRICHLET
!              Boundary condition at rmin:
!                 BC_DIRICHLET = homogeneous Dirichlet, psi_m(rmin) = 0
!                 BC_REGULAR   = finite/regular value allowed at rmin
!
!  rmax   in   real(dp), optional
!              Right boundary. If rmax is not present, it is assumed that
!              rmax = +infinity.
!
!  bc_right in integer, optional, default BC_DIRICHLET
!              Boundary condition at rmax:
!                 BC_DIRICHLET = homogeneous Dirichlet, psi_m(rmax) = 0
!                 BC_REGULAR   = finite/regular value allowed at rmax
!              With complex scaling (phi present) the right boundary must be
!              BC_DIRICHLET; specifying anything else returns an error.
!
!  N      in   integer
!              Number of collocation grid points.
!              Increasing N usually increases accuracy; for analytic V and B,
!              the convergence is exponential.
!
!  mu     in   real(dp)
!              "Effective mass" in the Schrodinger equation.
!              The kinetic term is -(hbar^2/(2*m)) d^2/dr^2.
!              Here the input mu is defined as mu = m / hbar^2, so the
!              coefficient becomes -(1/(2*mu)) d^2/dr^2.
!
!  JRot   in   integer, optional, default JRot = 0
!              Rotational quantum number J. If JRot is present and JRot /= 0,
!              the centrifugal term J*(J+1)/(2*mu*r*r) is added to each
!              diagonal element V_mm.
!              (The user may also include it in V explicitly; the result is the same.)
!
!  nEq    in   integer, optional, default nEq = 1
!              Number of coupled equations (channels). If not present, nEq = 1.
!
!  nEf    out  integer, optional
!              On successful return: the number of eigenvalues actually returned.
!              On a critical error: nEf = 0.
!              On a warning: nEf may be nonzero.
!              If nEf < requested nE, this usually means that the spurious
!              eigenvalue filter removed some states. This is especially
!              important for resonant states, where many spurious eigenvalues
!              may appear near genuine resonances.
!              The filter can be turned off by setting f_thr < 0.
!              The filter is off by default if no complex scaling is used
!              (that is, if phi is not present).
!
!  E      out  real(dp), dimension(nEf), optional
!              Output array with energies.
!              For complex problems, it contains the real parts of the eigenvalues.
!              For complex scaling (phi present), the eigenvalues are complex:
!              Eig = E - i*Gamma/2.
!
!  W      out  real(dp), dimension(nEf), optional
!              Output array with minus twice the imaginary parts of the eigenvalues.
!              For a complex-scaled Hamiltonian, W is the resonance width Gamma.
!
!----------------------------------------------------------------------
!  Selecting eigenvalues by target / in an interval / by indices
!  If none of these parameters (including nE) is present, then qsolver computes all
!  N*nEq eigenvalues of the Hamiltonian matrix. The subroutine may return
!  fewer eigenvalues than requested if the spurious-eigenvalue filter is enabled
!  (it is enabled by default when resonant states are requested).
!
!  Eigenvalue sorting behavior:
!  By default, if the target value (E0) is not provided, the returned
!  eigenstates are sorted by energy in ascending order. If E0 is provided,
!  the eigenvalues are sorted by increasing distance from the target,
!  |E - E0| (or |Ec-Ec0| in the complex case, where Ec0 = E0 - i*W0/2).
!  This default behavior can be changed with the sortE and sortW keywords.
!----------------------------------------------------------------------
!
!  nE     in   integer, optional
!              On entry: number of eigenvalues to be found.
!              If nE is not present and E0 is present, it is assumed that nE = 1.
!              If nE is present and E0 is not, then the lowest nE eigenvalues
!              are requested. This is equivalent to setting imin = 1 and imax = nE.
!
!              If nE < 0, then all N*nEq eigenvalues of the Hamiltonian
!              matrix are computed (this can be expensive).
!
!              If nE = 0, qsolver returns immediately with nEf = 0.
!
!  E0     in   real(dp), optional
!              If E0 is present, qsolver tries to compute nE eigenvalues
!              closest to the target value E = E0 - i*W0/2 (see W0).
!
!  W0     in   real(dp), optional, default W0 = 0.0_dp
!              Target width parameter used together with E0:
!              the target eigenvalue is E0 - i*W0/2. If W0 is not present,
!              then W0 = 0.
!
!  Emin   in   real(dp), optional
!              Lower bound of the energy interval to be searched (Emin < Emax).
!              If an energy interval is provided, then for Hermitian problems
!              the LAPACK solver can effectively search only for eigenvalues in this interval.
!              For non-Hermitian problems, the bounds are used for sorting/filtering.
!
!  Emax   in   real(dp), optional
!              Upper bound of the energy interval to be searched.
!
!  Wmax   in   real(dp), optional
!              Upper bound for W = -2*Im(E).
!              Eigenvalues with W > Wmax are excluded from the returned list.
!
!  imin   in   integer, optional
!              Lower index of the eigenvalues to be returned.
!              Together with imax, it defines the index range imin..imax.
!              For Hermitian problems without E0, qsolver may use a fast
!              LAPACK index-range algorithm to compute only the requested
!              eigenvalues. Otherwise, the range imin..imax is taken from
!              the ordered list described above.
!              If some states from the requested range are removed by
!              Emin/Emax, Wmax, or by the spurious-state filter, fewer
!              eigenvalues may be returned. The spurious-state filter can
!              be turned off by setting f_thr < 0.
!              1 <= imin <= imax <= N*nEq.
!
!  imax   in   integer, optional
!              Upper index of the eigenvalues to be returned.
!              Together with imin, it defines the index range imin..imax.
!              For Hermitian problems without E0, qsolver may use a fast
!              LAPACK index-range algorithm to compute only the requested
!              eigenvalues. Otherwise, the range imin..imax is taken from
!              the ordered list described above.
!              If some states from the requested range are removed by
!              Emin/Emax, Wmax, or by the spurious-state filter, fewer
!              eigenvalues may be returned. The spurious-state filter can
!              be turned off by setting f_thr < 0.
!              1 <= imin <= imax <= N*nEq.
!
!  sortE  in   integer, optional
!              Energy sorting mode:
!                0 = do not sort
!                1 = ascending
!                2 = descending
!              For complex problems, sorting by energy means sorting by Re(E).
!              If E0 is present and sortE is not provided, target sorting
!              by |E - E0| is used instead.
!              If E0 is not present and sortE is not provided, the default is 1.
!
!  sortW  in   integer, optional
!              Width sorting mode for complex problems:
!                0 = do not sort
!                1 = ascending by W = -2*Im(E)
!                2 = descending by W = -2*Im(E)
!              Ignored for Hermitian problems.
!              If E0 is present and sortW is not provided, target sorting
!              may be used (see above).
!              If E0 is not present and sortW is not provided, the default is 0.
!              sortE and sortW cannot both be nonzero.
!
!----------------------------------------------------------------------
!  Potential matrix V and derivative-coupling matrix B
!----------------------------------------------------------------------
!  V      in   real(dp), optional, V(r, k, m)
!              k, m :: integer
!              r    :: real(dp)
!              Function V(r, k, m) returns V_km(r) for real r on (rmin, rmax).
!              By default, it is assumed that V(k, m) = V(m, k) (symmetric),
!              which enables faster solvers (if complex scaling is not used).
!              If V is not symmetric, set is_V_sym = .false.
!              If V is not present, it is assumed that V = 0.
!
!  Vc     in   complex(dp), optional, Vc(r, k, m)
!              k, m :: integer
!              r    :: complex(dp)
!              Same as V, but complex-valued.
!              For complex scaling, it must be valid for complex r.
!              If complex scaling is not used, then it is assumed that
!              V(k, m) = conjg(V(m, k)) (Hermitian), which enables faster solvers.
!              If Vc is not Hermitian, set is_V_sym = .false.
!              If Vc is absent, qsolver may extrapolate V
!              into the complex plane when complex scaling is used.
!
!  V1     in   real(dp), optional, V1(r)
!              r :: real(dp)
!              Same as V for nEq = 1 (single-channel case).
!
!  V1c    in   complex(dp), optional, V1c(r)
!              r :: complex(dp)
!              Same as Vc for nEq = 1 (single-channel case).
!
!  B      in   real(dp), optional, B(r, k, m)
!              k, m :: integer
!              r    :: real(dp)
!              Function B(r, k, m) returns the first-derivative coupling matrix B_km(r).
!              If B is not present, it is assumed that B = 0.
!              By default, B is assumed to enter the Schrodinger equation
!              in Hermitian form. If B enters as a simple B*d/dr term,
!              set is_B_sym = .false.
!
!  Bc     in   complex(dp), optional, Bc(r, k, m)
!              k, m :: integer
!              r    :: complex(dp)
!              Same as B, but complex-valued.
!              For complex scaling, it must be valid for complex r.
!              If Bc is absent, qsolver may extrapolate B
!              into the complex plane when complex scaling is used.
!
!  B1     in   real(dp), optional, B1(r)
!              r :: real(dp)
!              Same as B for nEq = 1 (single-channel case).
!
!  B1c    in   complex(dp), optional, B1c(r)
!              r :: complex(dp)
!              Same as Bc for nEq = 1 (single-channel case).
!
!  G     in   real(dp), optional, G(r, k)
!              k :: integer
!              r :: real(dp)
!              Diagonal matrix function G_k(r) in the generalized equation
!              H*psi = E*G*psi for the multi-channel case.
!              If G is not present, it is assumed that G_k(r) = 1.
!
!  Gc    in   complex(dp), optional, Gc(r, k)
!              k :: integer
!              r :: complex(dp)
!              Same as G, but complex-valued.
!              For complex scaling, it must be valid for complex r.
!              If Gc is absent, qsolver may extrapolate G
!              into the complex plane when complex scaling is used.
!
!  G1    in   real(dp), optional, G1(r)
!              r :: real(dp)
!              Same as G for the single-channel case.
!
!  G1c   in   complex(dp), optional, G1c(r)
!              r :: complex(dp)
!              Same as Gc for the single-channel case.
!
!  is_V_sym  in   logical, optional, default .true.
!                 If .true., V_km(r) is assumed symmetric/Hermitian
!                 (V_km = conjg(V_mk) in the complex case).
!
!  is_B_sym  in   logical, optional, default .true.
!                 If .true., B is assumed to be used in Hermitian form.
!
!----------------------------------------------------------------------
!  Complex scaling (resonances)
!----------------------------------------------------------------------
!  phi    in   real(dp), optional
!              If phi is present, complex scaling is applied.
!              phi is the complex-scaling angle.
!              If phi is not present, no complex scaling is used.
!
!  rc     in   real(dp), optional
!              Exterior complex-scaling radius rc: scaling starts at r = rc.
!              If phi is present but rc is not present, simple complex scaling
!              is used (equivalent to rc = rmin).
!              If rmin is not present, rc must be specified when complex scaling is used.
!
!----------------------------------------------------------------------
!  Mapping parameters
!  Mapping is typically used to transform an infinite interval to a finite one,
!  which facilitates an efficient discretization of the problem.
!  It also helps to concentrate collocation points in the region where the wave
!  function oscillates rapidly, which can greatly reduce the required number
!  of collocation points (problem dimension).
!----------------------------------------------------------------------
!  mtype  in   character(:), optional
!              Type of coordinate transformation (mapping).
!              If mtype is not present or mtype = 'none', no mapping is used.
!
!                 Available one-parameter mappings x = x(y, beta):
!
!                 mtype    x(y, beta)
!                 -------  ------------------------------------------
!                 none     x = y
!                 atan     x = (2/pi) * atan(beta*y)
!                 rat      x = (y**beta - 1) / (y**beta + 1)
!                 tanh     x = tanh(beta*y)
!                 dexp     x = tanh( sinh(beta*y) )
!                 exp      x = 1 - exp( -y**beta )
!
!                 Available two-parameter mappings x = x(y, beta, ye):
!
!                 mtype    x(y, beta, ye)
!                 -------  -----------------------------------------------
!                 none     x = y
!                 atan     x = (2/pi) * atan( beta*(y/ye - 1) )
!                 rat      x = ((y/ye)**beta - 1) / ((y/ye)**beta + 1)
!                 tanh     x = tanh( beta*(y/ye - 1) )
!                 dexp     x = tanh( sinh( beta*(y/ye - 1) ) )
!                 exp      x = 1 - exp( -(y/ye)**beta )
!
!  beta   in   real(dp), optional
!              Mapping parameter. For diatomic potentials with long-range
!              behavior V(r) ~ -Cn/r^n, for 'atan' and 'rat'
!              a useful initial choice is beta = n/2 - 1.
!              For complex-scaled (resonance) problems, smaller beta often works better.
!
!  ye     in   real(dp), optional
!              Mapping parameter for two-parameter mappings. If ye is not present,
!              a one-parameter mapping is applied. In the single-channel case, best results
!              (faster convergence) are usually obtained when ye is close to the
!              equilibrium distance (minimum of V(r)).
!
!----------------------------------------------------------------------
!  Spurious eigenvalue filter (mainly for resonances)
!----------------------------------------------------------------------
!  f_thr  in   real(dp), optional, default f_thr = 1.0e-6
!              Spurious-eigenvalue threshold. By default, it is enabled
!              only for resonant-state calculations (phi keyword is present).
!              The program expands the wavefunction into a Legendre series.
!              For a physical eigenstate, |c_3last|/|c_max| < f_thr,
!              where c_3last is the average of the last three coefficients;
!              otherwise the eigenstate is excluded.
!              To enable the filter when complex scaling is not used, simply set
!              f_thr to a positive value.
!              If f_thr < 0, then eigenvalue filtering is disabled.
!
!----------------------------------------------------------------------
!  Grid points and quadrature weights
!----------------------------------------------------------------------
!
!  rv     out  real(dp), dimension(N), optional
!              Grid points r_i corresponding to the collocation nodes
!              (use together with qw for integrals; see Remarks below).
!  rvc    out  complex(dp), dimension(N), optional
!              Complex grid points r_i corresponding to the collocation nodes.
!              A complex grid appears when complex scaling is used to calculate resonant states.
!
!  qw     out  real(dp), dimension(N), optional
!              Quadrature weights for numerical integration with the obtained
!              wavefunctions on the qsolver grid.
!  qwc    out  complex(dp), dimension(N), optional
!              Complex quadrature weights for numerical integration with the obtained
!              wavefunctions on the qsolver grid. Complex weights appear when complex
!              scaling is used to calculate resonant states, because
!              Int ... dr = Int ... (dr/dy)*(dy/dx) dx, where dr/dy is complex outside rc.
!
!----------------------------------------------------------------------
!  Output eigenfunctions (right/ket and left/bra, real/complex)
!  To compute matrix elements f_km = <psi_k|f(r)|psi_m>,
!  use code such as
!
!  (real problems)
!     f_km = 0.0_dp
!     do j = 1, nEq
!        f_km = f_km + dot_product(braW(:, j, k), f(rv(:)) * WF(:, j, m))
!     end do
!
!  (complex problems)
!     f_km = (0.0_dp, 0.0_dp)
!     do j = 1, nEq
!        f_km = f_km + dot_product(braWc(:, j, k), f(rvc(:)) * WFc(:, j, m))
!     end do
!
!----------------------------------------------------------------------
!  WF     out  real(dp), dimension(N, nEq, nEf), optional
!              WF(i,m,k) contains wavefunction psi_m at grid point i
!              for the k-th state.
!              For real Hermitian problems, the normalization is
!              sum_m sum_i WF(i,m,k)**2 * qw(i) = 1.
!              Equivalently:
!              sum_m dot_product(braW(:,m,k), WF(:,m,k)) = 1.
!              For complex or non-symmetric problems, if WF is provided, it
!              contains the real part of the right eigenvector.
!
!  WF1    out  real(dp), dimension(N, nEf), optional
!              Same as WF for the single-channel case (nEq = 1).
!
!  WFc    out  complex(dp), dimension(N, nEq, nEf), optional
!              Complex wavefunctions (right eigenvectors, ket-vectors).
!              Normalizations:
!              - No complex scaling (phi absent):
!                sum_m sum_i conjg(WFc(i,m,k))*WFc(i,m,k)*qw(i) = 1.
!              - Complex scaling case (phi present):
!                right vectors are c-normalized:
!                sum_m sum_i WFc(i,m,k)**2*qwc(i) = 1.
!              In either case:
!              sum_m dot_product(braWc(:,m,k), WFc(:,m,k)) = 1.
!
!  WF1c   out  complex(dp), dimension(N, nEf), optional
!              Same as WFc for the single-channel case.
!
!  braW   out  real(dp), dimension(N, nEq, nEf), optional
!              Weighted bra vectors for real wavefunctions.
!              They are returned such that expectation values can be computed as
!              sum_m dot_product(braW(:,m,k), f(rv(:))*WF(:,m,k)) without explicit
!              quadrature weights.
!              For complex or non-Hermitian problems, use braWc/braW1c.
!
!  braW1  out  real(dp), dimension(N, nEf), optional
!              Same as braW for the single-channel case.
!
!  braWc  out  complex(dp), dimension(N, nEq, nEf), optional
!              Weighted bra vectors.
!              They are returned such that expectation values can be computed as
!              sum_m dot_product(braWc(:,m,k), f(rvc(:))*WFc(:,m,k)) without explicit
!              quadrature weights.
!
!  braW1c out  complex(dp), dimension(N, nEf), optional
!              Same as braWc for the single-channel case.
!
!----------------------------------------------------------------------
!  Energy / width derivatives
!----------------------------------------------------------------------
!  ndE    in   integer, optional, default = 1
!              Number of parameters p for which derivatives are requested.
!              Derivative subroutines return all ndE derivatives in an array
!              whose index p corresponds to derivative with respect to parameter p.
!              If ndE = 0, derivative calculations are disabled completely.
!              If ndE < 0, qsolver returns an input error.
!
!  dV     in   subroutine, optional, dV(r, k, m, dVdp)
!              real(dp), intent(in)  :: r
!              integer,  intent(in)  :: k, m
!              real(dp), intent(out) :: dVdp(:)
!              dVdp(p) = dV_km(r)/dp_p for real-valued multi-channel V.
!
!  dVc    in   subroutine, optional, dVc(r, k, m, dVdp)
!              complex(dp), intent(in)  :: r
!              integer,     intent(in)  :: k, m
!              complex(dp), intent(out) :: dVdp(:)
!              Same as dV, but for complex-valued multi-channel V.
!              For complex scaling, it must be valid for complex r.
!              If dVc is absent, qsolver may extrapolate dV
!              into the complex plane when complex scaling is used.
!
!  dV1    in   subroutine, optional, dV1(r, dVdp)
!              real(dp), intent(in)  :: r
!              real(dp), intent(out) :: dVdp(:)
!              Same as dV for the single-channel case.
!
!  dV1c   in   subroutine, optional, dV1c(r, dVdp)
!              complex(dp), intent(in)  :: r
!              complex(dp), intent(out) :: dVdp(:)
!              Same as dVc for the single-channel case.
!
!  dB     in   subroutine, optional, dB(r, k, m, dBdp)
!              real(dp), intent(in)  :: r
!              integer,  intent(in)  :: k, m
!              real(dp), intent(out) :: dBdp(:)
!              dBdp(p) = dB_km(r)/dp_p for real-valued multi-channel B.
!
!  dBc    in   subroutine, optional, dBc(r, k, m, dBdp)
!              complex(dp), intent(in)  :: r
!              integer,     intent(in)  :: k, m
!              complex(dp), intent(out) :: dBdp(:)
!              Same as dB, but for complex-valued multi-channel B.
!              For complex scaling, it must be valid for complex r.
!              If dBc is absent, qsolver may extrapolate dB
!              into the complex plane when complex scaling is used.
!
!  dB1    in   subroutine, optional, dB1(r, dBdp)
!              real(dp), intent(in)  :: r
!              real(dp), intent(out) :: dBdp(:)
!              Same as dB for the single-channel case.
!
!  dB1c   in   subroutine, optional, dB1c(r, dBdp)
!              complex(dp), intent(in)  :: r
!              complex(dp), intent(out) :: dBdp(:)
!              Same as dBc for the single-channel case.
!
!  dG     in   subroutine, optional, dG(r, k, dGdp)
!              real(dp), intent(in)  :: r
!              integer,  intent(in)  :: k
!              real(dp), intent(out) :: dGdp(:)
!              dGdp(p) = dG_k(r)/dp_p for real-valued multi-channel G.
!              If G depends on differentiation parameters, dG or dGc
!              should be provided. If dG/dGc is absent, G is assumed
!              independent of these parameters.
!
!  dGc    in   subroutine, optional, dGc(r, k, dGdp)
!              complex(dp), intent(in)  :: r
!              integer,     intent(in)  :: k
!              complex(dp), intent(out) :: dGdp(:)
!              Same as dG, but for complex-valued multi-channel G.
!              For complex scaling, it must be valid for complex r.
!              If dGc is absent, qsolver may extrapolate dG
!              into the complex plane when complex scaling is used.
!
!  dG1    in   subroutine, optional, dG1(r, dGdp)
!              real(dp), intent(in)  :: r
!              real(dp), intent(out) :: dGdp(:)
!              Same as dG for the single-channel case.
!
!  dG1c   in   subroutine, optional, dG1c(r, dGdp)
!              complex(dp), intent(in)  :: r
!              complex(dp), intent(out) :: dGdp(:)
!              Same as dGc for the single-channel case.
!
!  dE     out  real(dp), dimension(ndE, nEf), optional
!              On exit, dE(p,q) = dE_q/dp_p for eigenstate q.
!              For bound states, this is the energy derivative.
!              For resonances, this is the derivative of the resonance energy.
!
!  dW     out  real(dp), dimension(ndE, nEf), optional
!              On exit, dW(p,q) = dGamma_q/dp_p for eigenstate q.
!              For bound states, Gamma = 0 and dW is returned as zero.
!
!======================================================================
!  Error handling
!======================================================================
!  err_mode    in   integer, optional, default 1
!                   There are 3 error modes:
!                   err_mode = 1
!                      print all errors and warnings on screen
!                   err_mode = 2
!                      do not print any errors or warnings on screen
!                   err_mode = 3
!                      print all errors and warnings and halt.
!
!  err_code    out  integer, optional
!                   A decimal integer qpkji, where q, p, k, j, i are digits:
!                   q - number of ARPACK errors
!                   p - number of LAPACK errors
!                   k - number of internal errors
!                   j - number of allocation errors
!                   i - number of input errors
!
!  err_message out  character(:), allocatable, optional
!                   Deferred-length allocatable character variable.
!                   On exit, it is allocated and contains all accumulated
!                   error and warning messages joined into a single string,
!                   separated by '; '. If there are no messages, an empty
!                   string is returned.
!
!---------------------------------------------------------------------
!  Compilation example (Linux / macOS)
!----------------------------------------------------------------------
!  For the supplied test program qtest.f90, the modules qsolver_mod and
!  potentials must be available before qtest.f90 is compiled. One simple
!  command line is:
!  gfortran qsolver.f90 qpotentials.f90 qtest.f90 -llapack -lblas -larpack -o qtest
!  (add -L... if LAPACK/ARPACK are not in the default linker path)
!======================================================================

module qsolver_mod
   implicit none

   private

   ! qsolver, its version string, the working real kind, and boundary-condition constants are public
   public :: qsolver, qsolver_version, dp, BC_DIRICHLET, BC_REGULAR

   !================================================================
   ! Parameters
   !================================================================
   character(len=*), parameter :: qsolver_version = '1.0'
   integer, parameter  :: dp = selected_real_kind(15, 307)
   real(dp), parameter :: pi = acos(-1.0_dp)

   !================================================================
   ! Boundary condition constants
   !================================================================
   integer, parameter :: BC_DIRICHLET = 0
   integer, parameter :: BC_REGULAR   = 1

   ! spectral grid kinds (collocation node families on [-1,1])
   integer, parameter :: GRID_LEGENDRE = 1
   integer, parameter :: GRID_RADAU_M1 = 2   ! Gauss-Radau, endpoint -1 included
   integer, parameter :: GRID_RADAU_P1 = 3   ! Gauss-Radau, endpoint +1 included
   integer, parameter :: GRID_LOBATTO  = 4
   integer, parameter :: GRID_PLAIN    = 5   ! plain mesh: build only allocates storage, caller fills nodes

   !================================================================
   ! Error mode constants
   !================================================================
   integer, parameter :: ERRMODE_PRINTALL = 1   ! print errors+warnings, continue
   integer, parameter :: ERRMODE_NOWARN   = 2   ! silent, continue
   integer, parameter :: ERRMODE_STOP     = 3   ! print and halt immediately
   
   !================================================================
   ! Error code constants (decimal digit encoding: errcode = q p k j i)
   !   i - input errors       (+1)
   !   j - allocation errors  (+10)
   !   k - internal errors    (+100)
   !   p - LAPACK errors      (+1000)
   !   q - ARPACK errors      (+10000)
   ! Each digit is capped at 9.
   !================================================================
   integer, parameter :: ERR_INPUT    = 1
   integer, parameter :: ERR_ALLOC    = 10
   integer, parameter :: ERR_INTERNAL = 100
   integer, parameter :: ERR_LAPACK   = 1000
   integer, parameter :: ERR_ARPACK   = 10000
   
   ! maximum error/warning messages
   integer, parameter :: MAXERRMSG = 16

   !================================================================
   ! Error handling type
   !================================================================
   type err_type
       integer            :: mode    = ERRMODE_PRINTALL
       integer            :: errcode = 0
       integer            :: nerr    = 0
       integer            :: nwarn   = 0
       integer            :: nmsg    = 0
       logical            :: is_truncated = .false.
       character(len=256) :: messages(MAXERRMSG)
   contains
       procedure :: init         => qerr_init
       procedure :: clear        => qerr_clear
       procedure :: push_msg     => qerr_push_msg
       procedure :: push_err     => qerr_push_err
       procedure :: push_warning => qerr_push_warning
       procedure :: has_error    => qerr_has_error
       procedure :: finalize     => qerr_finalize
   end type err_type

   !===================================================
   ! Mapping class - defines mapping parameters and functions
   ! needed to transform equation to finite domain
   !===================================================
   type :: mapping_type
      integer  :: mn        ! mapping number
      real(dp) :: beta, ye  ! mapping parameters
      type(err_type), pointer :: err => null()
   contains
      procedure :: init => init_mapping
      procedure :: yx
      procedure :: xyInf
      procedure :: xyMInf
      procedure :: dydx
      procedure :: dxdyx
      procedure :: Fx
      procedure :: Qx
      procedure :: xy_r => xy       ! real version procedure is named xy
      procedure :: xy_c => xyc      ! complex version procedure is named xyc
      generic   :: xy => xy_r, xy_c ! now xy overloads two functions
   end type mapping_type

   !===================================================
   ! Wrapper class for real functions and derivative subroutines
   ! For example.
   ! Vf%set_par(m, k)  - set function parameters 
   ! Vf%f(r)           - evaluate scalar function at point r
   ! dVf%fv(r, dVdp)  - evaluate all derivatives at point r
   !===================================================
   type :: func_type
       type(mapping_type), pointer :: mapping => null()
       type(err_type),     pointer :: err     => null()

      integer :: n_par   ! number of parameters
      integer :: i, j, k ! integer parameter values
      integer :: n_vec = 0 ! length of vector returned by derivative subroutines

      logical :: is_associated = .false. ! is function associated

      ! pointers to real functions
      procedure(func1), pointer, nopass :: p1 => null() ! 1-arg function
      procedure(func2), pointer, nopass :: p2 => null() ! 2-arg function
      procedure(func3), pointer, nopass :: p3 => null() ! 3-arg function
      procedure(func4), pointer, nopass :: p4 => null() ! 4-arg function

      ! pointers to real derivative subroutines
      procedure(vfunc1), pointer, nopass :: v1 => null() ! f(x, df/dp)
      procedure(vfunc2), pointer, nopass :: v2 => null() ! f(x, i, df/dp)
      procedure(vfunc3), pointer, nopass :: v3 => null() ! f(x, i, j, df/dp)
   contains   
      ! initialization of real functions and derivative subroutines
      procedure :: init_func1  ! 1-arg function
      procedure :: init_func2  ! 2-arg function
      procedure :: init_func3  ! 3-arg function
      procedure :: init_func4  ! 4-arg function
      procedure :: init_vfunc1 ! derivative subroutine f(x, df/dp)
      procedure :: init_vfunc2 ! derivative subroutine f(x, i, df/dp)
      procedure :: init_vfunc3 ! derivative subroutine f(x, i, j, df/dp)
      procedure :: init_func   ! disassociate function

      ! setting function parameters
      procedure :: set_par1 => set_func_par1
      procedure :: set_par2 => set_func_par2
      procedure :: set_par3 => set_func_par3

      ! generic name (set_func_par) that covers both
      generic   :: set_par => set_par1, set_par2, set_par3

      ! setting mapping function
      procedure :: set_m       => set_func_mapping
      procedure :: unset_m     => unset_func_mapping
      generic   :: set_mapping => set_m, unset_m

      ! evaluate function or derivative vector
      procedure :: fr => func
      procedure :: fc => func_c
      generic   :: f  => fr, fc
      procedure :: fv => func_vec

      ! is function associated
      procedure :: is => is_func

   end type func_type

   type :: cfunc_type
       type(mapping_type), pointer :: mapping => null()
       type(err_type),     pointer :: err     => null()

      integer :: n_par     ! number of parameters
      integer :: i, j, k   ! integer parameter values
      integer :: n_vec = 0 ! length of vector returned by derivative subroutines

      logical :: is_associated = .false. ! is function associated?

      ! pointers to complex functions
      procedure(cfunc1), pointer, nopass :: p1 => null() ! 1-arg function
      procedure(cfunc2), pointer, nopass :: p2 => null() ! 2-arg function
      procedure(cfunc3), pointer, nopass :: p3 => null() ! 3-arg function
      procedure(cfunc4), pointer, nopass :: p4 => null() ! 4-arg function

      ! pointers to complex derivative subroutines
      procedure(cvfunc1), pointer, nopass :: v1 => null() ! f(x, df/dp)
      procedure(cvfunc2), pointer, nopass :: v2 => null() ! f(x, i, df/dp)
      procedure(cvfunc3), pointer, nopass :: v3 => null() ! f(x, i, j, df/dp)
   contains   
      ! initialization of complex functions and derivative subroutines
      procedure :: init_cfunc1  ! 1-arg function
      procedure :: init_cfunc2  ! 2-arg function
      procedure :: init_cfunc3  ! 3-arg function
      procedure :: init_cfunc4  ! 4-arg function
      procedure :: init_cvfunc1 ! derivative subroutine f(x, df/dp)
      procedure :: init_cvfunc2 ! derivative subroutine f(x, i, df/dp)
      procedure :: init_cvfunc3 ! derivative subroutine f(x, i, j, df/dp)
      procedure :: init_cfunc   ! disassociate function

      ! setting function parameters
      procedure :: set_par1 => set_cfunc_par1
      procedure :: set_par2 => set_cfunc_par2
      procedure :: set_par3 => set_cfunc_par3

      ! generic name (set_func_par) that covers both
      generic   :: set_par => set_par1, set_par2, set_par3

      ! evaluate function or derivative vector
      procedure :: fc => cfunc
      procedure :: fr => cfunc_r
      generic   :: f  => fc, fr
      procedure :: fv => cfunc_vec

      ! is function associated
      procedure :: is => is_cfunc

   end type cfunc_type

   !===================================================
   ! A spectral grid: collocation nodes together with the
   ! matching quadrature weights, kept as one object.
   !===================================================
   type :: grid_type
      integer :: kind = 0               ! GRID_LEGENDRE / GRID_RADAU_* / GRID_LOBATTO
      real(dp), allocatable :: x(:)     ! collocation nodes (on [-1,1] until set_scale)
      real(dp), allocatable :: qw(:)    ! quadrature weights
      real(dp), allocatable :: bary(:)  ! barycentric weights (scale-invariant)
      integer  :: lo = 0, hi = 0        ! active-node index range (formal boundaries excluded)
      real(dp) :: xmin                  ! current interval endpoints, set by build
      real(dp) :: xmax
      real(dp) :: length                ! = xmax - xmin (kept for convenience)
   contains
      procedure :: build     => grid_build
      procedure :: set_scale => grid_set_scale
   end type grid_type

   !===================================================
   ! A differentiation operator attached to a grid: holds D1, D2
   ! and a reference to the grid they are built on.
   !===================================================
   type :: diffmat_type
      type(grid_type), pointer :: grid => null()
      real(dp), allocatable :: D1(:,:), D2(:,:)
      real(dp), allocatable :: DH1(:,:), DH2(:,:)   ! Hermite-interpolation matrices
   contains
      procedure :: build         => diffmat_build
      procedure :: build_sym     => diffmat_build_sym
      procedure :: build_hermite => diffmat_build_hermite
      procedure :: set_scale     => diffmat_set_scale
   end type diffmat_type

   !===================================================
   ! Interfaces for function and derivative-subroutine shapes
   !===================================================
   abstract interface
      function func1(x) result(y)
         import :: dp
         real(dp), intent(in) :: x
         real(dp)             :: y
      end function func1

      function func2(x,k) result(y)
         import :: dp
         real(dp), intent(in) :: x
         integer, intent(in)  :: k
         real(dp)             :: y
      end function func2

      function func3(x, i, j) result(y)
         import :: dp
         real(dp), intent(in) :: x
         integer, intent(in)  :: i, j
         real(dp)             :: y
      end function func3

      function func4(x, i, j, k) result(y)
         import :: dp
         real(dp), intent(in) :: x
         integer, intent(in)  :: i, j, k
         real(dp)             :: y
      end function func4

      subroutine vfunc1(x, dydp)
         import :: dp
         real(dp), intent(in)  :: x
         real(dp), intent(out) :: dydp(:)
      end subroutine vfunc1

      subroutine vfunc2(x, i, dydp)
         import :: dp
         real(dp), intent(in)  :: x
         integer,  intent(in)  :: i
         real(dp), intent(out) :: dydp(:)
      end subroutine vfunc2

      subroutine vfunc3(x, i, j, dydp)
         import :: dp
         real(dp), intent(in)  :: x
         integer,  intent(in)  :: i, j
         real(dp), intent(out) :: dydp(:)
      end subroutine vfunc3

      function cfunc1(x) result(y)
         import :: dp
         complex(dp), intent(in) :: x
         complex(dp)             :: y
      end function cfunc1

      function cfunc2(x,k) result(y)
         import :: dp
         complex(dp), intent(in) :: x
         integer, intent(in)     :: k
         complex(dp)             :: y
      end function cfunc2

      function cfunc3(x, i, j) result(y)
         import :: dp
         complex(dp), intent(in) :: x
         integer, intent(in)     :: i, j
         complex(dp)             :: y
      end function cfunc3

      function cfunc4(x, i, j, k) result(y)
         import :: dp
         complex(dp), intent(in) :: x
         integer, intent(in)     :: i, j, k
         complex(dp)             :: y
      end function cfunc4

      subroutine cvfunc1(x, dydp)
         import :: dp
         complex(dp), intent(in)  :: x
         complex(dp), intent(out) :: dydp(:)
      end subroutine cvfunc1

      subroutine cvfunc2(x, i, dydp)
         import :: dp
         complex(dp), intent(in)  :: x
         integer,     intent(in)  :: i
         complex(dp), intent(out) :: dydp(:)
      end subroutine cvfunc2

      subroutine cvfunc3(x, i, j, dydp)
         import :: dp
         complex(dp), intent(in)  :: x
         integer,     intent(in)  :: i, j
         complex(dp), intent(out) :: dydp(:)
      end subroutine cvfunc3
   end interface

   !===================================================
   ! Explicit interfaces for BLAS/LAPACK/ARPACK
   !===================================================
   interface

      function dlamch(cmach) result(res)
         import :: dp
         character(len=1), intent(in) :: cmach
         real(dp) :: res
      end function dlamch

      subroutine dcopy(n, dx, incx, dy, incy)
         import :: dp
         integer, intent(in)       :: n, incx, incy
         real(dp), intent(in)      :: dx(*)
         real(dp), intent(inout)   :: dy(*)
      end subroutine dcopy

      subroutine zcopy(n, zx, incx, zy, incy)
         import :: dp
         integer, intent(in)        :: n, incx, incy
         complex(dp), intent(in)    :: zx(*)
         complex(dp), intent(inout) :: zy(*)
      end subroutine zcopy

      subroutine dsyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, &
         & isuppz, work, lwork, iwork, liwork, info)
         import :: dp
         character(len=1), intent(in) :: jobz, range, uplo
         integer, intent(in)          :: n, lda, il, iu, ldz, lwork, liwork
         real(dp), intent(inout)      :: a(lda,*)
         real(dp), intent(in)         :: vl, vu, abstol
         real(dp), intent(out)        :: w(*), z(ldz,*), work(*)
         integer, intent(out)         :: m, info
         integer, intent(out)         :: isuppz(*), iwork(*)
      end subroutine dsyevr

      subroutine zheevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, &
         & isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)
         import :: dp
         character(len=1), intent(in) :: jobz, range, uplo
         integer, intent(in)          :: n, lda, il, iu, ldz, lwork, lrwork, liwork
         complex(dp), intent(inout)   :: a(lda,*)
         complex(dp), intent(out)     :: z(ldz,*), work(*)
         real(dp), intent(in)         :: vl, vu, abstol
         real(dp), intent(out)        :: w(*), rwork(*)
         integer, intent(out)         :: m, info
         integer, intent(out)         :: isuppz(*), iwork(*)
      end subroutine zheevr

      subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, &
         & info)
         import :: dp
         character(len=1), intent(in) :: jobvl, jobvr
         integer, intent(in)          :: n, lda, ldvl, ldvr, lwork
         real(dp), intent(inout)      :: a(lda,*)
         real(dp), intent(out)        :: wr(*), wi(*), vl(ldvl,*), vr(ldvr,*), work(*)
         integer, intent(out)         :: info
      end subroutine dgeev

      subroutine zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, &
         & info)
         import :: dp
         character(len=1), intent(in) :: jobvl, jobvr
         integer, intent(in)          :: n, lda, ldvl, ldvr, lwork
         complex(dp), intent(inout)   :: a(lda,*)
         complex(dp), intent(out)     :: w(*), vl(ldvl,*), vr(ldvr,*), work(*)
         real(dp), intent(out)        :: rwork(*)
         integer, intent(out)         :: info
      end subroutine zgeev

      subroutine dgetrf(m, n, a, lda, ipiv, info)
         import :: dp
         integer, intent(in)          :: m, n, lda
         real(dp), intent(inout)      :: a(lda,*)
         integer, intent(out)         :: ipiv(*), info
      end subroutine dgetrf

      subroutine zgetrf(m, n, a, lda, ipiv, info)
         import :: dp
         integer, intent(in)          :: m, n, lda
         complex(dp), intent(inout)   :: a(lda,*)
         integer, intent(out)         :: ipiv(*), info
      end subroutine zgetrf

      subroutine dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp
         character(len=1), intent(in) :: trans
         integer, intent(in)          :: n, nrhs, lda, ldb
         real(dp), intent(in)         :: a(lda,*)
         integer, intent(in)          :: ipiv(*)
         real(dp), intent(inout)      :: b(ldb,*)
         integer, intent(out)         :: info
      end subroutine dgetrs

      subroutine zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp
         character(len=1), intent(in) :: trans
         integer, intent(in)          :: n, nrhs, lda, ldb
         complex(dp), intent(in)      :: a(lda,*)
         integer, intent(in)          :: ipiv(*)
         complex(dp), intent(inout)   :: b(ldb,*)
         integer, intent(out)         :: info
      end subroutine zgetrs

      subroutine dsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
         import :: dp
         character(len=1), intent(in) :: uplo
         integer, intent(in)          :: n, lda, lwork
         real(dp), intent(inout)      :: a(lda,*), work(*)
         integer, intent(out)         :: ipiv(*), info
      end subroutine dsytrf

      subroutine zhetrf(uplo, n, a, lda, ipiv, work, lwork, info)
         import :: dp
         character(len=1), intent(in) :: uplo
         integer, intent(in)          :: n, lda, lwork
         complex(dp), intent(inout)   :: a(lda,*), work(*)
         integer, intent(out)         :: ipiv(*), info
      end subroutine zhetrf

      subroutine zsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
         import :: dp
         character(len=1), intent(in) :: uplo
         integer, intent(in)          :: n, lda, lwork
         complex(dp), intent(inout)   :: a(lda,*), work(*)
         integer, intent(out)         :: ipiv(*), info
      end subroutine zsytrf

      subroutine dsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp
         character(len=1), intent(in) :: uplo
         integer, intent(in)          :: n, nrhs, lda, ldb
         real(dp), intent(in)         :: a(lda,*)
         integer, intent(in)          :: ipiv(*)
         real(dp), intent(inout)      :: b(ldb,*)
         integer, intent(out)         :: info
      end subroutine dsytrs

      subroutine zhetrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp
         character(len=1), intent(in) :: uplo
         integer, intent(in)          :: n, nrhs, lda, ldb
         complex(dp), intent(in)      :: a(lda,*)
         integer, intent(in)          :: ipiv(*)
         complex(dp), intent(inout)   :: b(ldb,*)
         integer, intent(out)         :: info
      end subroutine zhetrs

      subroutine zsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp
         character(len=1), intent(in) :: uplo
         integer, intent(in)          :: n, nrhs, lda, ldb
         complex(dp), intent(in)      :: a(lda,*)
         integer, intent(in)          :: ipiv(*)
         complex(dp), intent(inout)   :: b(ldb,*)
         integer, intent(out)         :: info
      end subroutine zsytrs

      subroutine dsaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, &
         & workd, workl, lworkl, info)
         import :: dp
         integer, intent(inout)       :: ido
         character(len=1), intent(in) :: bmat
         integer, intent(in)          :: n, nev, ncv, ldv, lworkl
         character(len=2), intent(in) :: which
         real(dp), intent(inout)      :: tol, resid(*), v(ldv,*), workd(*), workl(*)
         integer, intent(inout)       :: iparam(*), ipntr(*), info
      end subroutine dsaupd

      subroutine dseupd(rvec, howmny, select, d, z, ldz, sigma, bmat, n, which, nev, tol, &
         & resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info)
         import :: dp
         logical, intent(in)          :: rvec
         character(len=1), intent(in) :: howmny, bmat
         logical, intent(inout)       :: select(*)
         integer, intent(in)          :: ldz, n, nev, ncv, ldv, lworkl
         real(dp), intent(inout)      :: d(*), z(ldz,*), sigma, tol, resid(*), v(ldv,*), &
         & workd(*), workl(*)
         character(len=2), intent(in) :: which
         integer, intent(inout)       :: iparam(*), ipntr(*), info
      end subroutine dseupd

      subroutine dnaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, &
         & workd, workl, lworkl, info)
         import :: dp
         integer, intent(inout)       :: ido
         character(len=1), intent(in) :: bmat
         integer, intent(in)          :: n, nev, ncv, ldv, lworkl
         character(len=2), intent(in) :: which
         real(dp), intent(inout)      :: tol, resid(*), v(ldv,*), workd(*), workl(*)
         integer, intent(inout)       :: iparam(*), ipntr(*), info
      end subroutine dnaupd

      subroutine dneupd(rvec, howmny, select, dr, di, z, ldz, sigmar, sigmai, workev, bmat, &
         & n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, &
         & lworkl, info)
         import :: dp
         logical, intent(in)          :: rvec
         character(len=1), intent(in) :: howmny, bmat
         logical, intent(inout)       :: select(*)
         integer, intent(in)          :: ldz, n, nev, ncv, ldv, lworkl
         real(dp), intent(inout)      :: dr(*), di(*), z(ldz,*), sigmar, sigmai, workev(*), &
                                    &    tol, resid(*), v(ldv,*), workd(*), workl(*)
         character(len=2), intent(in) :: which
         integer, intent(inout)       :: iparam(*), ipntr(*), info
      end subroutine dneupd

      subroutine znaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, &
         & workd, workl, lworkl, rwork, info)
         import :: dp
         integer, intent(inout)       :: ido
         character(len=1), intent(in) :: bmat
         integer, intent(in)          :: n, nev, ncv, ldv, lworkl
         character(len=2), intent(in) :: which
         real(dp), intent(inout)      :: tol, rwork(*)
         complex(dp), intent(inout)   :: resid(*), v(ldv,*), workd(*), workl(*)
         integer, intent(inout)       :: iparam(*), ipntr(*), info
      end subroutine znaupd
        
      subroutine zneupd(rvec, howmny, select, d, z, ldz, sigma, workev, bmat, n, which, nev, &
         & tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, &
         & info)
         import :: dp
         logical, intent(in)          :: rvec
         character(len=1), intent(in) :: howmny, bmat
         logical, intent(inout)       :: select(*)
         integer, intent(in)          :: ldz, n, nev, ncv, ldv, lworkl
         complex(dp), intent(inout)   :: d(*), z(ldz,*), sigma, workev(*), resid(*), v(ldv,*), &
                                    &    workd(*), workl(*)
         real(dp), intent(inout)      :: tol, rwork(*)
         character(len=2), intent(in) :: which
         integer, intent(inout)       :: iparam(*), ipntr(*), info
      end subroutine zneupd

   end interface

   contains

   !================================================================
   ! The main subroutine for solving the Schrodinger equation
   !================================================================

   subroutine qsolver(rmin, rmax, &  ! integration limits
      &     bc_left, bc_right,    &  ! boundary conditions
      &     mu,                   &  ! effective mass mu = m/h^2
      &     JRot,                 &  ! rotation quantum number
      &     nEq,                  &  ! number of equations 
      &     N,                    &  ! number of mesh points
      &     V, V1, Vc, V1c,       &  ! matrix of potentials
      &     B, B1, Bc, B1c,       &  ! matrix of B-coupling functions
      &     G, Gc, G1, G1c,       &  ! diagonal matrix G(r) in H*psi = E*G*psi
      &     E,                    &  ! array with calculated energies 
      &     W,                    &  ! array with calculated resonance widths
      &     nE, nEf,              &  ! number of eigenvalues requested and found
      &     rc, phi,              &  ! complex scaling parameters
      &     mtype, ye, beta,      &  ! mapping parameters 
      &     is_V_sym,             &  ! is V matrix symmetric/Hermitian
      &     is_B_sym,             &  ! is B in Hermitian form
      &     E0, W0,               &  ! find energies and widths closest to E0 and W0
      &     Emin, Emax,           &  ! energy-range bounds
      &     Wmax,                 &  ! eigenvalues with W > Wmax will be ignored
      &     imin, imax,           &  ! index-range bounds
      &     sortE, sortW,         &  ! sorting modes
      &     WF, WF1, WFc, WF1c,   &  ! wavefunctions
      &     braW, braW1,          &  ! weighted bra vectors for real wavefunctions
      &     braWc, braW1c,        &  ! weighted bra vectors
      &     rv, rvc,              &  ! real and complex mesh points
      &     qw, qwc,              &  ! real and complex quadrature weights
      &     dV, dVc, dV1, dV1c,   &  ! potential derivatives
      &     dB, dBc, dB1, dB1c,   &  ! B-coupling derivatives
      &     dG, dGc, dG1, dG1c,   &  ! G-function derivatives
      &     ndE,                  &  ! number of derivative parameters
      &     dE, dW,               &  ! calculated energy and width derivatives 
      &     f_thr,                &  ! eigenvalue filter threshold
      &     err_mode,             &  ! error mode
      &     err_code,             &  ! error code
      &     err_message           &  ! error message
      &                           )

      implicit none

      !================================================================
      ! Scalar arguments
      !================================================================
      integer, intent(in)            :: N       
      integer, optional, intent(in)  :: JRot    
      integer, optional, intent(in)  :: nEq     
      integer, optional, intent(in)  :: nE
      integer, optional, intent(out) :: nEf
      integer, optional, intent(in)  :: imin, imax
      integer, optional, intent(in)  :: ndE

      logical, optional, intent(in)  :: is_V_sym  
      logical, optional, intent(in)  :: is_B_sym  
      integer, optional, intent(in)  :: sortE, sortW
      integer, optional, intent(in)  :: bc_left, bc_right

      real(dp), intent(in)           :: mu
      real(dp), optional, intent(in) :: rmin, rmax

      real(dp), optional, intent(in) :: E0, W0

      real(dp), optional, intent(in) :: Emin, Emax, Wmax

      real(dp), optional, intent(in) :: phi, rc
      real(dp), optional, intent(in) :: f_thr

      real(dp), optional, intent(in)     :: beta, ye
      character(*), optional, intent(in) :: mtype

      integer,      optional, intent(in)               :: err_mode
      integer,      optional, intent(out)              :: err_code
      character(:), optional, intent(out), allocatable :: err_message

      !================================================================
      ! Arrays arguments
      !================================================================
      real(dp),    optional, dimension(:),     intent(out) :: E, W
      real(dp),    optional, dimension(:),     intent(out) :: rv, qw
      real(dp),    optional, dimension(:,:),   intent(out) :: dE, dW
      real(dp),    optional, dimension(:,:,:), intent(out) :: WF
      real(dp),    optional, dimension(:,:),   intent(out) :: WF1
      real(dp),    optional, dimension(:,:,:), intent(out) :: braW
      real(dp),    optional, dimension(:,:),   intent(out) :: braW1
      complex(dp), optional, dimension(:),     intent(out) :: rvc, qwc
      complex(dp), optional, dimension(:,:,:), intent(out) :: WFc, braWc
      complex(dp), optional, dimension(:,:),   intent(out) :: WF1c, braW1c

      !================================================================
      ! Functions
      !================================================================
      ! one parameter functions f(x)
      procedure( func1), optional :: V1, B1, G1     ! real(dp)
      procedure(cfunc1), optional :: V1c, B1c, G1c  ! complex(dp)

      ! two-parameter diagonal matrix functions f(x, i)
      procedure( func2), optional :: G         ! real(dp)
      procedure(cfunc2), optional :: Gc        ! complex(dp)

      ! three parameter functions f(x, i, j)
      procedure( func3), optional :: V,  B     ! real(dp)
      procedure(cfunc3), optional :: Vc, Bc    ! complex(dp)

      ! derivative subroutines f(x, df/dp)
      procedure( vfunc1), optional :: dV1,  dB1,  dG1  ! real(dp)
      procedure(cvfunc1), optional :: dV1c, dB1c, dG1c ! complex(dp)

      ! derivative subroutines f(x, i, df/dp)
      procedure( vfunc2), optional :: dG   ! real(dp)
      procedure(cvfunc2), optional :: dGc  ! complex(dp)

      ! derivative subroutines f(x, i, j, df/dp)
      procedure( vfunc3), optional :: dV,  dB    ! real(dp)
      procedure(cvfunc3), optional :: dVc, dBc   ! complex(dp)

      !================================================================
      ! Local scalars
      !================================================================
     
      ! object that holds mapping parameters
      type(mapping_type), target ::  mapping

      ! objects those wrap input functions
      type(func_type)  :: Vf, dVf, Bf, dBf, Gf, dGf
      type(cfunc_type) :: Vcf, dVcf, Bcf, dBcf, Gcf, dGcf

      ! error handling (thread-local for this call)
      type(err_type), target :: err

      real(dp) :: inv2mu           ! 1/(2*mu)
      integer  :: l_nEq            ! local number of equations
      integer  :: l_ndE            ! local number of energy derivatives
      integer  :: l_nE             ! local number of eigenvalues to compute
      integer  :: l_nEf            ! local number of eigenvalues found
      integer  :: l_bc_left        ! local left boundary condition
      integer  :: l_bc_right       ! local right boundary condition
      integer  :: nH               ! the order of discretized Hamiltonian matrix
      real(dp) :: xmin, xmax       ! boundary points in transformed coordinates
      real(dp) :: l_Emin, l_Emax   ! local energy range
      integer  :: l_imin, l_imax   ! local index range
      real(dp) :: l_Wmax           ! local maximal resonance width value
      real(dp) :: yc, xc           ! complex radius 
      integer  :: i, j, k, m, p, q ! indices

      complex(dp) :: Ec0
      complex(dp) :: drdyc

      ! number of grid points on each subintervals
      integer     :: N1, N2

      ! diagonalization method arpack or lapack
      logical     :: is_arpack

      ! local eigenvalue filter threshold
      real(dp)    :: l_f_thr

      ! temporary variables
      real(dp)    :: t0, t1, t2
      complex(dp) :: c0, c1

      ! error status
      integer     :: istat

      logical :: is_efilter          ! is eigenvalue filtering enabled?
      logical :: is_quasi            ! is complex scaling used?
      logical :: is_ecs              ! is exterior complex scaling used?
      logical :: is_complex          ! is problem real or complex?
      logical :: is_complex_extrapol ! is extrapolation to the complex plane used?
      logical :: is_hermitian        ! is Hamiltonian Hermitian?
      logical :: is_V                ! is V present?
      logical :: is_B                ! is B coupling present?
      logical :: is_G                ! is diagonal matrix G present?
      logical :: is_G_complex        ! is diagonal matrix G complex-valued?
      logical :: is_B_sym_l          ! is B present in Hermitian form?
      logical :: is_V_sym_l          ! is V Hermitian (symmetric)?
      logical :: is_D_sym            ! do D matrices have symmetry?
      logical :: is_WF               ! is WF argument present?
      logical :: is_Centr            ! is centrifugal term present?
      logical :: is_lWF              ! is lWF argument present?
      logical :: is_complex_WF       ! is wavefunction complex?
      logical :: is_dE               ! is dE or dW argument present?
      logical :: is_dV               ! is dV argument present?
      logical :: is_dB               ! is dB argument present?
      logical :: is_dG               ! is dG argument present?
      logical :: need_WF             ! is WF need to be computed?
      logical :: need_lWF            ! is lWF need to be computed?
      logical :: is_erange           ! is erange
      logical :: is_irange           ! is index range requested
      logical :: is_irange_lapack    ! use LAPACK index range directly
      logical :: use_target_sort     ! use target-distance sorting
      integer :: l_sortE_mode        ! sort mode for energy (0/1/2)
      integer :: l_sortW_mode        ! sort mode for width (0/1/2)

      !================================================================
      ! Local arrays
      !================================================================
      ! vectors with calculated energies
      complex(dp), dimension(:), allocatable :: Ecv, Etcv
      real(dp),    dimension(:), allocatable :: Ev

      ! eigenvalue sorting indices
      integer,    dimension(:), allocatable :: Eindv, Etindv

      ! is Ev(i) true eigenvalue array
      logical,    dimension(:), allocatable :: is_Ev

      ! vectors contained grid points
      type(grid_type) :: xv             ! physical mesh as a grid object (GRID_PLAIN)
      real(dp),    dimension(:), allocatable :: yv
      type(grid_type), target :: grid_inner, grid_outer   ! ECS inner/outer grids (nodes + weights)
      complex(dp), dimension(:), allocatable :: l_rcv

      ! vectors containing function values used in the mapping
      real(dp),    dimension(:), allocatable :: Fv,  Qv
      real(dp),    dimension(:), allocatable :: gv,  ginvv
      complex(dp), dimension(:), allocatable :: Fcv, Qcv

      ! diagonal matrix G(r) on the grid and its inverse square root
      real(dp),    dimension(:,:), allocatable :: Gm
      complex(dp), dimension(:,:), allocatable :: Gcm
      real(dp),    dimension(:,:), allocatable :: Gsqrtim
      complex(dp), dimension(:,:), allocatable :: Gsqrticm

      type(grid_type), target :: g0   ! non-ECS spectral grid (nodes + weights)

      ! symmetrization vector for symmetric differential matrices
      real(dp),    dimension(:), allocatable :: bary

      ! differentiation matrices
      complex(dp), dimension(:,:), allocatable :: D1cm, D2cm
      real(dp),    dimension(:,:), allocatable :: D1m, D2m

      ! real and complex Hamiltonian matrix
      real(dp),    dimension(:,:), allocatable :: Hm
      complex(dp), dimension(:,:), allocatable :: Hcm

      ! eigenvector array
      real(dp),    dimension(:,:), allocatable :: Wm
      complex(dp), dimension(:,:), allocatable :: Wlcm, Wrcm


      ! arrays with potentials and coupling functions
      real(dp),    dimension(:,:,:), allocatable :: Vm, Bm
      complex(dp), dimension(:,:,:), allocatable :: Vcm, Bcm 

      ! arrays with derivatives of potentials and coupling functions
      real(dp),    dimension(:,:,:,:), allocatable :: dVm, dBm
      complex(dp), dimension(:,:,:,:), allocatable :: dVcm, dBcm
      real(dp),    dimension(:,:,:),   allocatable :: dGm
      complex(dp), dimension(:,:,:),   allocatable :: dGcm

      !================================================================
      ! Validate the arguments and set up the run: initialize the error
      ! handler and the coordinate mapping, check input consistency, and
      ! derive all the control flags, sizes and parameters (is_V/is_B/...,
      ! nH, boundary conditions, eigenvalue range, complex-scaling data).
      !================================================================
      call validate_inputs()
      if (err%has_error()) goto 100

      !=========================================================
      ! Now let's calculate!
      !=========================================================

      !=========================================================
      ! Build the collocation grids and differentiation operators, scale them
      ! onto the physical interval(s), assemble the physical mesh xv and the
      ! working first/second-derivative matrices (D1cm/D2cm or D1m/D2m).
      !=========================================================
      call build_diff_matrices()
      if (err%has_error()) goto 100

      !=========================================================
      ! Evaluate on the grid the potentials and couplings V/B/G (and their
      ! parameter derivatives), and the mapping factors gv/ginvv/F/Q.
      !=========================================================
      call sample_grid_functions()
      if (err%has_error()) goto 100

      !=========================================================
      ! Assemble the Hamiltonian (Hm or Hcm) from the differentiation
      ! matrices, potentials and couplings, and reduce the generalized
      ! problem to a standard one via the G^{-1/2} transform.
      !=========================================================
      call build_hamiltonian()
      if (err%has_error()) goto 100

      !=========================================================
      ! Diagonalize the Hamiltonian (LAPACK or ARPACK) and recover the
      ! right (and, when needed, left) eigenvectors.
      !=========================================================
      call solve_eigenproblem()
      if (err%has_error()) goto 100

      if (l_nEf .eq. 0) then 
         if (present(nEf)) nEf = 0
         goto 100
      endif

      !=================================================
      ! Flag spurious eigenvalues. The keep/discard mask is_Ev is allocated
      ! and set here (all true); the optional spectral-decay test then clears
      ! the spurious entries. sort_eigenpairs reads is_Ev afterwards.
      !=================================================
      call filter_spurious_states()
      if (err%has_error()) goto 100

      !=================================================
      ! Sort the eigenvalues (with their index map Eindv) and drop spurious
      ! and out-of-range states, leaving l_nEf wanted pairs.
      !=================================================
      call sort_eigenpairs()
      if (err%has_error()) goto 100

      !=================================================
      ! Energy derivatives dE (and width derivatives dW) with respect to the
      ! model parameters; builds the derivative Hamiltonian dHm/dHcm.
      !=================================================
      call compute_energy_derivatives()
      if (err%has_error()) goto 100

      !=========================================================
      ! Fill the caller's optional outputs: energies, widths, wavefunctions,
      ! weighted bra vectors, grid points, quadrature weights and counts.
      !=========================================================
      call assemble_outputs()
      if (err%has_error()) goto 100

      !=================================================
      ! That's all, folks!
      !=================================================

      !=================================================
      ! But one thing remains to do. It is error handling
      !=================================================
      100 continue

      call err%finalize(err_code, err_message)

   contains

      !=================================================
      ! Internal procedure (see the call site for what this phase does);
      ! shares qsolver state by host association.
      !=================================================
      subroutine validate_inputs()
         if (present(nEf)) nEf = 0

         !================================================================
         ! 1. Error handler and helper-object wiring
         !================================================================
         if (present(err_mode)) then
            if (err_mode < ERRMODE_PRINTALL .or. err_mode > ERRMODE_STOP) then
               call err%init(ERRMODE_PRINTALL)
               call err%push_err(ERR_INPUT, 'err_mode must be 1, 2, or 3')
               return
            endif
            call err%init(err_mode)
         else
            call err%init(ERRMODE_PRINTALL)
         end if

         ! attach error handler to helper objects
         mapping%err => err
         Vf%err      => err
         dVf%err     => err
         Bf%err      => err
         dBf%err     => err
         Gf%err      => err
         dGf%err     => err
         Vcf%err     => err
         dVcf%err    => err
         Bcf%err     => err
         dBcf%err    => err
         Gcf%err     => err
         dGcf%err    => err

         !================================================================
         ! 2. Coordinate mapping and integration domain [xmin, xmax]
         !================================================================
         if (.not. present(mtype)) then 
            ! no mapping
            call mapping%init()

         else if (present(beta) .and. present(ye)) then
            ! two-parametric mapping
            call mapping%init(mtype, beta, ye)   

         else if (present(beta)) then
            ! one-parametric mapping
            call mapping%init(mtype, beta)   

         else 
            ! no mapping
            call mapping%init(mtype)   
         endif

         if (err%has_error()) return

         if (present(rmin)) then 
            xmin = mapping%xy(rmin)
         else
            xmin = mapping%xyMInf() 
         endif

         if (present(rmax)) then 
            xmax = mapping%xy(rmax)
         else
            xmax = mapping%xyInf() 
         endif

         if (err%has_error()) return

         if (xmax .le. xmin) then
            call err%push_err(ERR_INPUT, 'rmax must be greater than rmin')
            return
         endif

         !================================================================
         ! 3. Boundary conditions
         !================================================================
         l_bc_left = BC_DIRICHLET
         if (present(bc_left)) l_bc_left = bc_left

         l_bc_right = BC_DIRICHLET
         if (present(bc_right)) l_bc_right = bc_right

         if (l_bc_left /= BC_DIRICHLET .and. l_bc_left /= BC_REGULAR) then
            call err%push_err(ERR_INPUT, 'bc_left has wrong value')
            return
         endif

         if (l_bc_right /= BC_DIRICHLET .and. l_bc_right /= BC_REGULAR) then
            call err%push_err(ERR_INPUT, 'bc_right has wrong value')
            return
         endif

         if (l_bc_left == BC_REGULAR .and. .not. present(rmin)) then
            call err%push_err(ERR_INPUT, 'bc_left = BC_REGULAR requires finite rmin')
            return
         endif

         if (l_bc_right == BC_REGULAR .and. .not. present(rmax)) then
            call err%push_err(ERR_INPUT, 'bc_right = BC_REGULAR requires finite rmax')
            return
         endif

         !================================================================
         ! 4. Problem dimensions and Schrodinger parameters
         !================================================================
         if (abs(mu) <= 0.5_dp / huge(mu)) then
            call err%push_err(ERR_INPUT, 'mu is too small: 1/(2*mu) would overflow')
            return
         end if

         inv2mu = 1/(2*mu)

         if (present(nEq)) then
            if (nEq < 1) then
               call err%push_err(ERR_INPUT, 'Number of equations must be positive')
               return
            endif
            l_nEq = nEq
         else
            l_nEq = 1
         endif

         if (N < 1) then
            call err%push_err(ERR_INPUT, 'N must be positive')
            return
         endif

         nH = N*l_nEq ! The order of the Hamiltonian matrix

         if (present(ndE)) then
            if (ndE < 0) then
               call err%push_err(ERR_INPUT, 'ndE must not be negative')
               return
            endif
            l_ndE = ndE
         else
            l_ndE = 1
         endif

         is_dE = (present(dE) .or. present(dW)) .and. (l_ndE > 0)

         !================================================================
         ! 5. Run mode: centrifugal, resonance (complex scaling), filter
         !================================================================
         if (present(JRot)) then
            if (JRot .ne. 0) then 
               is_Centr = .true. 
            else
               is_Centr = .false.
            endif
         else
            is_Centr = .false.
         endif

         is_quasi = present(phi)
      
         if (present(rc) .and. .not. is_quasi) then
            call err%push_err(ERR_INPUT, 'rc can be used only together with phi (complex scaling)')
            return
         endif

         is_ecs = present(rc)

         if (is_ecs .and. l_bc_right /= BC_DIRICHLET) then
            call err%push_err(ERR_INPUT, &
               &              'Complex scaling requires Dirichlet boundary at the right endpoint')
            return
         endif

         ! complex-scaling data: rotation factor and junction point xc/yc
         if (is_quasi) then
            j  = count((/ present(V), present(V1), present(B), present(B1), present(G), present(G1), &
               &           present(dV), present(dV1), present(dB), present(dB1), present(dG), present(dG1) /))
            if (j > 0) then 
               is_complex_extrapol = .true.
            else
               is_complex_extrapol = .false.
            endif

            drdyc = cmplx(cos(phi), sin(phi), dp)
         else
            drdyc = cmplx(1.0_dp, 0.0_dp, dp)
         endif

         xc = xmin

         if (is_ecs) then 
            yc = rc
            xc = mapping%xy(rc)
            if (xc <= xmin) then
               call err%push_err(ERR_INPUT, 'Parameter rc must be greater than rmin')
               return
            endif
            if (xc >= xmax) then
               call err%push_err(ERR_INPUT, 'Parameter rc must be less than rmax')
               return
            endif
         elseif (.not. present(rmin) .and. is_quasi) then
            call err%push_err(ERR_INPUT, 'If rmin = -infinity then rc must be specified if complex scaling used')
            return
         endif

         if (is_quasi .or. present(f_thr)) then
            is_efilter = .true. 
            if (present(f_thr)) then
               l_f_thr = f_thr
               if (f_thr < 0.0_dp) is_efilter = .false.
            else
               l_f_thr = 1.0e-6_dp ! default value
            endif
         else
            is_efilter = .false. 
         endif

         !================================================================
         ! 6. Input functions: selection, consistency and wrapper binding
         !================================================================
         j = count((/ present(V1), present(V1c), present(B1), present(B1c), present(G1), present(G1c)/))

         if (j > 0 .and. l_nEq > 1) then
            call err%push_err(ERR_INPUT, 'Number of equations is inconsistent with number of functions')
            return
         endif

         j = count((/ present(V), present(Vc), present(V1), present(V1c)/))

         if (j > 1) then
            call err%push_err(ERR_INPUT, 'Cannot decide which potential to use')
            return
         endif

         is_V = (j == 1) 

         j = count((/ present(B), present(Bc), present(B1), present(B1c)/))

         if (j > 1) then
            call err%push_err(ERR_INPUT, 'Cannot decide which B-function to use')
            return
         endif

         is_B = (j == 1)

         j = count((/ present(G), present(Gc), present(G1), present(G1c) /))

         if (j > 1) then
            call err%push_err(ERR_INPUT, 'Cannot decide which G-function to use')
            return
         endif

         is_G = (j == 1)
         is_G_complex = present(Gc) .or. present(G1c)

         j = count((/ present(dV), present(dVc), present(dV1), present(dV1c)/))

         if (is_dE .and. present(dE) .and. is_V .and. j .eq. 0) then
            call err%push_err(ERR_INPUT, 'Energy derivatives cannot be computed without providing &
            & derivatives of V')
            return
         endif

         if (is_dE .and. present(dW) .and. is_V .and. j .eq. 0) then
            call err%push_err(ERR_INPUT, 'Resonance widths derivatives cannot be computed without providing &
            & derivatives of V')
            return
         endif

         if (is_dE .and. j > 1) then
            call err%push_err(ERR_INPUT, 'Cannot decide which dV-function to use')
            return
         endif

         is_dV = is_dE .and. (j == 1) 

         j = count((/ present(dB), present(dBc), present(dB1), present(dB1c)/))

         if (is_dE .and. present(dE) .and. is_B .and. j .eq. 0) then
            call err%push_err(ERR_INPUT, 'Energy derivatives cannot be computed without providing &
            & derivatives of B-coupling functions')
            return
         endif

         if (is_dE .and. present(dW) .and. is_B .and. j .eq. 0) then
            call err%push_err(ERR_INPUT, 'Resonance widths derivatives cannot be computed without providing &
            & derivatives of B-coupling functions')
            return
         endif

         if (is_dE .and. j > 1) then
            call err%push_err(ERR_INPUT, 'Cannot decide which dB-function to use')
            return
         endif

         is_dB = is_dE .and. (j == 1)  

         j = count((/ present(dG), present(dGc), present(dG1), present(dG1c)/))

         if (is_dE .and. j > 1) then
            call err%push_err(ERR_INPUT, 'Cannot decide which dG-function to use')
            return
         endif

         is_dG = is_dE .and. (j == 1)

         j = count((/ present(Vc), present(V1c), present(Bc), present(B1c), present(Gc), present(G1c) /))

         is_complex = (j > 0) .or. is_quasi

         if (.not. is_complex) then
            if (is_dE .and. (present(dVc) .or. present(dV1c))) then
               call err%push_warning('Complex dV-function provided for a real problem; only the real part will be used')
            endif
            if (is_dE .and. (present(dBc) .or. present(dB1c))) then
               call err%push_warning('Complex dB-function provided for a real problem; only the real part will be used')
            endif
            if (is_dE .and. (present(dGc) .or. present(dG1c))) then
               call err%push_warning('Complex dG-function provided for a real problem; only the real part will be used')
            endif
         endif

         ! bind the chosen functions to their wrappers
         ! real functions
         if (present(V1))  call Vf%init_func1(V1)
         if (present(V))   call Vf%init_func3(V)

         if (present(B1))  call Bf%init_func1(B1)
         if (present(B))   call Bf%init_func3(B)

         if (present(G1))  call Gf%init_func1(G1)
         if (present(G))   call Gf%init_func2(G)

         if (is_dE .and. present(dV1)) call dVf%init_vfunc1(dV1)
         if (is_dE .and. present(dV))  call dVf%init_vfunc3(dV)
         dVf%n_vec = l_ndE

         if (is_dE .and. present(dB1)) call dBf%init_vfunc1(dB1)
         if (is_dE .and. present(dB))  call dBf%init_vfunc3(dB)
         dBf%n_vec = l_ndE

         if (is_dE .and. present(dG1)) call dGf%init_vfunc1(dG1)
         if (is_dE .and. present(dG))  call dGf%init_vfunc2(dG)
         dGf%n_vec = l_ndE

         ! complex ones
         if (present(V1c))  call Vcf%init_cfunc1(V1c)
         if (present(Vc))   call Vcf%init_cfunc3(Vc)
                         
         if (present(B1c))  call Bcf%init_cfunc1(B1c)
         if (present(Bc))   call Bcf%init_cfunc3(Bc)

         if (present(G1c))  call Gcf%init_cfunc1(G1c)
         if (present(Gc))   call Gcf%init_cfunc2(Gc)

         if (is_dE .and. present(dV1c)) call dVcf%init_cvfunc1(dV1c)
         if (is_dE .and. present(dVc))  call dVcf%init_cvfunc3(dVc)
         dVcf%n_vec = l_ndE

         if (is_dE .and. present(dB1c)) call dBcf%init_cvfunc1(dB1c)
         if (is_dE .and. present(dBc))  call dBcf%init_cvfunc3(dBc)
         dBcf%n_vec = l_ndE

         if (is_dE .and. present(dG1c)) call dGcf%init_cvfunc1(dG1c)
         if (is_dE .and. present(dGc))  call dGcf%init_cvfunc2(dGc)
         dGcf%n_vec = l_ndE

         !================================================================
         ! 7. Matrix character: symmetry and Hermiticity
         !================================================================
         ! defaults V symmetry and B in symmetric form
         is_B_sym_l = .true.
         is_V_sym_l = .true.

         ! override if provided
         if (present(is_B_sym)) is_B_sym_l = is_B_sym
         if (present(is_V_sym)) is_V_sym_l = is_V_sym

         ! quasi case overrides symmetry if complex functions are provided
         if (is_quasi) then 
            if (present(Vc)) is_V_sym_l = .false.
         endif

         is_hermitian = (.not. is_quasi) .and. is_B_sym_l .and. is_V_sym_l .and. (.not. is_G_complex)

         ! Regular endpoint matrices are not symmetrized yet; use the general
         ! eigenvalue path instead of Hermitian LAPACK for these boundary cases.
         if (l_bc_left == BC_REGULAR .or. l_bc_right == BC_REGULAR) is_hermitian = .false.

         !================================================================
         ! 8. Requested outputs (wavefunctions / bra vectors)
         !================================================================
         j = count((/ present(WF1), present(WF1c), present(braW1), present(braW1c)/))
         if (j > 0 .and. l_nEq > 1) then
            call err%push_err(ERR_INPUT, 'Single-channel wavefunction output is valid only for single-channel case')
            return
         endif

         is_WF  = present(WF) .or. present(WF1) .or. present(WFc) .or. present(WF1c) .or. &
         &        present(braW) .or. present(braW1)

         is_lWF = present(braWc) .or. present(braW1c)

         need_WF  = is_WF .or. is_dE .or. is_efilter .or. is_lWF 
         need_lWF = (.not. is_hermitian) .and. (is_lWF .or. is_dE) 

         is_complex_WF = is_complex .or. (.not. is_hermitian)

         if (present(braW) .and. is_complex_WF) then
            call err%push_err(ERR_INPUT, 'braW is real; use braWc for complex or non-Hermitian problems')
            return
         endif
         if (present(braW1) .and. is_complex_WF) then
            call err%push_err(ERR_INPUT, 'braW1 is real; use braW1c for complex or non-Hermitian problems')
            return
         endif

         !================================================================
         ! 9. Eigensolver configuration: sorting, count, ranges, backend
         !================================================================
         if (present(sortE)) then
            l_sortE_mode = sortE
         else
            if (present(E0)) then
               l_sortE_mode = 0
            else
               l_sortE_mode = 1
            endif
         endif

         if (present(sortW)) then
            l_sortW_mode = sortW
         else
            l_sortW_mode = 0
         endif

         if (.not. present(sortE) .and. present(sortW)) then
            l_sortE_mode = 0
         endif

         if (l_sortE_mode < 0 .or. l_sortE_mode > 2) then
            call err%push_err(ERR_INPUT, 'sortE must be 0, 1, or 2')
            return
         endif

         if (l_sortW_mode < 0 .or. l_sortW_mode > 2) then
            call err%push_err(ERR_INPUT, 'sortW must be 0, 1, or 2')
            return
         endif

         if (is_hermitian) then
            if (l_sortW_mode /= 0) then
               call err%push_warning('sortW is ignored for Hermitian problems')
            endif
            l_sortW_mode = 0
            use_target_sort = present(E0) .and. (.not. present(sortE))
         else
            if (l_sortE_mode /= 0 .and. l_sortW_mode /= 0) then
               call err%push_err(ERR_INPUT, 'sortE and sortW cannot both be nonzero')
               return
            endif
            use_target_sort = present(E0) .and. (.not. present(sortE)) .and. (.not. present(sortW))
         endif

         ! setting number of requested eigenvalues
         if (present(nE)) then 
            if (nE > 0) then
               l_nE = nE
            elseif (nE .eq. 0) then
               if (present(nEf)) nEf = 0
               return
            else 
               l_nE = nH
            endif   
         elseif (present(E0)) then
            l_nE = 1
         else
            l_nE = nH
         endif

         ! Initialisation of eigenvalue range
         is_irange = present(imax) .or. present(imin)
         is_erange = present(Emax) .or. present(Emin)

         if (present(imin)) then 
            l_imin = imin
         else 
            l_imin = 1
         endif

         if (present(imax)) then 
            l_imax = imax
         else 
            l_imax = min(nH, l_imin + l_nE - 1)
         endif

         if (l_imin > l_imax) then
            call err%push_err(ERR_INPUT, 'imin must be less than or equal to imax')
            return
         endif

         if (is_irange) then
            if (l_imin < 1) then
               call err%push_err(ERR_INPUT, 'imin is out of range')
               return
            endif
            if (l_imax > nH) then
               call err%push_err(ERR_INPUT, 'imax is out of range')
               return
            endif
            ! Ensure enough eigenvalues are computed before applying
            ! the final [imin, imax] index slice.
            l_nE = max(l_nE, l_imax)
         endif

         is_irange_lapack = is_hermitian .and. (.not. present(E0)) .and. (.not. is_erange) .and. &
         &                  (is_irange .or. (l_nE > 0 .and. l_nE < nH))

         if (present(Emin)) then 
            l_Emin = Emin
         else
            l_Emin = -huge(1.0_dp)
         endif

         if (present(Emax)) then 
            l_Emax = Emax
         else
            l_Emax =  huge(1.0_dp)
         endif

         if (l_Emax <= l_Emin) then
            call err%push_err(ERR_INPUT, 'Emin must be less than Emax')
            return
         endif

         if (present(Wmax)) then 
            l_Wmax = Wmax
         else
            l_Wmax =  huge(1.0_dp)
         endif

         !=========================================================
         ! Use target ARPACK only when the corresponding wrapper accepts
         ! the requested nev; otherwise fall back to LAPACK target solve.
         !=========================================================
         if (present(E0)) then
            if ((.not. is_complex) .and. is_hermitian) then
               is_arpack = l_nE < nH
            else
               is_arpack = l_nE < nH - 1
            endif
         else
            is_arpack = .false.
         endif

         !================================================================
         ! 10. Target eigenvalue
         !================================================================
         if (present(E0)) then
            if (present(W0)) then
               Ec0 = cmplx(E0, -W0/2, dp)
            else
               Ec0 = cmplx(E0, 0.0_dp, dp)
            endif
         endif      
      end subroutine validate_inputs

      !=================================================
      ! Internal procedure (see the call site for what this phase does);
      ! shares qsolver state by host association.
      !=================================================
      subroutine build_diff_matrices()
         type(diffmat_type) :: dm0, dm_inner, dm_outer  ! non-ECS + ECS differentiation operators
         real(dp) :: qxc                                ! mapping correction Q(xc) at the ECS junction

         if (is_complex) then
            allocate(D1cm(1:N,1:N), D2cm(1:N,1:N), stat=istat)
         else
            allocate(D1m(1:N,1:N), D2m(1:N,1:N), stat=istat)
         endif
         if (alloc_failed(istat, 'Allocation failed in the derivative matrix creation block', err)) return

         !=========================================================
         !  Interval splitting if exterior complex scaling
         !=========================================================
         if (is_ecs) then
            N1 = int(N*(xc-xmin)/(xmax-xmin))
            N2 = N - N1 + 1
            if (N1 < 2) then
               call err%push_err(ERR_INPUT, 'There is no points inside region of real r. Increase rc or remove it')
               return
            endif
            if (N2 < 2) then
               call err%push_err(ERR_INPUT, 'There is no points inside region of complex r. Decrease rc')
               return
            endif
         endif

         ! Physical mesh xv: a plain grid holding the N active mesh points 1:N
         ! plus the formal boundaries at 0 and N+1. Allocated once here; the
         ! nodes are filled below from the (already scaled) spectral grids, the
         ! same way for the ECS and non-ECS cases.
         call xv%build(N + 1, GRID_PLAIN, err)
         if (err%has_error()) return

         if (.not. is_ecs) then
            !=========================================================
            ! Non-ECS: build the differentiation operator once on a single
            ! interval, dispatching the grid kind by the boundary conditions.
            !=========================================================
            ! Pick the collocation family from the (left, right) boundary-condition
            ! pair, encoded as one selector (BC_DIRICHLET = 0, BC_REGULAR = 1).
            select case (2*l_bc_left + l_bc_right)
            case (2*BC_DIRICHLET + BC_DIRICHLET)        ! both Dirichlet -> Gauss-Lobatto
               call g0%build(N + 1, GRID_LOBATTO, err)
               call dm0%build_sym(g0, err)
               bary = g0%bary
               is_D_sym = .true.
            case (2*BC_DIRICHLET + BC_REGULAR)          ! Dirichlet at left -> Gauss-Radau(-1)
               call g0%build(N + 1, GRID_RADAU_M1, err)
               call dm0%build(g0, err)
               is_D_sym = .false.
            case (2*BC_REGULAR + BC_DIRICHLET)          ! Dirichlet at right -> Gauss-Radau(+1)
               call g0%build(N + 1, GRID_RADAU_P1, err)
               call dm0%build(g0, err)
               is_D_sym = .false.
            case default                                ! neither -> Gauss-Legendre
               call g0%build(N + 1, GRID_LEGENDRE, err)
               call dm0%build(g0, err)
               is_D_sym = .false.
            end select
            if (err%has_error()) return
 
            !=========================================================
            ! Move the grid and the operator onto [xmin, xmax] in one call: the
            ! quadrature weights and D1/D2 pick up their scale factors here.
            !=========================================================
            call dm0%set_scale(xmin, xmax, err)
            if (err%has_error()) return

            ! g0 now sits on [xmin, xmax], so its nodes are the mesh points.
            xv%x = g0%x
         endif

         if (is_quasi) then

            if (is_ecs) then
               !=========================================================
               ! Creating first and second derivative matrices in case of exterior
               ! complex scaling. Integration interval is divided by two subintervals.
               !
               ! Left subinterval [xmin; xc]:
               !   - BC_DIRICHLET at xmin: Gauss-Lobatto grid; grid_inner%x(0) = xmin is
               !     a node that is dropped from the solution vector (Dirichlet).
               !   - BC_REGULAR  at xmin: Gauss-Radau grid (endpoint = +1); grid_inner%x(0)
               !     is only a formal boundary, not a node. All N1 active nodes
               !     grid_inner%x(1:N1) enter the solution vector.
               !
               ! Right subinterval [xc; xmax] uses Gauss-Lobatto in both cases;
               ! the right boundary is constrained to Dirichlet (checked above).
               !=========================================================
               if (l_bc_left == BC_DIRICHLET) then
                  call grid_inner%build(N1, GRID_LOBATTO, err)
               else
                  call grid_inner%build(N1, GRID_RADAU_P1, err)
               endif
               call dm_inner%build(grid_inner, err)
               call dm_inner%build_hermite(err)

               call grid_outer%build(N2, GRID_LOBATTO, err)
               call dm_outer%build(grid_outer, err)
               if (err%has_error()) return

               !=========================================================
               ! Move each subinterval operator (and its grid) onto the physical
               ! subinterval. set_scale folds the chain-rule factors into D1/D2
               ! and the Hermite matrices DH1/DH2, and rescales the quadrature
               ! weights, so the assembly below contains no explicit 1/s factors.
               !=========================================================
               call dm_inner%set_scale(xmin, xc,   err)
               call dm_outer%set_scale(xc,   xmax, err)
               if (err%has_error()) return

               !=========================================================
               ! Physical mesh: the two grids now sit on their subintervals, so
               ! the mesh points are read straight from the scaled grid nodes --
               ! inner active nodes 1:N1 (grid_inner%x(N1) = xc is the junction),
               ! outer interior nodes map to N1+1:N -- with the formal boundaries
               ! at 0 and N+1.
               !=========================================================
               xv%x(0)      = xmin
               xv%x(N+1)    = xmax
               xv%x(1:N1)   = grid_inner%x(1:N1)
               xv%x(N1+1:N) = grid_outer%x(1:N2-1)

               !=========================================================
               ! Outer-side physical derivative d/dr at the junction xc, expressed
               ! on the solution values: drdyc-rotated outer derivative plus the
               ! mapping correction Qx. dm_outer%D1(0, :) is the (already scaled)
               ! outer first-derivative row at the junction node.
               !=========================================================
               qxc = mapping%Qx(xc)   ! loop-invariant junction correction (hoisted)
               D1cm(1:N1, 1:N1) = dm_inner%DH1(1:N1, 1:N1)
               do i = 1, N1
                  D1cm(i, N1) = D1cm(i, N1)                              &
                              + dm_inner%DH1(i, N1+1) *                  &
                                 (conjg(drdyc) * dm_outer%D1(0, 0)       &
                              - qxc * (1 - conjg(drdyc)))

                  D1cm(i, N1+1:N) = conjg(drdyc)*dm_outer%D1(0, 1:N2-1)*dm_inner%DH1(i, N1+1)
               enddo

               D1cm(N1+1:N, N1:N) = dm_outer%D1(1:N2-1, 0:N2-1)*conjg(drdyc)

               !=========================================================
               ! D1cm lower-left corner block is the only piece that no
               ! assignment above touches; zero it explicitly.
               !=========================================================
               D1cm(N1+1:N, 1:N1-1) = (0.0_dp, 0.0_dp)

               if (l_bc_left == BC_DIRICHLET) then
                  !=========================================================
                  ! Simplified form, valid for the Gauss-Lobatto inner grid
                  ! (bc_left = BC_DIRICHLET).
                  !
                  ! Feature of the GL Hermite differentiation matrix: the
                  ! extra-column contribution dm_inner%DH2(i, N1+1) vanishes for
                  ! i < N1, and for those rows dm_inner%DH2(i, 1:N1) coincides with
                  ! dm_inner%D2(i, 1:N1). Only the matching row i = N1 needs the
                  ! Hermite-augmented entries.
                  !=========================================================
                  D2cm(1:N1-1, 1:N1) = dm_inner%D2(1:N1-1, 1:N1)

                  D2cm(N1, 1:N1) = dm_inner%DH2(N1, 1:N1)

                  D2cm(N1, N1) = D2cm(N1, N1)                               &
                               + dm_inner%DH2(N1, N1+1) *                   &
                                  (conjg(drdyc) * dm_outer%D1(0, 0)         &
                               - qxc * (1 - conjg(drdyc)))

                  D2cm(N1, N1+1:N) = conjg(drdyc)*dm_outer%D1(0, 1:N2-1)*dm_inner%DH2(N1, N1+1)

                  !=========================================================
                  ! Rows 1:N1-1, columns N1+1:N are not touched in this
                  ! branch (the GL simplification leaves them out); zero
                  ! them explicitly.
                  !=========================================================
                  D2cm(1:N1-1, N1+1:N) = (0.0_dp, 0.0_dp)
               else
                  !=========================================================
                  ! Universal form, used for the Gauss-Radau inner grid
                  ! (bc_left = BC_REGULAR). The GL-specific simplification
                  ! above does NOT hold here, so dm_inner%DH2(i, N1+1) must be
                  ! carried for every row i = 1:N1.
                  !=========================================================
                  D2cm(1:N1, 1:N1) = dm_inner%DH2(1:N1, 1:N1)

                  do i = 1, N1
                     D2cm(i, N1) = D2cm(i, N1)                              &
                                 + dm_inner%DH2(i, N1+1) *                  &
                                    (conjg(drdyc) * dm_outer%D1(0, 0)       &
                                 - qxc * (1 - conjg(drdyc)))

                     D2cm(i, N1+1:N) = conjg(drdyc)*dm_outer%D1(0, 1:N2-1)*dm_inner%DH2(i, N1+1)
                  enddo
               endif

               D2cm(N1+1:N, N1:N) = dm_outer%D2(1:N2-1, 0:N2-1)*conjg(drdyc)**2

               !=========================================================
               ! D2cm lower-left corner block (common to both branches) is
               ! not touched by any assignment above; zero it explicitly.
               !=========================================================
               D2cm(N1+1:N, 1:N1-1) = (0.0_dp, 0.0_dp)

               is_D_sym = .false.

            else
               !=========================================================
               ! Creating complex scaled first and second derivative matrices
               ! in case of single interval. dm0 was already moved onto
               ! [xmin, xmax] by set_scale, so only the complex-rotation
               ! factor drdyc remains here.
               !=========================================================

               D1cm(1:N,1:N) = conjg(drdyc)*dm0%D1(1:N,1:N)
               D2cm(1:N,1:N) = conjg(drdyc)**2*dm0%D2(1:N,1:N)
            endif
         else
            !=========================================================
            ! No complex coordinate scaling: dm0 was already moved onto
            ! [xmin, xmax] by set_scale, so its D1/D2 are copied straight into
            ! the working matrices -- complex D1cm/D2cm if the Hamiltonian is
            ! complex (e.g. a complex potential), otherwise real D1m/D2m.
            !=========================================================
            if (is_complex) then
               D1cm(1:N,1:N) = dm0%D1(1:N,1:N)
               D2cm(1:N,1:N) = dm0%D2(1:N,1:N)
            else
               D1m(1:N,1:N)  = dm0%D1(1:N,1:N)
               D2m(1:N,1:N)  = dm0%D2(1:N,1:N)
            endif

         endif
      end subroutine build_diff_matrices

      !=================================================
      ! Evaluate, on the physical grid, every per-node array the Hamiltonian
      ! needs: potentials and couplings V/B/G (and their parameter
      ! derivatives), plus the mapping factors gv, ginvv and F/Q. Internal
      ! procedure: shares qsolver's grid, mapping and arrays by host
      ! association.
      !=================================================
      subroutine sample_grid_functions()
         complex(dp), allocatable :: xcv(:)        ! complex grid points for extrapolation
         real(dp),    allocatable :: drv(:)        ! real derivative-callback values
         complex(dp), allocatable :: dcv(:)        ! complex derivative-callback values

         ! physical grid in the mapped y-coordinate (needed below to evaluate
         ! the potentials/couplings; also reused later by assemble_outputs)
         allocate(yv(1:N), stat=istat)
         if (alloc_failed(istat, 'Allocation failed for yv', err)) return
         do i = 1, N
            yv(i) = mapping%yx(xv%x(i))            ! grid in y-coordinates
         enddo

         if (is_complex) then
            !=========================================================
            ! Allocating complex matrices with functions values on grid
            !=========================================================
            if (is_V .or. is_Centr) then
               allocate(Vcm(1:N+1, l_nEq, l_nEq), stat=istat)
               if (alloc_failed(istat, 'Allocation failed for Vcm', err)) return
            endif
            if (is_B) then
               allocate(Bcm(1:N+1, l_nEq, l_nEq), stat=istat)
               if (alloc_failed(istat, 'Allocation failed for Bcm', err)) return
            endif
            if (is_G) then
               allocate(Gcm(1:N+1, l_nEq), stat=istat)
               if (alloc_failed(istat, 'Allocation failed for Gcm', err)) return
            endif
            if (is_dV) then
               allocate(dVcm(1:N+1, l_ndE, l_nEq, l_nEq), stat=istat)
               if (alloc_failed(istat, 'Allocation failed for dVcm', err)) return
            endif
            if (is_dB) then
               allocate(dBcm(1:N+1, l_ndE, l_nEq, l_nEq), stat=istat)
               if (alloc_failed(istat, 'Allocation failed for dBcm', err)) return
            endif
            if (is_dG) then
               allocate(dGcm(1:N+1, l_ndE, l_nEq), stat=istat)
               if (alloc_failed(istat, 'Allocation failed for dGcm', err)) return
            endif
         else
            !=========================================================
            ! Allocating real matrices with functions values on grid 
            !=========================================================
            if (is_V .or. is_Centr) then
               allocate(Vm(1:N+1,l_nEq,l_nEq), stat=istat)
               if (alloc_failed(istat, 'Allocation failed for Vm', err)) return
            endif
            if (is_B) then
               allocate(Bm(1:N+1,l_nEq,l_nEq), stat=istat)
               if (alloc_failed(istat, 'Allocation failed for Bm', err)) return
            endif
            if (is_G) then
               allocate(Gm(1:N+1, l_nEq), stat=istat)
               if (alloc_failed(istat, 'Allocation failed for Gm', err)) return
            endif
            if (is_dV) then
               allocate(dVm(1:N+1, l_ndE, l_nEq, l_nEq), stat=istat)
               if (alloc_failed(istat, 'Allocation failed for dVm', err)) return
            endif
            if (is_dB) then
               allocate(dBm(1:N+1, l_ndE, l_nEq, l_nEq), stat=istat)
               if (alloc_failed(istat, 'Allocation failed for dBm', err)) return
            endif
            if (is_dG) then
               allocate(dGm(1:N+1, l_ndE, l_nEq), stat=istat)
               if (alloc_failed(istat, 'Allocation failed for dGm', err)) return
            endif
         endif

         if (is_Centr .and. .not. is_V) then
            if (is_complex) then
               Vcm = (0.0_dp, 0.0_dp)
            else
               Vm = 0.0_dp
            end if
         endif

         if (is_dE) then
            allocate(drv(l_ndE), dcv(l_ndE), stat=istat)
            if (alloc_failed(istat, 'Allocation failed for derivative callback values', err)) return
         endif

         if (is_quasi) then
            !=========================================================
            ! Case of compex scaled Hamiltonian
            !=========================================================

            !=========================================================
            ! Complex r grid
            !=========================================================
            allocate(l_rcv(1:N), stat=istat)
            if (alloc_failed(istat, 'Allocation failed for l_rcv', err)) return

            if (is_ecs) then 
               l_rcv(1:N1)   = cmplx(yv(1:N1), 0.0_dp, dp) 
               l_rcv(N1+1:N) = yc*(1 - drdyc) + yv(N1+1:N)*drdyc
            else
               l_rcv(1:N)    = rmin*(1 - drdyc) + yv(1:N)*drdyc
            endif 

            !=========================================================
            ! Complex x points for extrapolation to the complex plane
            !=========================================================
            if (is_complex_extrapol) then
               allocate(xcv(1:N+1), stat=istat)
               if (alloc_failed(istat, 'Allocation failed for xcv', err)) return
               if (is_ecs) then
                  xcv(N1+1:N) = mapping%xy(l_rcv(N1+1:N))
               else
                  xcv(1:N) = mapping%xy(l_rcv(1:N))
               endif
               xcv(N+1)    = xmax
            endif

            if (Vf%is()) then
               !=========================================================
               ! Extrapolation real V functions to the complex plane
               !=========================================================
               do k = 1, l_nEq 
                  do m = merge(k, 1, is_V_sym_l), l_nEq
                     call Vf%set_par(k, m) 
                     if (is_ecs) then
                        do j = 1, N1
                           Vcm(j, k, m) = Vf%f(yv(j))
                        enddo
                     endif
                  enddo 
               enddo

               call Vf%set_mapping(mapping)

               do k = 1, l_nEq
                  do m = merge(k, 1, is_V_sym_l), l_nEq
                     call Vf%set_par(k, m)
                     if (is_ecs) then
                        call extrapol_cmplx(xcv(N1+1:N+1), Vcm(N1+1:N+1, k, m), xc, xmax, Vf, err)
                     else
                        call extrapol_cmplx(xcv(1:N+1), Vcm(1:N+1, k, m), xmin, xmax, Vf, err)
                     endif
                     if (err%has_error()) return
                  enddo
               enddo

               call Vf%set_mapping()
            endif

            if (dVf%is()) then
               !=========================================================
               ! Extrapolation real dV functions to the complex plane
               !=========================================================
               if (is_ecs) then
                  do k = 1, l_nEq 
                     do m = merge(k, 1, is_V_sym_l), l_nEq
                        call dVf%set_par(k, m)
                        do j = 1, N1
                           call dVf%fv(yv(j), drv)
                           dVcm(j, :, k, m) = cmplx(drv, 0.0_dp, dp)
                        enddo
                     enddo 
                  enddo
               endif

               call dVf%set_mapping(mapping)

               do k = 1, l_nEq 
                  do m = merge(k, 1, is_V_sym_l), l_nEq
                     call dVf%set_par(k, m)
                     if (is_ecs) then
                        call extrapol_cmplx_vec(xcv(N1+1:N+1), dVcm(N1+1:N+1, :, k, m), xc, xmax, dVf, err)
                     else
                        call extrapol_cmplx_vec(xcv(1:N+1), dVcm(1:N+1, :, k, m), xmin, xmax, dVf, err)
                     endif
                     if (err%has_error()) return
                  enddo
               enddo

               call dVf%set_mapping()
            endif

            if (Bf%is()) then
               !=========================================================
               ! Extrapolation real B-functions to the complex plane
               !=========================================================
               do k = 1, l_nEq 
                  do m = 1, l_nEq
                     call Bf%set_par(k, m)
                     if (is_ecs) then
                        do j = 1, N1
                           Bcm(j, k, m) = Bf%f(yv(j))
                        enddo
                     endif
                  enddo 
               enddo

               call Bf%set_mapping(mapping)

               do k = 1, l_nEq
                  do m = 1, l_nEq
                     call Bf%set_par(k, m)
                     if (is_ecs) then
                        call extrapol_cmplx(xcv(N1+1:N+1), Bcm(N1+1:N+1, k, m), xc, xmax, Bf, err)
                     else
                        call extrapol_cmplx(xcv(1:N+1), Bcm(1:N+1, k, m), xmin, xmax, Bf, err)
                     endif
                     if (err%has_error()) return
                  enddo
               enddo

               call Bf%set_mapping()
            endif

            if (dBf%is()) then
               !=========================================================
               ! Extrapolation real dB-functions to the complex plane
               !=========================================================
               if (is_ecs) then
                  do k = 1, l_nEq
                     do m = 1, l_nEq
                        call dBf%set_par(k, m)
                        do j = 1, N1
                           call dBf%fv(yv(j), drv)
                           dBcm(j, :, k, m) = cmplx(drv, 0.0_dp, dp)
                        enddo
                     enddo
                  enddo
               endif

               call dBf%set_mapping(mapping)

               do k = 1, l_nEq
                  do m = 1, l_nEq
                     call dBf%set_par(k, m)
                     if (is_ecs) then
                        call extrapol_cmplx_vec(xcv(N1+1:N+1), dBcm(N1+1:N+1, :, k, m), xc, xmax, dBf, err)
                     else
                        call extrapol_cmplx_vec(xcv(1:N+1), dBcm(1:N+1, :, k, m), xmin, xmax, dBf, err)
                     endif
                     if (err%has_error()) return
                  enddo
               enddo

               call dBf%set_mapping()
            endif

            if (Vcf%is()) then
               !=========================================================
               ! Calculation of complex V functions at complex r
               !=========================================================
               do k = 1, l_nEq
                  do m = merge(k, 1, is_V_sym_l), l_nEq
                     call Vcf%set_par(k, m)
                     do j = 1, N
                        Vcm(j, k, m) = Vcf%f(l_rcv(j))
                     enddo 
                  enddo 
               enddo
            endif

            if (dVcf%is()) then
               !=========================================================
               ! Calculation of complex dV functions at complex r
               !=========================================================
               do k = 1, l_nEq
                  do m = merge(k, 1, is_V_sym_l), l_nEq
                     call dVcf%set_par(k, m)
                     do j = 1, N
                        call dVcf%fv(l_rcv(j), dcv)
                        dVcm(j, :, k, m) = dcv
                     enddo 
                  enddo 
               enddo
            endif

            if (Bcf%is()) then
               !=========================================================
               ! Calculation of complex B functions at complex r
               !=========================================================
               do k = 1, l_nEq
                  do m = 1, l_nEq
                     call Bcf%set_par(k, m)
                     do j = 1, N
                        Bcm(j, k, m) = Bcf%f(l_rcv(j))
                     enddo 
                  enddo 
               enddo
            endif

            if (dBcf%is()) then
               !=========================================================
               ! Calculation of complex dB functions at complex r
               !=========================================================
               do k = 1, l_nEq
                  do m = 1, l_nEq
                     call dBcf%set_par(k, m)
                     do j = 1, N
                        call dBcf%fv(l_rcv(j), dcv)
                        dBcm(j, :, k, m) = dcv
                     enddo 
                  enddo 
               enddo
            endif

            if (Gf%is()) then
               !=========================================================
               ! Extrapolation real diagonal G functions to the complex plane
               !=========================================================
               do m = 1, l_nEq
                  call Gf%set_par(m)
                  if (is_ecs) then
                     do j = 1, N1
                        Gcm(j, m) = cmplx(Gf%f(yv(j)), 0.0_dp, dp)
                     enddo
                  endif
               enddo

               call Gf%set_mapping(mapping)

               do m = 1, l_nEq
                  call Gf%set_par(m)
                  if (is_ecs) then
                     call extrapol_cmplx(xcv(N1+1:N+1), Gcm(N1+1:N+1, m), xc, xmax, Gf, err)
                  else
                     call extrapol_cmplx(xcv(1:N+1), Gcm(1:N+1, m), xmin, xmax, Gf, err)
                  endif
                  if (err%has_error()) return
               enddo

               call Gf%set_mapping()
            elseif (Gcf%is()) then
               !=========================================================
               ! Calculation of complex diagonal G functions at complex r
               !=========================================================
               do m = 1, l_nEq
                  call Gcf%set_par(m)
                  do j = 1, N
                     Gcm(j, m) = Gcf%f(l_rcv(j))
                  enddo
               enddo
            endif

            if (dGf%is()) then
               !=========================================================
               ! Extrapolation real dG functions to the complex plane
               !=========================================================
               if (is_ecs) then
                  do m = 1, l_nEq
                     call dGf%set_par(m)
                     do j = 1, N1
                        call dGf%fv(yv(j), drv)
                        dGcm(j, :, m) = cmplx(drv, 0.0_dp, dp)
                     enddo
                  enddo
               endif

               call dGf%set_mapping(mapping)

               do m = 1, l_nEq
                  call dGf%set_par(m)
                  if (is_ecs) then
                     call extrapol_cmplx_vec(xcv(N1+1:N+1), dGcm(N1+1:N+1, :, m), xc, xmax, dGf, err)
                  else
                     call extrapol_cmplx_vec(xcv(1:N+1), dGcm(1:N+1, :, m), xmin, xmax, dGf, err)
                  endif
                  if (err%has_error()) return
               enddo

               call dGf%set_mapping()
            elseif (dGcf%is()) then
               !=========================================================
               ! Calculation of complex dG functions at complex r
               !=========================================================
               do m = 1, l_nEq
                  call dGcf%set_par(m)
                  do j = 1, N
                     call dGcf%fv(l_rcv(j), dcv)
                     dGcm(j, :, m) = dcv
                  enddo
               enddo
            endif

            !=========================================================
            ! Applying V symmetry if possible
            !=========================================================
            if (is_V .and. is_V_sym_l) then 
               do k = 1, l_nEq 
                  do m = k+1, l_nEq
                     do j = 1, N
                        Vcm(j, m, k) = Vcm(j, k, m)
                     enddo 
                  enddo 
               enddo
            endif

            !=========================================================
            ! Applying dV symmetry if possible
            !=========================================================
            if (is_dV .and. is_V_sym_l) then 
               do p = 1, l_ndE
                  do k = 1, l_nEq 
                     do m = k+1, l_nEq
                        do j = 1, N
                           dVcm(j, p, m, k) = dVcm(j, p, k, m)
                        enddo 
                     enddo 
                  enddo 
               enddo
            endif

            !=========================================================
            ! Adding centrifugal part to diagonal part of Vm
            !=========================================================
            if (is_Centr) then
               do j = 1, N
                  c0 = JRot*(JRot + 1)/(2*mu*l_rcv(j)**2) 
                  do m = 1, l_nEq
                     Vcm(j, m, m) = Vcm(j, m, m) + c0
                  enddo
               enddo
            endif

         elseif (is_complex) then
            !=========================================================
            ! Complex problem without complex scaling
            !=========================================================

            if (Vf%is()) then
               !=========================================================
               ! Calculation of real V functions at real r-grid 
               !=========================================================
               do k = 1, l_nEq
                  do m = merge(k, 1, is_V_sym_l), l_nEq
                     call Vf%set_par(k, m)
                     do j = 1, N
                        Vcm(j, k, m) = Vf%f(yv(j))
                     enddo 
                  enddo
               enddo

            elseif (Vcf%is()) then 
               !=========================================================
               ! Calculation of complex V function at real r-grid 
               !=========================================================
               do k = 1, l_nEq
                  do m = merge(k, 1, is_V_sym_l), l_nEq
                     call Vcf%set_par(k, m)
                     do j = 1, N
                        Vcm(j, k, m) = Vcf%f(yv(j))
                     enddo 
                  enddo
               enddo
            endif

            if (dVf%is()) then
               !=========================================================
               ! Calculation of real dV functions at real r-grid 
               !=========================================================
               do k = 1, l_nEq
                  do m = merge(k, 1, is_V_sym_l), l_nEq
                     call dVf%set_par(k, m)
                     do j = 1, N
                        call dVf%fv(yv(j), drv)
                        dVcm(j, :, k, m) = cmplx(drv, 0.0_dp, dp)
                     enddo 
                  enddo 
               enddo

            elseif (dVcf%is()) then 
               !=========================================================
               ! Calculation of complex dV function at real r-grid 
               !=========================================================
               do k = 1, l_nEq
                  do m = merge(k, 1, is_V_sym_l), l_nEq
                     call dVcf%set_par(k, m)
                     do j = 1, N
                        call dVcf%fv(cmplx(yv(j), 0.0_dp, dp), dcv)
                        dVcm(j, :, k, m) = dcv
                     enddo 
                  enddo
               enddo
            endif

            if (Bf%is()) then
               !=========================================================
               ! Calculation of real B functions at real r-grid 
               !=========================================================
               do k = 1, l_nEq
                  do m = 1, l_nEq
                     call Bf%set_par(k, m)
                     do j = 1, N
                        Bcm(j, k, m) = Bf%F(yv(j))
                     enddo 
                  enddo
               enddo

            elseif (Bcf%is()) then 
               !=========================================================
               ! Calculation of complex B functions at real r-grid 
               !=========================================================
               do k = 1, l_nEq
                  do m = 1, l_nEq
                     call Bcf%set_par(k, m)
                     do j = 1, N
                        Bcm(j, k, m) = Bcf%f(yv(j))
                     enddo 
                  enddo
               enddo
            endif

            if (dBf%is()) then
               !=========================================================
               ! Calculation of real dB functions at real r-grid 
               !=========================================================
               do k = 1, l_nEq
                  do m = 1, l_nEq
                     call dBf%set_par(k, m)
                     do j = 1, N
                        call dBf%fv(yv(j), drv)
                        dBcm(j, :, k, m) = cmplx(drv, 0.0_dp, dp)
                     enddo 
                  enddo 
               enddo

            elseif (dBcf%is()) then 
               !=========================================================
               ! Calculation of complex dB functions at real r-grid 
               !=========================================================
               do k = 1, l_nEq
                  do m = 1, l_nEq
                     call dBcf%set_par(k, m)
                     do j = 1, N
                        call dBcf%fv(cmplx(yv(j), 0.0_dp, dp), dcv)
                        dBcm(j, :, k, m) = dcv
                     enddo 
                  enddo 
               enddo
            endif

            if (Gf%is()) then
               !=========================================================
               ! Calculation of real diagonal G functions at real r-grid
               !=========================================================
               do m = 1, l_nEq
                  call Gf%set_par(m)
                  do j = 1, N
                     Gcm(j, m) = cmplx(Gf%f(yv(j)), 0.0_dp, dp)
                  enddo
               enddo
            elseif (Gcf%is()) then
               !=========================================================
               ! Calculation of complex diagonal G functions at real r-grid
               !=========================================================
               do m = 1, l_nEq
                  call Gcf%set_par(m)
                  do j = 1, N
                     Gcm(j, m) = Gcf%f(cmplx(yv(j), 0.0_dp, dp))
                  enddo
               enddo
            endif

            if (dGf%is()) then
               !=========================================================
               ! Calculation of real dG functions at real r-grid
               !=========================================================
               do m = 1, l_nEq
                  call dGf%set_par(m)
                  do j = 1, N
                     call dGf%fv(yv(j), drv)
                     dGcm(j, :, m) = cmplx(drv, 0.0_dp, dp)
                  enddo
               enddo
            elseif (dGcf%is()) then
               !=========================================================
               ! Calculation of complex dG functions at real r-grid
               !=========================================================
               do m = 1, l_nEq
                  call dGcf%set_par(m)
                  do j = 1, N
                     call dGcf%fv(cmplx(yv(j), 0.0_dp, dp), dcv)
                     dGcm(j, :, m) = dcv
                  enddo
               enddo
            endif

            !=========================================================
            ! Applying V symmetry if possible
            !=========================================================
            if (is_V .and. is_V_sym_l) then 
               do k = 1, l_nEq 
                  do m = k+1, l_nEq
                     do j = 1, N
                        Vcm(j, m, k) = conjg(Vcm(j, k, m))
                     enddo 
                  enddo 
               enddo
            endif
            !=========================================================
            ! Applying dV symmetry if possible
            !=========================================================
            if (is_dV .and. is_V_sym_l) then 
               do p = 1, l_ndE
                  do k = 1, l_nEq 
                     do m = k+1, l_nEq
                        do j = 1, N
                           dVcm(j, p, m, k) = conjg(dVcm(j, p, k, m))
                        enddo 
                     enddo 
                  enddo 
               enddo
            endif

            !=========================================================
            ! Adding centrifugal part to diagonal part of Vm
            !=========================================================
            if (is_Centr) then
               do j = 1, N
                  t0 = JRot*(JRot + 1)/(2*mu*yv(j)**2) 
                  do m = 1, l_nEq
                     Vcm(j, m, m) = Vcm(j, m, m) + cmplx(t0, 0.0_dp, dp)
                  enddo
               enddo
            endif
         else  
            !=========================================================
            ! Real problem
            !=========================================================
            if (Vf%is()) then
               !=========================================================
               ! Calculation of real V functions at real r-grid 
               !=========================================================
               do k = 1, l_nEq
                  do m = merge(k, 1, is_V_sym_l), l_nEq
                     call Vf%set_par(k, m)
                     do j = 1, N
                        Vm(j, k, m) = Vf%f(yv(j))
                     enddo 
                  enddo
               enddo
            endif

            if (dVf%is()) then
               !=========================================================
               ! Calculation of real dV functions at real r-grid 
               !=========================================================
               do k = 1, l_nEq
                  do m = merge(k, 1, is_V_sym_l), l_nEq
                     call dVf%set_par(k, m)
                     do j = 1, N
                        call dVf%fv(yv(j), drv)
                        dVm(j, :, k, m) = drv
                     enddo 
                  enddo 
               enddo

            elseif (dVcf%is()) then 
               !========================================================
               ! Calculation of complex dV functions at real r-grid 
               ! convert it to real
               !=========================================================
               do k = 1, l_nEq
                  do m = merge(k, 1, is_V_sym_l), l_nEq
                     call dVcf%set_par(k, m)
                     do j = 1, N
                        call dVcf%fv(cmplx(yv(j), 0.0_dp, dp), dcv)
                        dVm(j, :, k, m) = real(dcv, dp)
                     enddo                                                        
                  enddo                                                           
               enddo
            endif

            if (Bf%is()) then
               !=========================================================
               ! Calculation of real B functions at real r-grid 
               !=========================================================
               do k = 1, l_nEq
                  do m = 1, l_nEq
                     call Bf%set_par(k, m)
                     do j = 1, N
                        Bm(j, k, m) = Bf%f(yv(j))
                     enddo 
                  enddo
               enddo
            endif

            if (dBf%is()) then
               !=========================================================
               ! Calculation of real dB functions at real r-grid 
               !=========================================================
               do k = 1, l_nEq
                  do m = 1, l_nEq
                     call dBf%set_par(k, m)
                     do j = 1, N
                        call dBf%fv(yv(j), drv)
                        dBm(j, :, k, m) = drv
                     enddo 
                  enddo
               enddo

            elseif (dBcf%is()) then 
               !=========================================================
               ! Calculation of complex dB functions at real r-grid 
               ! convert it to real
               !=========================================================
               do k = 1, l_nEq
                  do m = 1, l_nEq
                     call dBcf%set_par(k, m)
                     do j = 1, N
                        call dBcf%fv(cmplx(yv(j), 0.0_dp, dp), dcv)
                        dBm(j, :, k, m) = real(dcv, dp)
                     enddo 
                  enddo
               enddo
            endif

            if (Gf%is()) then
               !=========================================================
               ! Calculation of real diagonal G functions at real r-grid
               !=========================================================
               do m = 1, l_nEq
                  call Gf%set_par(m)
                  do j = 1, N
                     Gm(j, m) = Gf%f(yv(j))
                  enddo
               enddo
            endif

            if (dGf%is()) then
               !=========================================================
               ! Calculation of real dG functions at real r-grid
               !=========================================================
               do m = 1, l_nEq
                  call dGf%set_par(m)
                  do j = 1, N
                     call dGf%fv(yv(j), drv)
                     dGm(j, :, m) = drv
                  enddo
               enddo
            elseif (dGcf%is()) then
               !=========================================================
               ! Calculation of complex dG functions at real r-grid
               ! convert it to real
               !=========================================================
               do m = 1, l_nEq
                  call dGcf%set_par(m)
                  do j = 1, N
                     call dGcf%fv(cmplx(yv(j), 0.0_dp, dp), dcv)
                     dGm(j, :, m) = real(dcv, dp)
                  enddo
               enddo
            endif

            !=========================================================
            ! Applying V symmetry if possible
            !=========================================================
            if (is_V .and. is_V_sym_l) then 
               do k = 1, l_nEq 
                  do m = k + 1, l_nEq
                     do j = 1, N
                        Vm(j, m, k) = Vm(j, k, m)
                     enddo 
                  enddo 
               enddo
            endif

            !=========================================================
            ! Applying dV symmetry if possible
            !=========================================================
            if (is_dV .and. is_V_sym_l) then 
               do p = 1, l_ndE
                  do k = 1, l_nEq 
                     do m = k + 1, l_nEq
                        do j = 1, N
                           dVm(j, p, m, k) = dVm(j, p, k, m)
                        enddo 
                     enddo 
                  enddo 
               enddo
            endif

            !=========================================================
            ! Adding centrifugal part to diagonal part of Vm
            !=========================================================
            if (is_Centr) then
               do j = 1, N
                  t0 = JRot*(JRot + 1)/(2*mu*yv(j)**2) 
                  do m = 1, l_nEq
                     Vm(j, m, m) = Vm(j, m, m) + t0
                  enddo
               enddo
            endif
         endif

         if (err%has_error()) return

         !=========================================================
         ! Initializing arrays with functions needed for construction of Hamiltoninan matrix
         !=========================================================
         allocate(gv(N), ginvv(N), stat=istat)
         if (alloc_failed(istat, 'Allocation failed for gv/ginvv', err)) return

         gv(1:N)    = mapping%dydx(xv%x(1:N))
         ginvv(1:N) = mapping%dxdyx(xv%x(1:N))

         if (is_complex) then
            if (is_B) then
               allocate(Fcv(1:N), Qcv(1:N), stat=istat)
            else
               allocate(Fcv(1:N), stat=istat)
            endif
            if (alloc_failed(istat, 'Allocation failed for Fcv/Qcv', err)) return
            Fcv(1:N) = cmplx(mapping%Fx(xv%x(1:N)), 0.0_dp, dp)

            if (is_B) then 
               Qcv(1:N) = cmplx(mapping%Qx(xv%x(1:N)), 0.0_dp, dp)
            endif

            if (is_quasi) then
               if (is_ecs) then
                  Fcv(N1+1:N) = Fcv(N1+1:N)*conjg(drdyc)**2
                  if (is_B) Qcv(N1+1:N) = Qcv(N1+1:N)*conjg(drdyc)
               else
                  Fcv(1:N) = Fcv(1:N)*conjg(drdyc)**2
                  if (is_B) Qcv(1:N) = Qcv(1:N)*conjg(drdyc)
               endif
            endif
         else
            if (is_B) then
               allocate(Fv(1:N), Qv(1:N), stat=istat)
            else
               allocate(Fv(1:N), stat=istat)
            endif
            if (alloc_failed(istat, 'Allocation failed for Fv/Qv', err)) return
            Fv(1:N) = mapping%Fx(xv%x(1:N))
            if (is_B) then 
               Qv(1:N) = mapping%Qx(xv%x(1:N))
            endif
         endif
      end subroutine sample_grid_functions

      !=================================================
      ! Assemble the (real Hm or complex Hcm) Hamiltonian matrix from the
      ! differentiation matrices, potentials and couplings, then apply the
      ! G^{-1/2} similarity transform that reduces the generalized problem to
      ! a standard one. Internal procedure: shares qsolver's arrays and flags
      ! by host association.
      !=================================================
      subroutine build_hamiltonian()
         if (is_complex) then 
            !=========================================================
            ! Complex Hamiltonian Hcm (size nH x nH, where nH = N*l_nEq)
            ! Convention: if is_hermitian=.true., only the UPPER part of Hcm
            !             is required/used (lower part may remain unfilled).
            !=========================================================
            allocate(Hcm(nH,nH), stat=istat)
            if (alloc_failed(istat, 'Allocation failed for Hcm', err)) return

            Hcm = (0.0_dp, 0.0_dp)

            if (is_B) then 
               !=================================================
               ! Build B-dependent part of Hamiltonian matrix
               !=================================================

               if (is_B_sym_l) then
                  !=========================================================
                  ! Build raw upper part (B-symmetric construction)
                  ! Only the upper triangle/blocks are formed here:
                  ! off-diagonal blocks with k>m are fully computed (NxN)
                  ! diagonal blocks are computed only for j>=i
                  !=========================================================
                  do m = 1, l_nEq
                     ! inter-block: k > m
                     do k = m + 1, l_nEq
                        do i = 1, N
                           do j = 1, N
                              Hcm(i+(m-1)*N, j+(k-1)*N) = &
                              Bcm(i, m, k)*ginvv(j) - conjg(Bcm(j, k, m))*ginvv(i)
                           enddo
                        enddo
                     enddo

                     ! intra-block: k = m, j >= i
                     do i = 1, N
                        do j = i, N
                           Hcm(i+(m-1)*N, j+(m-1)*N) = &
                           Bcm(i, m, m)*ginvv(j) - conjg(Bcm(j, m, m))*ginvv(i)
                        enddo
                     enddo

                  enddo

                  if (is_hermitian) then
                     !=================================================
                     ! Hermitian final Hamiltonian: only the upper part is needed
                     ! Apply D1cm * inv2mu only to the upper part that was built above.
                     ! The lower triangle is intentionally left unfilled/unused.
                     !=================================================
                     do m = 1, l_nEq

                        ! scale off-diagonal blocks in the upper block triangle (k>m)
                        do k = m+1, l_nEq
                           Hcm(1+(m-1)*N:m*N, 1+(k-1)*N:k*N) = &
                           &              Hcm(1+(m-1)*N:m*N, 1+(k-1)*N:k*N) * D1cm(1:N,1:N) * inv2mu
                        enddo

                        ! scale diagonal blocks, upper triangle only (j>=i)
                        do i = 1, N
                           do j = i, N
                              Hcm(i+(m-1)*N, j+(m-1)*N) = Hcm(i+(m-1)*N, j+(m-1)*N) * D1cm(i,j) * inv2mu
                           enddo
                        enddo
                     enddo

                  else
                     !=================================================
                     ! Non-Hermitian final Hamiltonian: full matrix is needed
                     ! Complete the lower triangle from the raw skew-Hermitian relation,
                     ! then apply D1cm * inv2mu to all blocks.
                     !=================================================
                     do i = 1, nH-1
                        do j = i+1, nH
                           Hcm(j,i) = -conjg(Hcm(i,j))   ! raw is skew-Hermitian
                        enddo
                     enddo

                     do m = 1, l_nEq
                        do k = 1, l_nEq
                           Hcm(1+(m-1)*N:m*N, 1+(k-1)*N:k*N) = &
                           Hcm(1+(m-1)*N:m*N, 1+(k-1)*N:k*N) * D1cm(1:N,1:N) * inv2mu
                        enddo
                     enddo
                  endif
                  !=================================================
                  ! Add Qcv-dependent diagonal-in-i contribution (B-symmetric form)
                  ! - Hermitian mode: only upper block triangle (k>=m)
                  ! - Otherwise: all blocks (k=1..l_nEq)
                  !=================================================
                  do m = 1, l_nEq
                     do k = merge(m, 1, is_hermitian), l_nEq
                        do i = 1, N
                           Hcm(i+(m-1)*N, i+(k-1)*N) = Hcm(i+(m-1)*N, i+(k-1)*N) + &
                           Qcv(i) * (Bcm(i, m, k) + conjg(Bcm(i, k, m))) * ginvv(i) * inv2mu
                        enddo 
                     enddo
                  enddo

               else
                  !=================================================
                  ! No B symmetry: Hermiticity cannot be assumed
                  ! Build the full matrix directly (all blocks, all columns),
                  ! including D1cm and inv2mu during construction.
                  !=================================================
                  do m = 1, l_nEq
                     do k = 1, l_nEq
                        do j = 1, N
                           Hcm(1+(m-1)*N:m*N, j+(k-1)*N) = &
                           Bcm(1:N, m, k) * D1cm(1:N, j) * (inv2mu * ginvv(j))
                        enddo
                     enddo
                  enddo
                  !=================================================
                  ! Add Qcv-dependent diagonal-in-i contribution (no B-symmetry form)
                  ! Full matrix is updated (all block pairs m,k).
                  !=================================================
                  do m = 1, l_nEq
                     do k = 1, l_nEq
                        do i = 1, N
                           Hcm(i+(m-1)*N, i+(k-1)*N) = Hcm(i+(m-1)*N, i+(k-1)*N) + &
                           Qcv(i) * Bcm(i, m, k) * ginvv(i) * inv2mu
                        enddo
                     enddo
                  enddo

               endif
            endif

            !=================================================
            ! Build and scale diagonal-block contribution D2cm:
            ! 1) add diagonal shift Fcv to D2cm
            ! 2) scale by -ginvv(i)*ginvv(j)*inv2mu (elementwise)
            !=================================================
            do j = 1, N
               D2cm(j,j) = D2cm(j,j) + Fcv(j)
            end do

            do j = 1, N
               do i = 1, N
                  D2cm(i,j) = -D2cm(i,j)*ginvv(i)*ginvv(j)*inv2mu
               enddo 
            enddo

            if (is_hermitian) then 
               !=================================================
               ! Hermitian mode: only upper triangles of diagonal blocks
               ! are required/updated (j>=i).
               !=================================================
               do m = 1, l_nEq
                  do i = 1, N
                     do j = i, N
                        Hcm(i+(m-1)*N, j+(m-1)*N) = Hcm(i+(m-1)*N, j+(m-1)*N) + D2cm(i,j) 
                     enddo
                  enddo
               enddo

            else
               !=================================================
               ! Full-matrix mode: add D2cm to every diagonal block.
               !=================================================
               do k = 1, l_nEq
                  Hcm(1+(k-1)*N:k*N, 1+(k-1)*N:k*N) = Hcm(1+(k-1)*N:k*N, 1+(k-1)*N:k*N) + D2cm(1:N,1:N)
               enddo

            endif
            !=================================================
            ! Add V contribution on "diagonal-in-i" entries:
            ! - Hermitian mode: only upper block triangle (k>=m)
            ! - Otherwise: all blocks (k=1..l_nEq)
            !=================================================
            if (is_V .or. is_Centr) then 
               do m = 1, l_nEq
                  do k = merge(m, 1, is_hermitian), l_nEq
                     do i = 1, N
                        Hcm(i+(m-1)*N, i+(k-1)*N) = Hcm(i+(m-1)*N, i+(k-1)*N) + Vcm(i, m, k)
                     enddo 
                  enddo
               enddo
            endif

         else  
            !=================================================
            ! real case: same structure as complex case, with real arithmetic
            !=================================================
            allocate(Hm(nH,nH), stat=istat)
            if (alloc_failed(istat, 'Allocation failed for Hm', err)) return
            Hm = 0.0_dp

            if (is_B) then 
               !=================================================
               ! Build B-dependent part of Hamiltonian matrix
               !=================================================

               if (is_B_sym_l) then
                  !=================================================
                  ! Build raw upper part (real B-symmetric construction)
                  ! Only the upper triangle/blocks are formed here:
                  !   - off-diagonal blocks with k>m are fully computed (NxN)
                  !   - diagonal blocks are computed only for j>=i
                  !=================================================
                  do m = 1, l_nEq
                     ! inter-block: k > m
                     do k = m + 1, l_nEq
                        do i = 1, N
                           do j = 1, N
                              Hm(i+(m-1)*N, j+(k-1)*N) =  Bm(i, m, k)*ginvv(j) - Bm(j, k, m)*ginvv(i)
                           enddo
                        enddo
                     enddo

                     ! intra-block: k = m, j >= i
                     do i = 1, N
                        do j = i, N
                           Hm(i+(m-1)*N, j+(m-1)*N) =  Bm(i, m, m)*ginvv(j) - Bm(j, m, m)*ginvv(i)
                        enddo
                     enddo

                  enddo

                  if (is_hermitian) then
                     !=================================================
                     ! Hermitian final Hamiltonian (real): only upper part needed
                     ! Apply D1m * inv2mu only to the upper part built above.
                     ! The lower triangle is intentionally left unfilled/unused.
                     !=================================================
                     do m = 1, l_nEq

                        ! scale off-diagonal blocks in the upper block triangle (k>m)
                        do k = m+1, l_nEq
                           Hm(1+(m-1)*N:m*N, 1+(k-1)*N:k*N) = &
                           &              Hm(1+(m-1)*N:m*N, 1+(k-1)*N:k*N) * D1m(1:N,1:N) * inv2mu
                        enddo

                        ! scale diagonal blocks, upper triangle only (j>=i)
                        do i = 1, N
                           do j = i, N
                              Hm(i+(m-1)*N, j+(m-1)*N) = Hm(i+(m-1)*N, j+(m-1)*N) * D1m(i,j) * inv2mu
                           enddo
                        enddo

                     enddo

                  else
                     !=================================================
                     ! Non-Hermitian final Hamiltonian (real): full matrix needed
                     ! Complete the lower triangle from the raw skew-symmetric relation,
                     ! then apply D1m * inv2mu to all blocks.
                     !=================================================
                     do i = 1, nH-1
                        do j = i+1, nH
                           Hm(j,i) = -Hm(i,j)   ! raw is skew-symmetric
                        enddo
                     enddo

                     do m = 1, l_nEq
                        do k = 1, l_nEq
                           Hm(1+(m-1)*N:m*N, 1+(k-1)*N:k*N) = &
                           Hm(1+(m-1)*N:m*N, 1+(k-1)*N:k*N) * D1m(1:N,1:N) * inv2mu
                        enddo
                     enddo
                  endif

                  !=================================================
                  ! Add Qv-dependent diagonal-in-i contribution (real, B-symmetric form)
                  ! - Hermitian mode: only upper block triangle (k>=m)
                  ! - Otherwise: all blocks (k=1..l_nEq)
                  !=================================================
                  do m = 1, l_nEq
                     do k = merge(m, 1, is_hermitian), l_nEq
                        do i = 1, N
                           Hm(i+(m-1)*N, i+(k-1)*N) = Hm(i+(m-1)*N, i+(k-1)*N) + &
                           Qv(i) * (Bm(i, m, k) + Bm(i, k, m)) * ginvv(i) * inv2mu
                        enddo 
                     enddo
                  enddo

               else
                  !=================================================
                  ! No B symmetry (real): Hermiticity cannot be assumed
                  ! Build the full matrix directly (all blocks, all columns),
                  ! including D1m and inv2mu during construction.
                  !=================================================
                  do m = 1, l_nEq
                     do k = 1, l_nEq
                        do j = 1, N
                           Hm(1+(m-1)*N:m*N, j+(k-1)*N) = &
                           Bm(1:N, m, k) * D1m(1:N, j) * (inv2mu * ginvv(j))
                        enddo
                     enddo
                  enddo
                  !=================================================
                  ! Add Qv-dependent diagonal-in-i contribution (real, no B-symmetry form)
                  ! Full matrix is updated (all block pairs m,k).
                  !=================================================
                  do m = 1, l_nEq
                     do k = 1, l_nEq
                        do i = 1, N
                           Hm(i+(m-1)*N, i+(k-1)*N) = Hm(i+(m-1)*N, i+(k-1)*N) + &
                           Qv(i) * Bm(i, m, k) * ginvv(i) * inv2mu
                        enddo
                     enddo
                  enddo

               endif
            endif
            !=================================================
            ! Build and scale diagonal-block contribution D2m:
            ! 1) add diagonal shift Fv to D2m
            ! 2) scale by -ginvv(i)*ginvv(j)*inv2mu (elementwise)
            !=================================================
            do j = 1, N
               D2m(j,j) = D2m(j,j) + Fv(j)
            end do

            do j = 1, N
               do i = 1, N
                  D2m(i,j) = -D2m(i,j)*ginvv(i)*ginvv(j)*inv2mu
               enddo 
            enddo

            if (is_hermitian) then 
               !=================================================
               ! Hermitian mode: only upper triangles of diagonal blocks
               ! are required/updated (j>=i).
               !=================================================
               do m = 1, l_nEq
                  do i = 1, N
                     do j = i, N
                        Hm(i+(m-1)*N, j+(m-1)*N) = Hm(i+(m-1)*N, j+(m-1)*N) + D2m(i,j) 
                     enddo
                  enddo
               enddo

            else
               !=================================================
               ! Full-matrix mode: add D2m to every diagonal block.
               !=================================================
               do k = 1, l_nEq
                  Hm(1+(k-1)*N:k*N, 1+(k-1)*N:k*N) = Hm(1+(k-1)*N:k*N, 1+(k-1)*N:k*N) + D2m(1:N,1:N)
               enddo

            endif

            !=================================================
            ! Add V contribution on "diagonal-in-i" entries:
            ! - Hermitian mode: only upper block triangle (k>=m)
            ! - Otherwise: all blocks (k=1..l_nEq)
            !=================================================
            if (is_V .or. is_Centr) then 
               do m = 1, l_nEq
                  do k = merge(m, 1, is_hermitian), l_nEq
                     do i = 1, N
                        Hm(i+(m-1)*N, i+(k-1)*N) = Hm(i+(m-1)*N, i+(k-1)*N) + Vm(i, m, k)
                     enddo 
                  enddo
               enddo
            endif

         endif

         !=================================================
         ! Generalized eigenproblem H*psi = E*G*psi with diagonal G(r)
         ! is reduced to a standard one A*y = E*y by
         !    A = D*H*D,   D = G^(-1/2),   psi = D*y.
         ! Here G is diagonal in the collocation basis because it is
         ! diagonal in channel space and acts pointwise on the grid.
         !=================================================
         if (is_G) then
            if (is_complex) then
               allocate(Gsqrticm(N, l_nEq), stat=istat)
               if (alloc_failed(istat, 'Allocation failed for Gsqrticm', err)) return

               do m = 1, l_nEq
                  do i = 1, N
                     if (abs(Gcm(i, m)) <= tiny(1.0_dp)) then
                        call err%push_err(ERR_INPUT, 'Diagonal matrix G(r) must not vanish on the grid')
                        return
                     endif
                  enddo
                  Gsqrticm(1:N, m) = 1/sqrt(Gcm(1:N, m))
               enddo

               if (is_hermitian) then
                  do m = 1, l_nEq
                     q = 1 + (m-1)*N

                     do j = 1, N
                        Hcm(q:q+j-1, q+j-1) = Gsqrticm(1:j, m) * Hcm(q:q+j-1, q+j-1) * Gsqrticm(j, m)
                     enddo

                     do k = 1, m-1
                        p = 1 + (k-1)*N
                        do j = 1, N
                           Hcm(p:p+N-1, q+j-1) = Gsqrticm(1:N, k) * Hcm(p:p+N-1, q+j-1) * Gsqrticm(j, m)
                        enddo
                     enddo
                  enddo
               else
                  do m = 1, l_nEq
                     q = 1 + (m-1)*N
                     do k = 1, l_nEq
                        p = 1 + (k-1)*N
                        do j = 1, N
                           Hcm(p:p+N-1, q + j - 1) = Gsqrticm(1:N, k) * Hcm(p:p+N-1, q + j - 1) * Gsqrticm(j, m)
                        enddo
                     enddo
                  enddo
               endif

            else
               allocate(Gsqrtim(N, l_nEq), stat=istat)
               if (alloc_failed(istat, 'Allocation failed for Gsqrtim', err)) return

               do m = 1, l_nEq
                  if (any(Gm(1:N, m) <= 0.0_dp)) then
                     call err%push_err(ERR_INPUT, 'Real diagonal matrix G(r) must be positive on the grid')
                     return
                  endif
                  Gsqrtim(1:N, m) = 1.0_dp/sqrt(Gm(1:N, m))
               enddo

               if (is_hermitian) then
                  do m = 1, l_nEq
                     q = 1 + (m-1)*N

                     do j = 1, N
                        Hm(q:q+j-1, q+j-1) = Gsqrtim(1:j, m) * Hm(q:q+j-1, q+j-1) * Gsqrtim(j, m)
                     enddo

                     do k = 1, m-1
                        p = 1 + (k-1)*N
                        do j = 1, N
                           Hm(p:p+N-1, q+j-1) = Gsqrtim(1:N, k) * Hm(p:p+N-1, q+j-1) * Gsqrtim(j, m)
                        enddo
                     enddo
                  enddo
               else
                  do m = 1, l_nEq
                     q = 1 + (m-1)*N
                     do k = 1, l_nEq
                        p = 1 + (k-1)*N
                        do j = 1, N
                           Hm(p:p+N-1, q+j-1) = Gsqrtim(1:N, k) * Hm(p:p+N-1, q+j-1) * Gsqrtim(j, m)
                        enddo
                     enddo
                  enddo
               endif
            endif
         endif
      end subroutine build_hamiltonian

      !=================================================
      ! Solve the (real/complex, Hermitian/non-Hermitian) matrix eigenproblem
      ! with LAPACK or ARPACK, recover right/left eigenvectors, count the found
      ! pairs (l_nEf) and undo the G^{-1/2} similarity transform. Internal
      ! procedure: shares qsolver's matrices, flags and outputs by host
      ! association.
      !=================================================
      subroutine solve_eigenproblem()
         logical                    :: is_complex_sym   ! is Hamiltonian complex symmetric?
         real(dp),    allocatable   :: Htm(:,:)         ! transpose of Hm (left eigenvectors)
         complex(dp), allocatable   :: Htcm(:,:)        ! conjugate transpose of Hcm

         is_complex_sym = is_complex .and. (.not. is_hermitian) .and. &
         &                 is_B_sym_l .and. is_V_sym_l .and. is_D_sym

         !=================================================
         ! In ARPACK mode, left eigenvectors are found from the
         ! conjugate-transposed Hamiltonian, except for complex-symmetric
         ! problems where they are conjugates of right eigenvectors.
         !=================================================
         if (is_arpack .and. need_lWF) then
            if (is_complex) then       
               if (.not. is_complex_sym) then
                  allocate(Htcm(nH,nH), stat=istat)
                  if (alloc_failed(istat, 'Allocation failed for Htcm', err)) return
                  Htcm = transpose(conjg(Hcm))
               endif
            else
               allocate(Htm(nH,nH), stat=istat)
               if (alloc_failed(istat, 'Allocation failed for Htm', err)) return
               Htm = transpose(Hm)
            endif
         endif

         !=================================================
         ! Finding eigenvalues and eigenvectors of Hamiltonian matrix
         !=================================================
         if (is_complex) then 
            if (is_hermitian) then
               if (is_arpack) then
                  if (need_WF) then
                     call arpack_solve_cmplx_herm(Hcm, nH, nH, E0, Ev, l_nE, Wrcm, err)
                  else
                     call arpack_solve_cmplx_herm(Hcm, nH, nH, E0, Ev, l_nE, err = err)
                  endif
               else
                  if (need_WF) then
                     if (is_irange_lapack) then
                        call lapack_solve_cmplx_herm(Hcm, nH, nH, Ev, imin = l_imin, imax = l_imax, vec = Wrcm, err = err)
                     elseif (is_erange) then
                        call lapack_solve_cmplx_herm(Hcm, nH, nH, Ev, el = l_Emin, eu = l_Emax, vec = Wrcm, err = err)
                     else
                        call lapack_solve_cmplx_herm(Hcm, nH, nH, Ev, vec = Wrcm, err = err)
                     endif
                  else
                     if (is_irange_lapack) then
                        call lapack_solve_cmplx_herm(Hcm, nH, nH, Ev, imin = l_imin, imax = l_imax, err = err)
                     elseif (is_erange) then
                        call lapack_solve_cmplx_herm(Hcm, nH, nH, Ev, el = l_Emin, eu = l_Emax, err = err)
                     else
                        call lapack_solve_cmplx_herm(Hcm, nH, nH, Ev, err = err)
                     endif
                  endif
               endif
               if (err%has_error()) return
            else
               if (is_arpack) then 
                  !=========================================================
                  ! Complex non-Hermitian ARPACK branch.
                  !
                  ! Use the complex-symmetric shift-invert solver only if the
                  ! full discretized matrix is complex symmetric. Symmetric V/B
                  ! alone is not enough: the differentiation matrices must also
                  ! preserve symmetry (is_D_sym = .true.).
                  !=========================================================
                  if (is_complex_sym) then
                     if (need_WF) then
                        call arpack_solve_cmplx_symm_nonherm(Hcm, nH, nH, Ec0, Ecv, l_nE, Wrcm, err)
                     else
                        call arpack_solve_cmplx_symm_nonherm(Hcm, nH, nH, Ec0, Ecv, l_nE, err = err)
                     endif
                     if (err%has_error()) return
                     if (need_lWF) then 
                        allocate(Wlcm(size(Wrcm,1), size(Wrcm,2)), stat=istat)
                        if (alloc_failed(istat, 'Allocation failed for Wlcm in complex-symmetric ARPACK branch', err)) return
                        Wlcm = conjg(Wrcm)
                        allocate(Etcv(size(Ecv)), stat=istat)
                        if (alloc_failed(istat, 'Allocation failed for Etcv in complex-symmetric ARPACK branch', err)) return
                        Etcv = Ecv
                     endif
                  else
                     if (need_WF) then
                        call arpack_solve_cmplx_nonherm(Hcm, nH, nH, Ec0, Ecv, l_nE, Wrcm, err)
                     else
                        call arpack_solve_cmplx_nonherm(Hcm, nH, nH, Ec0, Ecv, l_nE, err = err)
                     endif
                     if (err%has_error()) return
                     if (need_lWF) then 
                        call arpack_solve_cmplx_nonherm(Htcm, nH, nH, conjg(Ec0), Etcv, l_nE, Wlcm, err)
                        if (err%has_error()) return
                        Etcv = conjg(Etcv)
                     endif
                  endif
               else
                  if (need_WF .and. need_lWF) then
                     if (is_complex_sym) then
                        ! Complex-symmetric case: left eigenvectors are conjugates
                        ! of right eigenvectors.
                        call lapack_solve_cmplx_nonherm(Hcm, nH, nH, Ecv, vecr = Wrcm, err = err)
                        if (err%has_error()) return
                        allocate(Wlcm(size(Wrcm,1), size(Wrcm,2)), stat=istat)
                        if (alloc_failed(istat, 'Allocation failed for Wlcm in complex-symmetric LAPACK branch', err)) return
                        Wlcm = conjg(Wrcm)
                     else
                        call lapack_solve_cmplx_nonherm(Hcm, nH, nH, Ecv, vecl = Wlcm, vecr = Wrcm, err = err)
                     endif
                  elseif (need_WF) then
                     call lapack_solve_cmplx_nonherm(Hcm, nH, nH, Ecv, vecr = Wrcm, err = err)
                  else
                     call lapack_solve_cmplx_nonherm(Hcm, nH, nH, Ecv, err = err)
                  endif
                  if (err%has_error()) return
               endif 
            endif
         else
            if (is_hermitian) then
               if (is_arpack) then
                  if (need_WF) then
                     call arpack_solve_real_sym(Hm, nH, nH, E0, Ev, l_nE, Wm, err)
                  else
                     call arpack_solve_real_sym(Hm, nH, nH, E0, Ev, l_nE, err = err)
                  endif
               else
                  if (need_WF) then
                     if (is_irange_lapack) then
                        call lapack_solve_real_sym(Hm, nH, nH, Ev, imin = l_imin, imax = l_imax, vec = Wm, err = err)
                     elseif (is_erange) then
                        call lapack_solve_real_sym(Hm, nH, nH, Ev, el = l_Emin, eu = l_Emax, vec = Wm, err = err)
                     else
                        call lapack_solve_real_sym(Hm, nH, nH, Ev, vec = Wm, err = err)
                     endif
                  else
                     if (is_irange_lapack) then
                        call lapack_solve_real_sym(Hm, nH, nH, Ev, imin = l_imin, imax = l_imax, err = err)
                     elseif (is_erange) then
                        call lapack_solve_real_sym(Hm, nH, nH, Ev, el = l_Emin, eu = l_Emax, err = err)
                     else
                        call lapack_solve_real_sym(Hm, nH, nH, Ev, err = err)
                     endif
                  endif
               endif
               if (err%has_error()) return
            else
               if (is_arpack) then 
                  if (need_WF) then
                     call arpack_solve_real_nonsym(Hm, nH, nH, E0, Ecv, l_nE, Wrcm, err)
                  else
                     call arpack_solve_real_nonsym(Hm, nH, nH, E0, Ecv, l_nE, err = err)
                  endif
                  if (err%has_error()) return
                  if (need_lWF) then 
                     call arpack_solve_real_nonsym(Htm, nH, nH, E0, Etcv, l_nE, Wlcm, err)
                     if (err%has_error()) return
                     Etcv = conjg(Etcv)
                  endif
               else
                  if (need_WF .and. need_lWF) then
                     call lapack_solve_real_nonsym(Hm, nH, nH, Ecv, vecl = Wlcm, vecr = Wrcm, err = err)
                  elseif (need_WF) then
                     call lapack_solve_real_nonsym(Hm, nH, nH, Ecv, vecr = Wrcm, err = err)
                  else
                     call lapack_solve_real_nonsym(Hm, nH, nH, Ecv, err = err)
                  endif
                  if (err%has_error()) return
               endif 
            endif
         endif

         !=================================================
         ! How many eigenvalues found?
         !=================================================
         if (is_hermitian) then
            l_nEf = size(Ev)
         else 
            l_nEf = size(Ecv)
         endif

         !=================================================
         ! Recover eigenvectors of the original generalized problem:
         !   psi = D*y,   lpsi = D^dag*l_y
         !=================================================
         if (is_G .and. need_WF) then
            if (is_complex) then
               if (allocated(Wrcm)) then
                  do j = 1, size(Wrcm, 2)
                     do m = 1, l_nEq
                        q = 1 + (m-1)*N
                        Wrcm(q:q+N-1, j) = Gsqrticm(1:N, m) * Wrcm(q:q+N-1, j)
                     enddo
                  enddo
               endif
               if (allocated(Wlcm)) then
                  do j = 1, size(Wlcm, 2)
                     do m = 1, l_nEq
                        q = 1 + (m-1)*N
                        Wlcm(q:q+N-1, j) = conjg(Gsqrticm(1:N, m)) * Wlcm(q:q+N-1, j)
                     enddo
                  enddo
               endif
            else
               if (allocated(Wm)) then
                  do j = 1, size(Wm, 2)
                     do m = 1, l_nEq
                        q = 1 + (m-1)*N
                        Wm(q:q+N-1, j) = Gsqrtim(1:N, m) * Wm(q:q+N-1, j)
                     enddo
                  enddo
               endif
               if (allocated(Wrcm)) then
                  do j = 1, size(Wrcm, 2)
                     do m = 1, l_nEq
                        q = 1 + (m-1)*N
                        Wrcm(q:q+N-1, j) = cmplx(Gsqrtim(1:N, m), 0.0_dp, dp) * Wrcm(q:q+N-1, j)
                     enddo
                  enddo
               endif
               if (allocated(Wlcm)) then
                  do j = 1, size(Wlcm, 2)
                     do m = 1, l_nEq
                        q = 1 + (m-1)*N
                        Wlcm(q:q+N-1, j) = cmplx(Gsqrtim(1:N, m), 0.0_dp, dp) * Wlcm(q:q+N-1, j)
                     enddo
                  enddo
               endif
            endif
         endif
      end subroutine solve_eigenproblem

      !=================================================
      ! Flag spurious ("ghost") eigenvalues via the Legendre spectral-decay
      ! test, setting is_Ev(:). Internal procedure: it shares qsolver's
      ! grids, eigenvectors and flags by host association.
      !=================================================
      subroutine filter_spurious_states()
         integer                  :: nfl, nfh   ! non-ECS projection node range
         real(dp),    allocatable :: ar(:)      ! |Legendre coefficients| (decay test)
         complex(dp), allocatable :: ac(:)      ! complex Legendre coefficients
         real(dp),    allocatable :: Wv(:)      ! combined real wavefunction on the mesh
         complex(dp), allocatable :: Wcv(:)     ! combined complex wavefunction on the mesh

         !=================================================
         ! Allocate the keep/discard mask; initially assume every state is a
         ! true energy. When the spectral filter is switched off, that is the
         ! final mask, so return right away.
         !=================================================
         allocate(is_Ev(l_nEf), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in eigenvalue filter block', err)) return
         is_Ev = .true.

         if (.not. is_efilter) return

         !=================================================
         ! Applying spectral decay filter to discard spurious (ghost) eigenvalues
         !
         ! In case of exterior complex scaling, the spectral filter is applied
         ! interval by interval. The outer interval [xc; xmax] is checked first
         ! if its amplitude is non-negligible; the inner interval is then
         ! checked only if it is also non-negligible.
         !=================================================
         !=================================================
         ! The Legendre coefficients are produced node-by-node by the
         ! grid-aware projection routines (legendre_coeffs_*), which
         ! take the grid itself (nodes + weights + interval) and the
         ! sampled values, so no separate reference grid/table is kept here.
         ! For a single interval, only the projection window has to be
         ! chosen: a Lobatto grid (both ends Dirichlet) spans its full node
         ! set g0%lo:g0%hi, including the zero boundary nodes; a
         ! Radau/Legendre grid uses the N interior solution nodes 1:N.
         !=================================================
         if (.not. is_ecs) then
            if (g0%kind == GRID_LOBATTO) then
               nfl = g0%lo;  nfh = g0%hi
            else
               nfl = 1;      nfh = N
            endif
         endif

         if (is_complex_WF) then
            allocate(ar(0:N+1), ac(0:N+1), Wcv(0:N+1), stat=istat)
         else
            allocate(ar(0:N+1), Wv(0:N+1), stat=istat)
         endif
         if (alloc_failed(istat, 'Allocation failed in eigenvalue filter block', err)) return

         do q = 1, l_nEf
            if (is_complex_WF) then
               Wcv = (0.0_dp, 0.0_dp)

               do m = 1, l_nEq
                  Wcv(1:N) = Wcv(1:N) + Wrcm(1+(m-1)*N:m*N, q)
               enddo

               if (is_D_sym) then
                  Wcv(1:N) = Wcv(1:N)/bary(1:N)
               endif

               !=================================================
               ! For exterior complex scaling, compare amplitudes on
               ! the inner and outer intervals and skip the spectral check on
               ! an interval whose amplitude is negligible relative to the
               ! larger one. A negligible outer tail is ignored; if the inner
               ! part is negligible after the outer check, keep the state.
               !=================================================
               if (is_ecs) then
                  t1 = maxval(abs(Wcv(1:N1)))
                  t2 = maxval(abs(Wcv(N1:N+1)))
                  t0 = max(t1, t2)
                  if (t0 <= 0.0_dp) then
                     is_Ev(q) = .false.
                     cycle
                  endif

                  if (t2/t0 > l_f_thr) then
                     call legendre_coeffs_cmplx(grid_outer, 0, N2, Wcv(N1:N+1), ac(0:N2), err)
                     if (err%has_error()) return
                     ar(0:N2) = abs(ac(0:N2))
                     if (decay_tail_exceeds(ar(0:N2), l_f_thr)) then
                        is_Ev(q) = .false.
                        cycle
                     endif
                  endif

                  if (t1/t0 < l_f_thr) cycle
               endif

               if (is_ecs) then
                  if (l_bc_left == BC_REGULAR) then
                     !=================================================
                     ! Active inner nodes are grid_inner%x(1:N1) (Radau, endpoint = +1);
                     ! the inner polynomial has degree N1-1, not N1.
                     !=================================================
                     call legendre_coeffs_cmplx(grid_inner, 1, N1, Wcv(1:N1), ac(0:N1-1), err)
                  else
                     call legendre_coeffs_cmplx(grid_inner, 0, N1, Wcv(0:N1), ac(0:N1), err)
                  endif
                  if (err%has_error()) return
                  ar(0:N1) = abs(ac(0:N1))
               else
                  call legendre_coeffs_cmplx(g0, nfl, nfh, Wcv(nfl:nfh), ac(0:nfh-nfl), err)
                  if (err%has_error()) return
                  ar(0:nfh-nfl) = abs(ac(0:nfh-nfl))
               endif
            else
               Wv = 0.0_dp

               do m = 1, l_nEq
                  Wv(1:N) = Wv(1:N) + Wm(1+(m-1)*N:m*N, q)
               enddo

               if (is_D_sym) then
                  Wv(1:N) = Wv(1:N) / bary(1:N)
               endif

               call legendre_coeffs_real(g0, nfl, nfh, Wv(nfl:nfh), ar(0:nfh-nfl), err)
               if (err%has_error()) return
               ar = abs(ar)
            endif

            !=================================================
            ! Highest Legendre index actually populated by the projection.
            ! ECS: for the Radau inner subinterval (BC_REGULAR at left) the
            ! polynomial has degree N1-1, otherwise N1. Non-ECS: the degree
            ! is set by the grid-derived projection window, nfh - nfl.
            !=================================================
            if (is_ecs) then
               if (l_bc_left == BC_REGULAR) then
                  k = N1 - 1
               else
                  k = N1
               endif
            else
               k = nfh - nfl
            endif

            if (decay_tail_exceeds(ar(0:k), l_f_thr)) then
               is_Ev(q) = .false.
               cycle
            endif
         enddo
      end subroutine filter_spurious_states

      !=================================================
      ! Internal procedure (see the call site for what this phase does);
      ! shares qsolver state by host association.
      !=================================================
      subroutine sort_eigenpairs()
         if (is_hermitian) then 
            if (use_target_sort) then
               call sort_real_energies(Ev, Eindv, E0 = E0, sort_mode = 1, err = err)
            else
               call sort_real_energies(Ev, Eindv, sort_mode = l_sortE_mode, err = err)
            endif
            if (err%has_error()) return
         else
            if (use_target_sort) then
               call sort_complex_energies(Ecv, Eindv, Ec0 = Ec0, err = err)
            else
               call sort_complex_energies(Ecv, Eindv, sortE = l_sortE_mode, sortW = l_sortW_mode, err = err)
            endif
            if (err%has_error()) return
         endif

         if (is_arpack .and. need_lWF) then
            if (size(Etcv) .ne. size(Ecv)) then
               call err%push_err(ERR_INTERNAL, 'Something wrong in arpack calculation of left eigenvectors')
               return
            endif
            if (use_target_sort) then
               call sort_complex_energies(Etcv, Etindv, Ec0 = Ec0, err = err)
            else
               call sort_complex_energies(Etcv, Etindv, sortE = l_sortE_mode, sortW = l_sortW_mode, err = err)
            endif
            if (err%has_error()) return
         endif

         !=================================================
         ! Removing spurious ones and which are out of boundaries from index array
         !=================================================
         m = 0
         do i = 1, l_nEf
            if (is_hermitian) then
               t0 = Ev(Eindv(i))
               if (is_Ev(Eindv(i)) .and. l_Emin <= t0 .and. t0 <= l_Emax) then
                  m = m + 1
                  Eindv(m) = Eindv(i)
                  if (is_arpack .and. need_lWF) Etindv(m) = Etindv(i)
               end if
            else
               t0 = real(Ecv(Eindv(i)), dp)
               t1 = -2*aimag(Ecv(Eindv(i)))
               if (is_Ev(Eindv(i)) .and. l_Emin <= t0 .and. t0 <= l_Emax .and. t1 <= l_Wmax) then
                  m = m + 1
                  Eindv(m) = Eindv(i)
                  if (is_arpack .and. need_lWF) Etindv(m) = Etindv(i)
               end if
            endif
         end do

         if (is_irange .and. .not. is_irange_lapack) then
            if (m < l_imin) then
               m = 0
            else
               k = min(m, l_imax) - l_imin + 1
               do i = 1, k
                  Eindv(i) = Eindv(l_imin + i - 1)
                  if (is_arpack .and. need_lWF) then
                     Etindv(i) = Etindv(l_imin + i - 1)
                  endif
               enddo
               m = k
            endif
         endif

         l_nEf = m
         if (present(nE)) then
            if (nE > 0) l_nEf = min(l_nEf, nE)
         elseif (present(E0)) then
            l_nEf = min(l_nEf, 1)
         endif
      end subroutine sort_eigenpairs

      !=================================================
      ! Internal procedure (see the call site for what this phase does);
      ! shares qsolver state by host association.
      !=================================================
      subroutine compute_energy_derivatives()
         complex(dp), allocatable :: dHcm(:,:)   ! complex derivative Hamiltonian
         real(dp),    allocatable :: dHm(:,:)     ! real derivative Hamiltonian
         real(dp),    allocatable :: tv(:)        ! real    work vector (dH * W)
         complex(dp), allocatable :: tcv(:)       ! complex work vector (dH * W)

         if (is_dE) then

            allocate(tcv(1:nH), tv(1:nH), stat=istat)
            if (alloc_failed(istat, 'Allocation failed for tcv/tv', err)) return

            if (present(dE)) then 
               if (size(dE, 1) < l_ndE .or. size(dE, 2) < l_nEf) then
                  call err%push_err(ERR_INPUT, 'dE array is too small')
                  return
               endif
               dE(1:l_ndE, 1:l_nEf) = 0.0_dp 
            endif

            if (present(dW)) then 
               if (size(dW, 1) < l_ndE .or. size(dW, 2) < l_nEf) then
                  call err%push_err(ERR_INPUT, 'dW array is too small')
                  return
               endif
               dW(1:l_ndE, 1:l_nEf) = 0.0_dp 
            endif

            if (is_dV) then 
               !=========================================================
               ! Calculation derivatives with respect to V functions
               ! dE(p, q) = <Wl_q| dH/da_p |Wr_q> / <Wl_q | Wr_q>
               ! Hermitian case: Wl = Wr
               !=========================================================
               if (is_complex) then 
                  do p = 1, l_ndE
                     do q = 1, l_nEf
                        !=========================================================
                        ! Calculation of dH*Wr vector
                        !=========================================================
                        do m = 1, l_nEq 
                           do i = 1, N 
                              tcv(i + (m-1)*N) = &
                              sum(dVcm(i, p, m, 1:l_nEq) * Wrcm(i:i+N*(l_nEq-1):N, Eindv(q)))
                           enddo
                        enddo

                        m = Eindv(q)
                        if (is_hermitian) then
                           if (is_G) then
                              c0 = (0.0_dp, 0.0_dp)
                              do k = 1, l_nEq
                                 c0 = c0 + dot_product(Wrcm(1+(k-1)*N:k*N, m), &
                                    &                   Gcm(1:N, k)*Wrcm(1+(k-1)*N:k*N, m))
                              enddo
                              t0 = real(c0, dp)
                           else
                              t0 = real(dot_product(Wrcm(:, m), Wrcm(:, m)), dp)
                           endif
                           t1 = real(dot_product(Wrcm(:, m), tcv(:)), dp)
                        elseif (is_arpack) then
                           if (is_G) then
                              c0 = (0.0_dp, 0.0_dp)
                              do k = 1, l_nEq
                                 c0 = c0 + dot_product(Wlcm(1+(k-1)*N:k*N, Etindv(q)), &
                                    &                   Gcm(1:N, k)*Wrcm(1+(k-1)*N:k*N, m))
                              enddo
                           else
                              c0 = dot_product(Wlcm(:, Etindv(q)), Wrcm(:, m))
                           endif
                           c1 = dot_product(Wlcm(:, Etindv(q)), tcv(:))
                        else
                           if (is_G) then
                              c0 = (0.0_dp, 0.0_dp)
                              do k = 1, l_nEq
                                 c0 = c0 + dot_product(Wlcm(1+(k-1)*N:k*N, m), &
                                    &                   Gcm(1:N, k)*Wrcm(1+(k-1)*N:k*N, m))
                              enddo
                           else
                              c0 = dot_product(Wlcm(:, m), Wrcm(:, m))
                           endif
                           c1 = dot_product(Wlcm(:, m), tcv(:))
                        endif
                        if (is_hermitian) then
                           if (present(dE)) dE(p, q) = t1/t0
                           if (present(dW)) dW(p, q) = 0.0_dp
                        else
                           if (present(dE)) dE(p, q) = real(c1/c0, dp)
                           if (present(dW)) dW(p, q) = -2*aimag(c1/c0)
                        endif
                     enddo
                  enddo
               else
                  do p = 1, l_ndE
                     do q = 1, l_nEf
                        !=========================================================
                        ! Calculation of dH*Wr vector
                        !=========================================================
                        if (is_hermitian) then
                           do m = 1, l_nEq 
                              do i = 1, N 
                                 tv(i + (m-1)*N) = &
                                 sum(dVm(i, p, m, 1:l_nEq)*Wm(i:i+N*(l_nEq-1):N, Eindv(q)))
                              enddo
                           enddo
                        else 
                           do m = 1, l_nEq 
                              do i = 1, N 
                                 tcv(i + (m-1)*N) = &
                                 sum(dVm(i, p, m, 1:l_nEq)*Wrcm(i:i+N*(l_nEq-1):N, Eindv(q)))
                              enddo
                           enddo
                        endif

                        m = Eindv(q)

                        if (is_hermitian) then
                           if (is_G) then
                              t0 = 0.0_dp
                              do k = 1, l_nEq
                                 t0 = t0 + dot_product(Wm(1+(k-1)*N:k*N, m), &
                                    &                  Gm(1:N, k)*Wm(1+(k-1)*N:k*N, m))
                              enddo
                           else
                              t0 = dot_product(Wm(:, m), Wm(:, m))
                           endif
                           t1 = dot_product(Wm(:, m), tv)
                        elseif (is_arpack) then
                           if (is_G) then
                              c0 = (0.0_dp, 0.0_dp)
                              do k = 1, l_nEq
                                 c0 = c0 + dot_product(Wlcm(1+(k-1)*N:k*N, Etindv(q)), &
                                    &                   Gm(1:N, k)*Wrcm(1+(k-1)*N:k*N, m))
                              enddo
                           else
                              c0 = dot_product(Wlcm(:, Etindv(q)), Wrcm(:, m))
                           endif
                           c1 = dot_product(Wlcm(:, Etindv(q)), tcv)
                        else
                           if (is_G) then
                              c0 = (0.0_dp, 0.0_dp)
                              do k = 1, l_nEq
                                 c0 = c0 + dot_product(Wlcm(1+(k-1)*N:k*N, m), &
                                    &                   Gm(1:N, k)*Wrcm(1+(k-1)*N:k*N, m))
                              enddo
                           else
                              c0 = dot_product(Wlcm(:, m), Wrcm(:, m))
                           endif
                           c1 = dot_product(Wlcm(:, m), tcv)
                        endif

                        if (is_hermitian) then
                           if (present(dE)) dE(p, q) = t1/t0
                           if (present(dW)) dW(p, q) = 0.0_dp
                        else
                           if (present(dE)) dE(p, q) = real(c1/c0, dp)
                           if (present(dW)) dW(p, q) = -2*aimag(c1/c0)
                        endif

                     enddo
                  enddo
               endif
            endif

            if (is_dB) then 
               !=========================================================
               ! Calculation derivatives with respect to B functions
               ! dE(p, q) = <Wl_q| dH/da_p |Wr_q> / <Wl_q | Wr_q>
               ! Hermitian: Wl = Wr
               !=========================================================
               if (is_complex) then
                  allocate(dHcm(nH, nH), stat=istat)
                  if (alloc_failed(istat, 'Allocation failed for dHcm', err)) return

                  do p = 1, l_ndE

                     dHcm = (0.0_dp, 0.0_dp)

                     if (is_B_sym_l) then
                        !========================
                        ! Build raw upper part (B-symmetric construction) for dH
                        !========================
                        do m = 1, l_nEq

                           ! inter-block: k > m
                           do k = m + 1, l_nEq
                              do i = 1, N
                                 do j = 1, N
                                    dHcm(i+(m-1)*N, j+(k-1)*N) = &
                                    dBcm(i, p, m, k)*ginvv(j) - conjg(dBcm(j, p, k, m))*ginvv(i)
                                 enddo
                              enddo
                           enddo

                           ! intra-block: k = m, j >= i
                           do i = 1, N
                              do j = i, N
                                 dHcm(i+(m-1)*N, j+(m-1)*N) = &
                                 dBcm(i, p, m, m)*ginvv(j) - conjg(dBcm(j, p, m, m))*ginvv(i)
                              enddo
                           enddo

                        enddo

                        if (is_hermitian) then
                           !=================================================
                           ! Hermitian: scale only the upper part that was built
                           !=================================================
                           do m = 1, l_nEq

                              do k = m + 1, l_nEq
                                 dHcm(1+(m-1)*N:m*N, 1+(k-1)*N:k*N) = &
                                 dHcm(1+(m-1)*N:m*N, 1+(k-1)*N:k*N) * D1cm(1:N, 1:N) * inv2mu
                              enddo

                              do i = 1, N
                                 do j = i, N
                                    dHcm(i+(m-1)*N, j+(m-1)*N) = &
                                    dHcm(i+(m-1)*N, j+(m-1)*N) * D1cm(i, j) * inv2mu
                                 enddo
                              enddo

                           enddo

                        else
                           !=================================================
                           ! Non-Hermitian: complete lower from raw skew-Hermitian,
                           ! then scale all blocks
                           !=================================================
                           do i = 1, nH - 1
                              do j = i + 1, nH
                                 dHcm(j, i) = -conjg(dHcm(i, j))
                              enddo
                           enddo

                           do m = 1, l_nEq
                              do k = 1, l_nEq
                                 dHcm(1+(m-1)*N:m*N, 1+(k-1)*N:k*N) = &
                                 dHcm(1+(m-1)*N:m*N, 1+(k-1)*N:k*N) * D1cm(1:N, 1:N) * inv2mu
                              enddo
                           enddo

                        endif
                        !=================================================
                        ! Add Qv-dependent diagonal-in-i contribution
                        !=================================================
                        do m = 1, l_nEq
                           do k = merge(m, 1, is_hermitian), l_nEq
                              do i = 1, N
                                 dHcm(i+(m-1)*N, i+(k-1)*N) = &
                                 dHcm(i+(m-1)*N, i+(k-1)*N) + &
                                 Qcv(i) * (dBcm(i, p, m, k) + conjg(dBcm(i, p, k, m))) * &
                                 ginvv(i) * inv2mu
                              enddo
                           enddo
                        enddo
                        if (is_hermitian) then
                           ! Complete Hermitian lower triangle for matmul
                           do i = 2, nH
                              do j = 1, i - 1
                                 dHcm(i, j) = conjg(dHcm(j, i))
                              enddo
                           enddo
                        endif

                     else
                        !=================================================
                        ! No B symmetry: build full dH directly
                        !=================================================
                        do m = 1, l_nEq
                           do k = 1, l_nEq
                              do j = 1, N
                                 dHcm(1+(m-1)*N:m*N, j+(k-1)*N) = &
                                 dBcm(1:N, p, m, k) * D1cm(1:N, j) * (inv2mu * ginvv(j))
                              enddo
                           enddo
                        enddo
                        ! Qv-dependent diagonal-in-i contribution (no B-symmetry form)
                        do m = 1, l_nEq
                           do k = 1, l_nEq
                              do i = 1, N
                                 dHcm(i+(m-1)*N, i+(k-1)*N) = &
                                 dHcm(i+(m-1)*N, i+(k-1)*N) + &
                                 Qcv(i) * dBcm(i, p, m, k) * ginvv(i) * inv2mu
                              enddo
                           enddo
                        enddo

                     endif
                     !=========================================================
                     ! Now compute dE(p, q) with dot_product 
                     !=========================================================
                     do q = 1, l_nEf

                        m = Eindv(q)
                        tcv = matmul(dHcm, Wrcm(:, m))

                        if (is_hermitian) then
                           if (is_G) then
                              c0 = (0.0_dp, 0.0_dp)
                              do k = 1, l_nEq
                                 c0 = c0 + dot_product(Wrcm(1+(k-1)*N:k*N, m), &
                                    &                   Gcm(1:N, k)*Wrcm(1+(k-1)*N:k*N, m))
                              enddo
                              t0 = real(c0, dp)
                           else
                              t0 = real(dot_product(Wrcm(:, m), Wrcm(:, m)), dp)
                           endif
                           t1 = real(dot_product(Wrcm(:, m), tcv(:)), dp)
                        elseif (is_arpack) then
                           if (is_G) then
                              c0 = (0.0_dp, 0.0_dp)
                              do k = 1, l_nEq
                                 c0 = c0 + dot_product(Wlcm(1+(k-1)*N:k*N, Etindv(q)), &
                                    &                   Gcm(1:N, k)*Wrcm(1+(k-1)*N:k*N, m))
                              enddo
                           else
                              c0 = dot_product(Wlcm(:, Etindv(q)), Wrcm(:, m))
                           endif
                           c1 = dot_product(Wlcm(:, Etindv(q)), tcv(:))
                        else
                           if (is_G) then
                              c0 = (0.0_dp, 0.0_dp)
                              do k = 1, l_nEq
                                 c0 = c0 + dot_product(Wlcm(1+(k-1)*N:k*N, m), &
                                    &                   Gcm(1:N, k)*Wrcm(1+(k-1)*N:k*N, m))
                              enddo
                           else
                              c0 = dot_product(Wlcm(:, m), Wrcm(:, m))
                           endif
                           c1 = dot_product(Wlcm(:, m), tcv(:))
                        endif
                        if (is_hermitian) then
                           if (present(dE)) dE(p, q) = dE(p, q) + t1/t0
                           if (present(dW)) dW(p, q) = dW(p, q) + 0.0_dp
                        else
                           if (present(dE)) dE(p, q) = dE(p, q) + real(c1/c0, dp)
                           if (present(dW)) dW(p, q) = dW(p, q) - 2*aimag(c1/c0)
                        endif
                     enddo
                  enddo
               else
                  !=========================================================
                  ! Case of real Hamiltonian
                  ! Wavefunctions steel can by complex
                  !=========================================================
                  allocate(dHm(nH, nH), stat=istat)
                  if (alloc_failed(istat, 'Allocation failed for dHm', err)) return

                  do p = 1, l_ndE

                     dHm = 0.0_dp

                     if (is_B_sym_l) then
                        !=========================================================
                        ! Build raw upper part (B-symmetric construction) for dH
                        !=========================================================
                        do m = 1, l_nEq

                           ! inter-block: k > m
                           do k = m + 1, l_nEq
                              do i = 1, N
                                 do j = 1, N
                                    dHm(i+(m-1)*N, j+(k-1)*N) = &
                                    dBm(i, p, m, k)*ginvv(j) - dBm(j, p, k, m)*ginvv(i)
                                 enddo
                              enddo
                           enddo

                           ! intra-block: k = m, j >= i
                           do i = 1, N
                              do j = i, N
                                 dHm(i+(m-1)*N, j+(m-1)*N) = &
                                 dBm(i, p, m, m)*ginvv(j) - dBm(j, p, m, m)*ginvv(i)
                              enddo
                           enddo

                        enddo

                        if (is_hermitian) then
                           !=================================================
                           ! Hermitian (real): scale only the upper part that was built
                           !=================================================
                           do m = 1, l_nEq

                              do k = m + 1, l_nEq
                                 dHm(1+(m-1)*N:m*N, 1+(k-1)*N:k*N) = &
                                 dHm(1+(m-1)*N:m*N, 1+(k-1)*N:k*N) * D1m(1:N, 1:N) * inv2mu
                              enddo

                              do i = 1, N
                                 do j = i, N
                                    dHm(i+(m-1)*N, j+(m-1)*N) = &
                                    dHm(i+(m-1)*N, j+(m-1)*N) * D1m(i, j) * inv2mu
                                 enddo
                              enddo

                           enddo

                        else
                           !=================================================
                           ! Non-Hermitian (real): complete lower from raw skew-symmetric,
                           ! then scale all blocks
                           !=================================================
                           do i = 1, nH - 1
                              do j = i + 1, nH
                                 dHm(j, i) = -dHm(i, j)
                              enddo
                           enddo

                           do m = 1, l_nEq
                              do k = 1, l_nEq
                                 dHm(1+(m-1)*N:m*N, 1+(k-1)*N:k*N) = &
                                 dHm(1+(m-1)*N:m*N, 1+(k-1)*N:k*N) * D1m(1:N, 1:N) * inv2mu
                              enddo
                           enddo

                        endif
                        !=================================================
                        ! Add Qv-dependent diagonal-in-i contribution (real, B-symmetric form)
                        !=================================================
                        do m = 1, l_nEq
                           do k = merge(m, 1, is_hermitian), l_nEq
                              do i = 1, N
                                 dHm(i+(m-1)*N, i+(k-1)*N) = dHm(i+(m-1)*N, i+(k-1)*N) + &
                                 Qv(i) * (dBm(i, p, m, k) + dBm(i, p, k, m)) * ginvv(i) * inv2mu
                              enddo
                           enddo
                        enddo

                        if (is_hermitian) then
                           ! Complete symmetric lower triangle for matmul
                           do i = 2, nH
                              do j = 1, i - 1
                                 dHm(i, j) = dHm(j, i)
                              enddo
                           enddo
                        endif

                     else
                        !=================================================
                        ! No B symmetry (real): build full dH directly
                        !=================================================
                        do m = 1, l_nEq
                           do k = 1, l_nEq
                              do j = 1, N
                                 dHm(1+(m-1)*N:m*N, j+(k-1)*N) = &
                                 dBm(1:N, p, m, k) * D1m(1:N, j) * (inv2mu * ginvv(j))
                              enddo
                           enddo
                        enddo

                        ! Qv-dependent diagonal-in-i contribution (real, no B-symmetry form)
                        do m = 1, l_nEq
                           do k = 1, l_nEq
                              do i = 1, N
                                 dHm(i+(m-1)*N, i+(k-1)*N) = dHm(i+(m-1)*N, i+(k-1)*N) + &
                                 Qv(i) * dBm(i, p, m, k) * ginvv(i) * inv2mu
                              enddo
                           enddo
                        enddo

                     endif
                     !=========================================================
                     ! Now compute dE(p, q)
                     !=========================================================
                     do q = 1, l_nEf

                        m = Eindv(q)
                        if (is_hermitian) then
                           tv  = matmul(dHm, Wm(:, m))
                        else
                           tcv = matmul(dHm, Wrcm(:, m))
                        endif

                        if (is_hermitian) then
                           if (is_G) then
                              t0 = 0.0_dp
                              do k = 1, l_nEq
                                 t0 = t0 + dot_product(Wm(1+(k-1)*N:k*N, m), &
                                    &                  Gm(1:N, k)*Wm(1+(k-1)*N:k*N, m))
                              enddo
                           else
                              t0 = dot_product(Wm(:, m), Wm(:, m))
                           endif
                           t1 = dot_product(Wm(:, m), tv)
                        elseif (is_arpack) then
                           if (is_G) then
                              c0 = (0.0_dp, 0.0_dp)
                              do k = 1, l_nEq
                                 c0 = c0 + dot_product(Wlcm(1+(k-1)*N:k*N, Etindv(q)), &
                                    &                   Gm(1:N, k)*Wrcm(1+(k-1)*N:k*N, m))
                              enddo
                           else
                              c0 = dot_product(Wlcm(:, Etindv(q)), Wrcm(:, m))
                           endif
                           c1 = dot_product(Wlcm(:, Etindv(q)), tcv(:))
                        else
                           if (is_G) then
                              c0 = (0.0_dp, 0.0_dp)
                              do k = 1, l_nEq
                                 c0 = c0 + dot_product(Wlcm(1+(k-1)*N:k*N, m), &
                                    &                   Gm(1:N, k)*Wrcm(1+(k-1)*N:k*N, m))
                              enddo
                           else
                              c0 = dot_product(Wlcm(:, m), Wrcm(:, m))
                           endif
                           c1 = dot_product(Wlcm(:, m), tcv(:))
                        endif
                        if (is_hermitian) then
                           if (present(dE)) dE(p, q) = dE(p, q) + t1/t0
                           if (present(dW)) dW(p, q) = dW(p, q) + 0.0_dp
                        else
                           if (present(dE)) dE(p, q) = dE(p, q) + real(c1/c0, dp)
                           if (present(dW)) dW(p, q) = dW(p, q) - 2*aimag(c1/c0)
                        endif
                     enddo
                  enddo
               endif
            endif

            if (is_dG) then
               !=========================================================
               ! Derivative contribution from the generalized matrix G:
               ! dE = <L|(dH - E*dG)|R> / <L|G|R>.
               ! This block adds only the -E*dG part.
               !=========================================================
               if (is_complex) then
                  do p = 1, l_ndE
                     do q = 1, l_nEf
                        m = Eindv(q)

                        do k = 1, l_nEq
                           tcv(1+(k-1)*N:k*N) = dGcm(1:N, p, k) * Wrcm(1+(k-1)*N:k*N, m)
                        enddo

                        if (is_hermitian) then
                           if (is_G) then
                              c0 = (0.0_dp, 0.0_dp)
                              do k = 1, l_nEq
                                 c0 = c0 + dot_product(Wrcm(1+(k-1)*N:k*N, m), &
                                    &                   Gcm(1:N, k)*Wrcm(1+(k-1)*N:k*N, m))
                              enddo
                              t0 = real(c0, dp)
                           else
                              t0 = real(dot_product(Wrcm(:, m), Wrcm(:, m)), dp)
                           endif
                           t1 = real(dot_product(Wrcm(:, m), tcv), dp)
                           if (present(dE)) dE(p, q) = dE(p, q) - Ev(m)*t1/t0
                           if (present(dW)) dW(p, q) = dW(p, q) + 0.0_dp
                        elseif (is_arpack) then
                           if (is_G) then
                              c0 = (0.0_dp, 0.0_dp)
                              do k = 1, l_nEq
                                 c0 = c0 + dot_product(Wlcm(1+(k-1)*N:k*N, Etindv(q)), &
                                    &                   Gcm(1:N, k)*Wrcm(1+(k-1)*N:k*N, m))
                              enddo
                           else
                              c0 = dot_product(Wlcm(:, Etindv(q)), Wrcm(:, m))
                           endif
                           c1 = -Ecv(m)*dot_product(Wlcm(:, Etindv(q)), tcv)/c0
                           if (present(dE)) dE(p, q) = dE(p, q) + real(c1, dp)
                           if (present(dW)) dW(p, q) = dW(p, q) - 2*aimag(c1)
                        else
                           if (is_G) then
                              c0 = (0.0_dp, 0.0_dp)
                              do k = 1, l_nEq
                                 c0 = c0 + dot_product(Wlcm(1+(k-1)*N:k*N, m), &
                                    &                   Gcm(1:N, k)*Wrcm(1+(k-1)*N:k*N, m))
                              enddo
                           else
                              c0 = dot_product(Wlcm(:, m), Wrcm(:, m))
                           endif
                           c1 = -Ecv(m)*dot_product(Wlcm(:, m), tcv)/c0
                           if (present(dE)) dE(p, q) = dE(p, q) + real(c1, dp)
                           if (present(dW)) dW(p, q) = dW(p, q) - 2*aimag(c1)
                        endif
                     enddo
                  enddo
               else
                  do p = 1, l_ndE
                     do q = 1, l_nEf
                        m = Eindv(q)

                        if (is_hermitian) then
                           do k = 1, l_nEq
                              tv(1+(k-1)*N:k*N) = dGm(1:N, p, k) * Wm(1+(k-1)*N:k*N, m)
                           enddo
                        else
                           do k = 1, l_nEq
                              tcv(1+(k-1)*N:k*N) = dGm(1:N, p, k) * Wrcm(1+(k-1)*N:k*N, m)
                           enddo
                        endif

                        if (is_hermitian) then
                           if (is_G) then
                              t0 = 0.0_dp
                              do k = 1, l_nEq
                                 t0 = t0 + dot_product(Wm(1+(k-1)*N:k*N, m), &
                                    &                  Gm(1:N, k)*Wm(1+(k-1)*N:k*N, m))
                              enddo
                           else
                              t0 = dot_product(Wm(:, m), Wm(:, m))
                           endif
                           t1 = dot_product(Wm(:, m), tv)
                           if (present(dE)) dE(p, q) = dE(p, q) - Ev(m)*t1/t0
                           if (present(dW)) dW(p, q) = dW(p, q) + 0.0_dp
                        elseif (is_arpack) then
                           if (is_G) then
                              c0 = (0.0_dp, 0.0_dp)
                              do k = 1, l_nEq
                                 c0 = c0 + dot_product(Wlcm(1+(k-1)*N:k*N, Etindv(q)), &
                                    &                   Gm(1:N, k)*Wrcm(1+(k-1)*N:k*N, m))
                              enddo
                           else
                              c0 = dot_product(Wlcm(:, Etindv(q)), Wrcm(:, m))
                           endif
                           c1 = -Ecv(m)*dot_product(Wlcm(:, Etindv(q)), tcv)/c0
                           if (present(dE)) dE(p, q) = dE(p, q) + real(c1, dp)
                           if (present(dW)) dW(p, q) = dW(p, q) - 2*aimag(c1)
                        else
                           if (is_G) then
                              c0 = (0.0_dp, 0.0_dp)
                              do k = 1, l_nEq
                                 c0 = c0 + dot_product(Wlcm(1+(k-1)*N:k*N, m), &
                                    &                   Gm(1:N, k)*Wrcm(1+(k-1)*N:k*N, m))
                              enddo
                           else
                              c0 = dot_product(Wlcm(:, m), Wrcm(:, m))
                           endif
                           c1 = -Ecv(m)*dot_product(Wlcm(:, m), tcv)/c0
                           if (present(dE)) dE(p, q) = dE(p, q) + real(c1, dp)
                           if (present(dW)) dW(p, q) = dW(p, q) - 2*aimag(c1)
                        endif
                     enddo
                  enddo
               endif
            endif
         endif
      end subroutine compute_energy_derivatives

      !=================================================
      ! Internal procedure (see the call site for what this phase does);
      ! shares qsolver state by host association.
      !=================================================
      subroutine assemble_outputs()
         complex(dp), allocatable :: qwcv(:)   ! complex quadrature weights (weighted bra)

         if (present(E)) then
            if (size(E) < l_nEf) then
               call err%push_err(ERR_INPUT, 'Size of E array is too small')
               return
            endif
            if (is_hermitian) then 
               do i = 1, l_nEf
                  E(i) = Ev(Eindv(i))
               enddo
            else
               do i = 1, l_nEf
                  E(i) = real(Ecv(Eindv(i)), dp)
               enddo
            endif
         endif

         !=========================================================
         ! Return obtained widths 
         !=========================================================
         if (present(W)) then 
            if (size(W) < l_nEf) then
               call err%push_err(ERR_INPUT, 'Size of W array is too small')
               return
            endif
            if (is_hermitian) then 
               W(1:l_nEf) = 0.0_dp
            else
               do i = 1, l_nEf
                  W(i) = -2*aimag(Ecv(Eindv(i)))
               enddo
            endif
         endif

         !=================================================
         ! Calculation of quadrature weights. The combined mesh weights are
         ! stored in xv%qw (xv is the physical mesh object); its storage was
         ! allocated together with the nodes by xv%build.
         !=================================================
         if (is_WF .or. is_lWF .or. present(qw) .or. present(qwc)) then
            ! LGL weights on mapped interval(s).
            if (is_ecs) then
               ! Two-subinterval case (real + complex-scaled part)
               if (N1 /= 0) then
                  xv%qw(1:N1-1) = grid_inner%qw(1:N1-1)
                  xv%qw(N1)     = grid_inner%qw(N1) + grid_outer%qw(0)
               endif
               xv%qw(N1+1:N) = grid_outer%qw(1:N2-1)
            else
               ! Single-interval case
               xv%qw(1:N) = g0%qw(1:N)
            endif

         endif

         if (is_WF .or. is_lWF) then
            !=========================================================
            ! Constants from LGL weights:
            !   symmetric single-interval case uses one constant t1
            !=========================================================
            if (is_hermitian .and. .not. is_D_sym) then
               call err%push_err(ERR_INTERNAL, 'Hermitian run with non-symmetric D matrices is not supported')
               return
            endif

            if (is_D_sym .and. is_ecs) then
               call err%push_err(ERR_INPUT, 'Symmetric differentiation matrix with split interval is not implemented')
               return
            endif

            if (is_D_sym) then
               ! LGL endpoint normalization constant on the mapped interval:
               ! (xmax-xmin)/(M(M-1)) for a Lobatto grid of M = g0%hi-g0%lo+1 nodes.
               t1 = (xmax - xmin)/(dble(g0%hi - g0%lo)*(g0%hi - g0%lo + 1))
            endif

            !=========================================================
            ! Normalize  wavefunctions
            !=========================================================
            if (is_hermitian) then
               if (is_complex_WF) then 
                  do i = 1, l_nEf 
                     m = Eindv(i)
                     t0 = 0.0_dp
                     !=========================================================
                     ! Normalize complex wavefunction such that <psi_i|psi_i> = 1
                     !=========================================================
                     do k = 1, l_nEq
                        associate (WFki => Wrcm(1+(k-1)*N:k*N, m))
                           t0 = t0 + real(dot_product(WFki, WFki), dp)
                        end associate
                     enddo
                     if (t0*t1 <= tiny(1.0_dp)) then
                        call err%push_err(ERR_INTERNAL, 'Zero wavefunction norm')
                        return
                     endif
                     Wrcm(1:nH, m) = Wrcm(1:nH, m)/sqrt(t0*t1)
                  enddo
               else
                  do i = 1, l_nEf 
                     q = Eindv(i)
                     t0 = 0.0_dp
                     !=========================================================
                     ! Normalize real wavefunction such that <psi_i|psi_i> = 1
                     !=========================================================
                     do k = 1, l_nEq
                        associate (WFki => Wm(1+(k-1)*N:k*N, q))
                           t0 = t0 + sum(WFki**2)
                        end associate
                     enddo
                     if (t0*t1 <= tiny(1.0_dp)) then
                        call err%push_err(ERR_INTERNAL, 'Zero wavefunction norm')
                        return
                     endif
                     Wm(1:nH, q) = Wm(1:nH, q)/sqrt(t0*t1)
                  enddo
               endif
            else
               if (is_ecs) then 
                  do i = 1, l_nEf 
                     q = Eindv(i)
                     c0 = (0.0_dp, 0.0_dp)
                     t2 = grid_outer%qw(0)
                     !=========================================================
                     ! Exterior complex scaling (ECS) case: rc is present.
                     !
                     ! For r > rc the integration contour is rotated into the
                     ! complex plane, so the quadrature weights become complex
                     ! in the exterior region.
                     ! Because of this, the usual Hermitian norm
                     !
                     !     <rpsi_i | rpsi_i> = integral(conjg(rpsi_i)*rpsi_i*dr)
                     !
                     ! is, in general, complex and therefore cannot be used for
                     ! normalization to unity.
                     !
                     ! Instead, for a complex-symmetric problem (is_V_sym = .true.
                     ! is_B_sym = .true.) left eigenfunction is a conjugate of the
                     ! right eigenfunction so norm becomes the so-called c-product:
                     !
                     !     integral(rpsi_i(r)**2*dr) = 1
                     !
                     ! For complex scaling problems qsolver always norms wavefunction 
                     ! this way.
                     !
                     !=========================================================
                     do k = 1, l_nEq
                        associate (WFki => Wrcm(1+(k-1)*N: N1 + (k-1)*N, q))
                           c0 = c0 + sum(xv%qw(1:N1) * WFki**2)
                           c0 = c0 + (drdyc - 1.0_dp) * t2 * WFki(N1)**2
                        end associate
                        associate (WFki => Wrcm(N1 + 1 + (k-1)*N: k*N, q))
                           c0 = c0 + sum(xv%qw(N1+1:N) * WFki**2) * drdyc
                        end associate
                     enddo
                     if (abs(c0) <= tiny(1.0_dp)) then
                        call err%push_err(ERR_INTERNAL, 'Zero wavefunction norm')
                        return
                     endif
                     Wrcm(1:nH, q) = Wrcm(1:nH, q) / sqrt(c0)

                     if (is_lWF) then
                        p = Eindv(i)
                        if (is_arpack) p = Etindv(i)
                        c0 = (0.0_dp, 0.0_dp)
                        !=========================================================
                        ! Exterior complex scaling case (rc is present).
                        ! Normalize complex left eigenfunction such that <lpsi_i|rpsi_i> = 1
                        !=========================================================
                        do k = 1, l_nEq
                           associate (lWFki => Wlcm(1+(k-1)*N: N1 + (k-1)*N, p), & 
                                   &  rWFki => Wrcm(1+(k-1)*N: N1 + (k-1)*N, q))
                              c0 = c0 + dot_product(lWFki, xv%qw(1:N1)*rWFki)
                              c0 = c0 + (drdyc - 1.0_dp) * t2 * conjg(lWFki(N1)) * rWFki(N1)
                           end associate
                           associate (lWFki => Wlcm(N1 + 1 + (k-1)*N: k*N, p), &
                                   &  rWFki => Wrcm(N1 + 1 + (k-1)*N: k*N, q))
                              c0 = c0 + dot_product(lWFki, xv%qw(N1+1:N)*rWFki)*drdyc
                           end associate
                        enddo
                        if (abs(c0) <= tiny(1.0_dp)) then
                           call err%push_err(ERR_INTERNAL, 'Zero left/right wavefunction overlap')
                           return
                        endif
                        Wlcm(1:nH, p) = Wlcm(1:nH, p)/conjg(c0)

                     endif
                  enddo
               else
                  do i = 1, l_nEf 
                     q = Eindv(i)
                     !=========================================================
                     ! Normalize complex right wavefunction.
                     !
                     ! For complex scaling (phi present): use c-norm.
                     ! Otherwise: use Hermitian norm with conjg(psi)*psi.
                     !=========================================================
                     if (is_quasi) then
                        c0 = (0.0_dp, 0.0_dp)
                        do k = 1, l_nEq
                           associate (rWFki => Wrcm(1+(k-1)*N:k*N, q))
                              if (is_D_sym) then
                                 c0 = c0 + t1*sum(rWFki**2)
                              else
                                 c0 = c0 + sum(xv%qw(1:N)*rWFki**2)
                              endif
                           end associate
                        enddo
                        c0 = c0*drdyc
                        if (abs(c0) <= tiny(1.0_dp)) then
                           call err%push_err(ERR_INTERNAL, 'Zero wavefunction norm')
                           return
                        endif
                        Wrcm(1:nH, q) = Wrcm(1:nH, q)/sqrt(c0)
                     else
                        t0 = 0.0_dp
                        do k = 1, l_nEq
                           associate (rWFki => Wrcm(1+(k-1)*N:k*N, q))
                              if (is_D_sym) then
                                 t0 = t0 + t1*real(dot_product(rWFki, rWFki), dp)
                              else
                                 t0 = t0 + real(dot_product(rWFki, xv%qw(1:N)*rWFki), dp)
                              endif
                           end associate
                        enddo
                        if (t0 <= tiny(1.0_dp)) then
                           call err%push_err(ERR_INTERNAL, 'Zero wavefunction norm')
                           return
                        endif
                        Wrcm(1:nH, q) = Wrcm(1:nH, q)/sqrt(t0)
                     endif

                     if (is_lWF) then
                        p = Eindv(i)
                        if (is_arpack) p = Etindv(i)
                        c0 = (0.0_dp, 0.0_dp)
                        !=========================================================
                        ! Normalize complex left wavefunction such that <lpsi_i|rpsi_i> = 1
                        !=========================================================
                        do k = 1, l_nEq
                           associate (lWFki => Wlcm(1+(k-1)*N:k*N, p), & 
                                      rWFki => Wrcm(1+(k-1)*N:k*N, q))
                              c0 = c0 + dot_product(lWFki, xv%qw(1:N)*rWFki)
                           end associate
                        enddo
                        if (is_quasi) then  
                           c0 = c0*drdyc
                           if (abs(c0) <= tiny(1.0_dp)) then
                              call err%push_err(ERR_INTERNAL, 'Zero left/right wavefunction overlap')
                              return
                           endif
                           Wlcm(1:nH, p) = Wlcm(1:nH, p)/conjg(c0)
                        else
                           if (abs(c0) <= tiny(1.0_dp)) then
                              call err%push_err(ERR_INTERNAL, 'Zero left/right wavefunction overlap')
                              return
                           endif
                           Wlcm(1:nH, p) = Wlcm(1:nH, p)/conjg(c0)
                        endif

                     endif
                  enddo
               endif
            endif

            !=========================================================
            ! Undo differential matrices symmetrization
            !=========================================================
            if (is_D_sym) then
               if (is_complex_WF) then 
                  do i = 1, l_nEf 
                     q = Eindv(i)
                     do k = 1, l_nEq
                        Wrcm(1+(k-1)*N:k*N, q) = Wrcm(1+(k-1)*N:k*N, q) / bary(1:N)
                     enddo
                  enddo
               else
                  do i = 1, l_nEf 
                     q = Eindv(i)
                     do k = 1, l_nEq
                        Wm(1+(k-1)*N:k*N, q) = Wm(1+(k-1)*N:k*N, q) / bary(1:N)
                     enddo
                  enddo
               endif

               if (is_lWF) then 
                  do i = 1, l_nEf 
                     q = Eindv(i)
                     if (is_arpack) q = Etindv(i)
                     do k = 1, l_nEq
                        Wlcm(1+(k-1)*N:k*N, q) = Wlcm(1+(k-1)*N:k*N, q) * bary(1:N)
                     enddo
                  enddo
               endif
            endif

            !=========================================================
            ! Extract original unmapped wavefunctions
            !=========================================================
            if (is_complex_WF) then 
               do i = 1, l_nEf 
                  q = Eindv(i)
                  do k = 1, l_nEq
                     Wrcm(1+(k-1)*N:k*N, q) = Wrcm(1+(k-1)*N:k*N, q) * sqrt(ginvv(1:N))
                  enddo
               enddo
            else
               do i = 1, l_nEf 
                  q = Eindv(i)
                  do k = 1, l_nEq
                     Wm(1+(k-1)*N:k*N, q) = Wm(1+(k-1)*N:k*N, q) * sqrt(ginvv(1:N))
                  enddo
               enddo
            endif

            if (is_lWF) then 
               do i = 1, l_nEf 
                  q = Eindv(i)
                  if (is_arpack) q = Etindv(i)
                  do k = 1, l_nEq
                     Wlcm(1+(k-1)*N:k*N, q) = Wlcm(1+(k-1)*N:k*N, q) / sqrt(ginvv(1:N)) 
                  enddo
               enddo
            endif
         endif
         
         !=========================================================
         ! Return obtained real wavefunctions
         !=========================================================
         if (present(WF)) then
            if (size(WF, 1) < N .or. size(WF, 2) < l_nEq .or. size(WF, 3) < l_nEf) then
               call err%push_err(ERR_INPUT, 'WF array is too small')
               return
            endif
            if (is_complex_WF) then 
               do i = 1, l_nEf 
                  do k = 1, l_nEq
                     WF(1:N, k, i) = real(Wrcm(1+(k-1)*N:k*N, Eindv(i)), dp)
                  enddo
               enddo
            else
               do i = 1, l_nEf 
                  do k = 1, l_nEq
                     WF(1:N, k, i) = Wm(1+(k-1)*N:k*N, Eindv(i))
                  enddo
               enddo
            endif
         endif

         !=========================================================
         ! Return obtained real wavefunctions for one equation case
         !=========================================================
         if (present(WF1)) then
            if (size(WF1, 1) < N .or. size(WF1, 2) < l_nEf) then
               call err%push_err(ERR_INPUT, 'WF1 array is too small')
               return
            endif
            if (is_complex_WF) then 
               do i = 1, l_nEf 
                  WF1(1:N, i) = real(Wrcm(1:N, Eindv(i)), dp)
               enddo
            else
               do i = 1, l_nEf 
                  WF1(1:N, i) = Wm(1:N, Eindv(i))
               enddo
            endif
         endif

         !=========================================================
         ! Return obtained complex wavefunctions
         !=========================================================
         if (present(WFc)) then
            if (size(WFc, 1) < N .or. size(WFc, 2) < l_nEq .or. size(WFc, 3) < l_nEf) then
               call err%push_err(ERR_INPUT, 'WFc array is too small')
               return
            endif
            if (is_complex_WF) then 
               do i = 1, l_nEf 
                  do k = 1, l_nEq
                     WFc(1:N, k, i) = Wrcm(1+(k-1)*N:k*N, Eindv(i))
                  enddo
               enddo
            else
               do i = 1, l_nEf 
                  do k = 1, l_nEq
                     WFc(1:N, k, i) = cmplx(Wm(1+(k-1)*N:k*N, Eindv(i)), 0.0_dp, dp)
                  enddo
               enddo
            endif
         endif

         !=========================================================
         ! Return obtained complex wavefunctions for one equation case
         !=========================================================
         if (present(WF1c)) then
            if (size(WF1c, 1) < N .or. size(WF1c, 2) < l_nEf) then
               call err%push_err(ERR_INPUT, 'WF1c array is too small')
               return
            endif
            if (is_complex_WF) then 
               do i = 1, l_nEf 
                  WF1c(1:N, i) = Wrcm(1:N, Eindv(i))
               enddo
            else
               do i = 1, l_nEf 
                  WF1c(1:N, i) = cmplx(Wm(1:N, Eindv(i)), 0.0_dp, dp)
               enddo
            endif
         endif

         !=========================================================
         ! Return weighted real bra vectors
         !=========================================================
         if (present(braW)) then
            if (size(braW, 1) < N .or. size(braW, 2) < l_nEq .or. size(braW, 3) < l_nEf) then
               call err%push_err(ERR_INPUT, 'braW array is too small')
               return
            endif
            do i = 1, l_nEf
               do k = 1, l_nEq
                  braW(1:N, k, i) = Wm(1+(k-1)*N:k*N, Eindv(i)) * xv%qw(1:N) * gv(1:N)
               enddo
            enddo
         endif

         !=========================================================
         ! Return weighted real bra vectors in case of single equation
         !=========================================================
         if (present(braW1)) then
            if (size(braW1, 1) < N .or. size(braW1, 2) < l_nEf) then
               call err%push_err(ERR_INPUT, 'braW1 array is too small')
               return
            endif
            do i = 1, l_nEf
               braW1(1:N, i) = Wm(1:N, Eindv(i)) * xv%qw(1:N) * gv(1:N)
            enddo
         endif

         !=========================================================
         ! Build complex quadrature weights for weighted bra vectors
         !=========================================================
         if (present(braWc) .or. present(braW1c)) then
            allocate(qwcv(1:N), stat=istat)
            if (alloc_failed(istat, 'Allocation failed for qwcv', err)) return

            qwcv(1:N) = cmplx(xv%qw(1:N)*gv(1:N), 0.0_dp, dp)
            if (is_quasi) then
               if (is_ecs) then
                  t1 = grid_outer%qw(0) * gv(N1)
                  qwcv(N1) = qwcv(N1) + t1*(drdyc - 1.0_dp)
                  qwcv(N1 + 1:N) = qwcv(N1 + 1:N)*drdyc
               else
                  qwcv(1:N) = qwcv(1:N)*drdyc
               endif
            endif
         endif

         !=========================================================
         ! Return obtained weighted complex bra vectors
         !=========================================================
         if (present(braWc)) then
            if (size(braWc, 1) < N .or. size(braWc, 2) < l_nEq .or. size(braWc, 3) < l_nEf) then
               call err%push_err(ERR_INPUT, 'braW array is too small')
               return
            endif
            if (.not. is_hermitian) then 
               do i = 1, l_nEf 
                  p = Eindv(i)
                  if (is_arpack) p = Etindv(i)
                  do k = 1, l_nEq
                     braWc(1:N, k, i) = Wlcm(1+(k-1)*N:k*N, p) * conjg(qwcv(1:N))
                  enddo
               enddo
            else
               if (is_complex_WF) then 
                  do i = 1, l_nEf 
                     do k = 1, l_nEq
                        braWc(1:N, k, i) = Wrcm(1+(k-1)*N:k*N, Eindv(i)) * conjg(qwcv(1:N))
                     enddo
                  enddo
               else
                  do i = 1, l_nEf 
                     do k = 1, l_nEq
                        braWc(1:N, k, i) = cmplx(Wm(1+(k-1)*N:k*N, Eindv(i)), 0.0_dp, dp) * conjg(qwcv(1:N))
                     enddo
                  enddo
               endif
            endif
         endif

         !=========================================================
         ! Return obtained weighted complex bra vectors in case of
         ! single equation
         !=========================================================
         if (present(braW1c)) then
            if (size(braW1c, 1) < N .or. size(braW1c, 2) < l_nEf) then
               call err%push_err(ERR_INPUT, 'braW array is too small')
               return
            endif
            if (.not. is_hermitian) then 
               do i = 1, l_nEf 
                  p = Eindv(i)
                  if (is_arpack) p = Etindv(i)
                  braW1c(1:N, i) = Wlcm(1:N, p) * conjg(qwcv(1:N))
               enddo
            else
               if (is_complex_WF) then 
                  do i = 1, l_nEf 
                     braW1c(1:N, i) = Wrcm(1:N, Eindv(i)) * conjg(qwcv(1:N))
                  enddo
               else
                  do i = 1, l_nEf 
                     braW1c(1:N, i) = cmplx(Wm(1:N, Eindv(i)), 0.0_dp, dp) * conjg(qwcv(1:N))
                  enddo
               endif
            endif
         endif

         !=================================================
         ! Return real grid points in r 
         ! If grid is complex then only real part is returned
         !=================================================
         if (present(rv)) then 
            if (size(rv) < N) then
               call err%push_err(ERR_INPUT, 'array rv is too small')
               return
            endif
            if (is_quasi) then 
               rv(1:N) = real(l_rcv(1:N), dp)
            else
               rv(1:N) = yv(1:N)
            endif
         endif

         !=================================================
         ! Return complex grid points in r 
         ! If grid is complex then only real part is returned
         !=================================================
         if (present(rvc)) then 
            if (size(rvc) < N) then
               call err%push_err(ERR_INPUT, 'array rvc is too small')
               return
            endif
            if (is_quasi) then 
               rvc(1:N) = l_rcv(1:N)
            else
               rvc(1:N) = yv(1:N)
            endif
         endif

         !=================================================
         ! Return real quadrature weights (mesh weights times the Jacobian g;
         ! xv%qw itself is left as the pure quadrature weights).
         !=================================================
         if (present(qw)) then
            if (size(qw) < N) then
               call err%push_err(ERR_INPUT, 'array qw is too small')
               return
            endif
            qw(1:N) = xv%qw(1:N)*gv(1:N)
         endif

         !=================================================
         ! Return complex quadrature weights
         !=================================================
         if (present(qwc)) then
            if (size(qwc) < N) then
               call err%push_err(ERR_INPUT, 'array qwc is too small')
               return
            endif
            qwc(1:N) = xv%qw(1:N)*gv(1:N)
            if (is_quasi) then 
               if (is_ecs) then 
                  t1 = grid_outer%qw(0) * gv(N1)
                  qwc(N1) = qwc(N1) + t1*(drdyc - 1.0_dp)
                  qwc(N1 + 1: N) = qwc(N1 + 1: N)*drdyc
               else
                  qwc(1:N) = qwc(1:N)*drdyc
               endif
            endif
         endif

         !=================================================
         ! Return number of eigenvalues found
         !=================================================
         if (present(nEf)) nEf = l_nEf
      end subroutine assemble_outputs
   end subroutine qsolver

   !=======================================================================
   ! Initialize the error handler.
   ! Must be called at the start of qsolver before any other err methods.
   ! Sets error mode and resets all counters and message buffer.
   !   err_mode = ERRMODE_PRINTALL (1) : print errors and warnings, continue
   !   err_mode = ERRMODE_NOWARN   (2) : silent, continue
   !   err_mode = ERRMODE_STOP     (3) : print and halt immediately on error
   !=======================================================================
   subroutine qerr_init(self, err_mode)
       class(err_type), intent(out) :: self
       integer,         intent(in)  :: err_mode
       self%mode    = err_mode
       self%errcode = 0
       self%nerr    = 0
       self%nwarn   = 0
       self%nmsg    = 0
       self%is_truncated = .false.
   end subroutine qerr_init
   
   !=======================================================================
   ! Reset all error state while preserving the current mode.
   ! Useful if the same err object is reused across multiple calls.
   !=======================================================================
   subroutine qerr_clear(self)
       class(err_type), intent(inout) :: self
       self%errcode = 0
       self%nerr    = 0
       self%nwarn   = 0
       self%nmsg    = 0
       self%is_truncated = .false.
   end subroutine qerr_clear
   
   !=======================================================================
   ! Internal helper. Appends a prefixed message to the message buffer
   ! and prints it unless mode = ERRMODE_NOWARN.
   ! If the buffer is full, the last slot is overwritten with a
   ! truncation notice, which is printed only once.
   !=======================================================================
   subroutine qerr_push_msg(self, prefix, msg)
       class(err_type),  intent(inout) :: self
       character(len=*), intent(in)    :: prefix, msg
       character(len=256)              :: line
       character(len=*), parameter     :: trunc_msg = '... (further messages suppressed)'
   
       line = trim(prefix) // trim(msg)
   
       if (self%nmsg < MAXERRMSG) then
           self%nmsg = self%nmsg + 1
           self%messages(self%nmsg) = line
           if (self%mode /= ERRMODE_NOWARN) write(*, '(A)') trim(line)
       elseif (.not. self%is_truncated) then
           self%messages(MAXERRMSG) = trunc_msg
           self%is_truncated = .true.
           if (self%mode /= ERRMODE_NOWARN) write(*, '(A)') trunc_msg
       end if
   end subroutine qerr_push_msg
   
   !=======================================================================
   ! Register an error with category code and message.
   ! Increments nerr and adds code to errcode (digit capped at 9).
   ! In mode ERRMODE_STOP halts the program immediately after printing.
   ! In modes ERRMODE_PRINTALL and ERRMODE_NOWARN accumulates and returns.
   !   code : one of ERR_INPUT, ERR_ALLOC, ERR_INTERNAL,
   !                 ERR_LAPACK, ERR_ARPACK
   !=======================================================================
   subroutine qerr_push_err(self, code, msg)
       class(err_type),  intent(inout) :: self
       integer,          intent(in)    :: code
       character(len=*), intent(in)    :: msg
   
       self%nerr = self%nerr + 1
   
       ! Saturate the relevant decimal digit at 9
       if (mod(self%errcode / code, 10) < 9) then
           self%errcode = self%errcode + code
       else
           call self%push_msg('WARNING: ', &
               'errcode digit saturated at 9 for category: ' // trim(msg))
       end if
   
       call self%push_msg('ERROR: ', msg)
   
       ! Mode 3: halt immediately on first error
       if (self%mode == ERRMODE_STOP) stop 1
   end subroutine qerr_push_err
   
   !=======================================================================
   ! Register a warning message.
   ! Increments nwarn and appends message to buffer.
   ! Printed unless mode = ERRMODE_NOWARN.
   ! Does not affect errcode or nerr.
   !=======================================================================
   subroutine qerr_push_warning(self, msg)
       class(err_type),  intent(inout) :: self
       character(len=*), intent(in)    :: msg
   
       self%nwarn = self%nwarn + 1
       call self%push_msg('WARNING: ', msg)
   end subroutine qerr_push_warning   
   !=======================================================================
   ! Returns .true. if at least one error was registered via push_err.
   ! Used to trigger early exit (goto 100) in qsolver.
   !=======================================================================
   function qerr_has_error(self) result(has)
       class(err_type), intent(in) :: self
       logical :: has
       has = (self%nerr > 0)
   end function qerr_has_error
   
   !=======================================================================
   ! Called once at the end of qsolver.
   ! Fills optional output arguments errcode and errmessage with all
   ! accumulated results so the caller can inspect them.
   ! Never stops — halting is push_err's responsibility (mode 3).
   ! In mode 3 this routine is never reached after an error.
   ! errmessage is allocatable: caller needs no size declaration,
   ! the string is allocated to exactly the required length.
   !=======================================================================
   subroutine qerr_finalize(self, errcode, errmessage)
       class(err_type),                         intent(in)  :: self
       integer,                       optional, intent(out) :: errcode
       character(len=:), allocatable, optional, intent(out) :: errmessage
       integer :: i
   
       if (present(errcode)) errcode = self%errcode
   
       if (present(errmessage)) then
           errmessage = ''
           do i = 1, self%nmsg
               if (i > 1) errmessage = errmessage // '; '
               errmessage = errmessage // trim(self%messages(i))
           end do
       end if
   end subroutine qerr_finalize

   !=======================================================================
   ! Allocation-status check. Returns .true. (after recording an ERR_ALLOC
   ! with the given message) when istat indicates a failed allocation, so a
   ! whole allocate guard collapses to one line:
   !     allocate(X(...), stat=istat)
   !     if (alloc_failed(istat, 'Allocation failed for X', err)) goto 100
   !=======================================================================
   logical function alloc_failed(istat, msg, err)
      integer,          intent(in)    :: istat
      character(len=*), intent(in)    :: msg
      class(err_type),  intent(inout) :: err

      alloc_failed = (istat /= 0)
      if (alloc_failed) call err%push_err(ERR_ALLOC, msg)
   end function alloc_failed

   !================================================================
   ! One or two-parametric mapping initialization
   ! Initialize mapping type and mapping parameters.
   !================================================================
   subroutine init_mapping(self, name, beta, ye)
      class(mapping_type), intent(inout) :: self
      character(*), optional, intent(in) :: name
      real(dp),     optional, intent(in) :: beta, ye

      if (.not. present(name)) then 
         self%mn = 0
         return
      else 
         select case (name)
            case('none')
               self%mn = 0
               return
            case ('atan')
               self%mn = 1
            case ('rat')
               self%mn = 2
            case ('tanh')
               self%mn = 3
            case ('dexp')
               self%mn = 4
            case ('exp')
               self%mn = 5
            case default 
               if (associated(self%err)) then
                  call self%err%push_err(ERR_INPUT, 'Unknown mapping type')
               endif
               self%mn = 0
               self%beta = 0.0_dp
               self%ye = 0.0_dp
               return
         end select
      endif

      if (.not. present(beta)) then
         if (associated(self%err)) then
            call self%err%push_err(ERR_INPUT, 'At least parameter beta must be specified if mapping is used')
         endif
         self%mn = 0
         self%beta = 0.0_dp
         self%ye = 0.0_dp
         return

      else if (present(ye)) then
         ! two-parameter mapping
         if (self%mn == 1 .or. self%mn == 3 .or. self%mn == 4) then
            ! mappings of the form  x = f(beta*(y/ye - 1)) implemented as f((beta/ye)*(y - ye))
            self%beta = beta/ye
            self%ye   = ye
         else
            ! mappings of the form x = f((y/ye)**beta)
            self%beta = beta
            self%ye   = ye
         endif

      else
         ! one-parameter mapping
         self%beta = beta
         if (self%mn == 2 .or. self%mn == 5) then
            ! 'rat' and 'exp' mappings use y/ye, so default scale is ye = 1
            self%ye = 1.0_dp
         else
            self%ye = 0.0_dp
         endif
      endif

   end subroutine

   !===================================================
   ! Functions initialization subroutines
   !===================================================
   !================================================================
   ! Version for f(x)
   ! Bind real function with signature f(x).
   !================================================================
   subroutine init_func1(self, f)
      class(func_type), intent(inout) :: self
      procedure(func1) :: f

      self%p1    => f
      self%n_par =  1
      self%is_associated = .true. 

   end subroutine init_func1

   !================================================================
   ! Version for f(x,k)
   ! Bind real function with signature f(x, k).
   !================================================================
   subroutine init_func2(self, f)
      class(func_type), intent(inout) :: self
      procedure(func2) :: f

      self%p2    => f
      self%n_par =  2
      self%is_associated = .true. 

   end subroutine init_func2

   !================================================================
   ! Version for f(x,i,j)
   ! Bind real function with signature f(x, i, j).
   !================================================================
   subroutine init_func3(self, f)
      class(func_type), intent(inout) :: self
      procedure(func3) :: f

      self%p3    => f
      self%n_par =  3
      self%is_associated = .true. 

   end subroutine init_func3

   !================================================================
   ! Bind real function with signature f(x, i, j, k).
   !================================================================
   subroutine init_func4(self, f)
      class(func_type), intent(inout) :: self
      procedure(func4) :: f

      self%p4    => f
      self%n_par =  4
      self%is_associated = .true. 

   end subroutine init_func4

   !================================================================
   ! Version for f(x, df/dp)
   ! Bind real derivative subroutine with signature f(x, df/dp).
   !================================================================
   subroutine init_vfunc1(self, f)
      class(func_type), intent(inout) :: self
      procedure(vfunc1) :: f

      self%v1    => f
      self%n_par =  1
      self%is_associated = .true. 

   end subroutine init_vfunc1

   !================================================================
   ! Version for f(x,i, df/dp)
   ! Bind real derivative subroutine with signature f(x, i, df/dp).
   !================================================================
   subroutine init_vfunc2(self, f)
      class(func_type), intent(inout) :: self
      procedure(vfunc2) :: f

      self%v2    => f
      self%n_par =  2
      self%is_associated = .true. 

   end subroutine init_vfunc2

   !================================================================
   ! Version for f(x,i,j, df/dp)
   ! Bind real derivative subroutine with signature f(x, i, j, df/dp).
   !================================================================
   subroutine init_vfunc3(self, f)
      class(func_type), intent(inout) :: self
      procedure(vfunc3) :: f

      self%v3    => f
      self%n_par =  3
      self%is_associated = .true. 

   end subroutine init_vfunc3

   !================================================================
   ! Reset real-function binding.
   !================================================================
   subroutine init_func(self)
      class(func_type), intent(inout) :: self

      self%is_associated = .false. 

   end subroutine init_func

   !================================================================
   ! Version for f(x)
   ! Bind complex function with signature f(x).
   !================================================================
   subroutine init_cfunc1(self, f)
      class(cfunc_type), intent(inout) :: self
      procedure(cfunc1) :: f

      self%p1    => f
      self%n_par =  1
      self%is_associated = .true. 

   end subroutine init_cfunc1

   !================================================================
   ! Version for f(x,k)
   ! Bind complex function with signature f(x, k).
   !================================================================
   subroutine init_cfunc2(self, f)
      class(cfunc_type), intent(inout) :: self
      procedure(cfunc2) :: f

      self%p2    => f
      self%n_par =  2
      self%is_associated = .true. 

   end subroutine init_cfunc2

   !================================================================
   ! Version for f(x,i,j)
   ! Bind complex function with signature f(x, i, j).
   !================================================================
   subroutine init_cfunc3(self, f)
      class(cfunc_type), intent(inout) :: self
      procedure(cfunc3) :: f

      self%p3    => f
      self%n_par =  3
      self%is_associated = .true. 

   end subroutine init_cfunc3

   !================================================================
   ! Version for f(x,i,j,k)
   ! Bind complex function with signature f(x, i, j, k).
   !================================================================
   subroutine init_cfunc4(self, f)
      class(cfunc_type), intent(inout) :: self
      procedure(cfunc4) :: f

      self%p4 => f
      self%n_par = 4
      self%is_associated = .true. 

   end subroutine init_cfunc4

   !================================================================
   ! Version for f(x, df/dp)
   ! Bind complex derivative subroutine with signature f(x, df/dp).
   !================================================================
   subroutine init_cvfunc1(self, f)
      class(cfunc_type), intent(inout) :: self
      procedure(cvfunc1) :: f

      self%v1    => f
      self%n_par =  1
      self%is_associated = .true. 

   end subroutine init_cvfunc1

   !================================================================
   ! Version for f(x,i, df/dp)
   ! Bind complex derivative subroutine with signature f(x, i, df/dp).
   !================================================================
   subroutine init_cvfunc2(self, f)
      class(cfunc_type), intent(inout) :: self
      procedure(cvfunc2) :: f

      self%v2    => f
      self%n_par =  2
      self%is_associated = .true. 

   end subroutine init_cvfunc2

   !================================================================
   ! Version for f(x,i,j, df/dp)
   ! Bind complex derivative subroutine with signature f(x, i, j, df/dp).
   !================================================================
   subroutine init_cvfunc3(self, f)
      class(cfunc_type), intent(inout) :: self
      procedure(cvfunc3) :: f

      self%v3    => f
      self%n_par =  3
      self%is_associated = .true. 

   end subroutine init_cvfunc3

   !================================================================
   ! Reset complex-function binding.
   !================================================================
   subroutine init_cfunc(self)
      class(cfunc_type), intent(inout) :: self

      self%is_associated = .false. 

   end subroutine init_cfunc

   !================================================================
   ! Subroutines for setting function parameters
   !================================================================
   ! Set 2-parameter function-index parameter.
   !================================================================
   subroutine set_func_par1(self, i)
      class(func_type), intent(inout) :: self
      integer, intent(in) :: i

      self%i = i

   end subroutine set_func_par1

   !================================================================
   ! Set 3-parameter function-index parameters.
   !================================================================
   subroutine set_func_par2(self, i, j)
      class(func_type), intent(inout) :: self
      integer, intent(in) :: i, j

      self%i = i
      self%j = j

   end subroutine set_func_par2

   !================================================================
   ! Set 4-parameter function-index parameters.
   !================================================================
   subroutine set_func_par3(self, i, j, k)
      class(func_type), intent(inout) :: self
      integer, intent(in) :: i, j, k

      self%i = i
      self%j = j
      self%k = k

   end subroutine set_func_par3

   !================================================================
   ! Set 2-parameter complex function-index parameter.
   !================================================================
   subroutine set_cfunc_par1(self, i)
      class(cfunc_type), intent(inout) :: self
      integer, intent(in) :: i

      self%i = i

   end subroutine set_cfunc_par1

   !================================================================
   ! Set 3-parameter complex function-index parameters.
   !================================================================
   subroutine set_cfunc_par2(self, i, j)
      class(cfunc_type), intent(inout) :: self
      integer, intent(in) :: i, j

      self%i = i
      self%j = j

   end subroutine set_cfunc_par2

   !================================================================
   ! Set 4-parameter complex function-index parameters.
   !================================================================
   subroutine set_cfunc_par3(self, i, j, k)
      class(cfunc_type), intent(inout) :: self
      integer, intent(in) :: i, j, k

      self%i = i
      self%j = j
      self%k = k

   end subroutine set_cfunc_par3

   !===================================================
   ! Setting functions mapping
   !===================================================
   ! Set optional coordinate mapping function.
   !================================================================
   subroutine set_func_mapping(self, mapping)
      class(func_type), intent(inout)        :: self
      type(mapping_type), target, intent(in) :: mapping

      self%mapping => mapping

   end subroutine

   !================================================================
   ! Remove coordinate mapping function.
   !================================================================
   subroutine unset_func_mapping(self)
      class(func_type), intent(inout) :: self

      self%mapping => null()

   end subroutine

   !===================================================
   ! Query function association state
   !===================================================
   function is_func(self)
      class(func_type), intent(in) :: self
      logical :: is_func
      
      is_func = self%is_associated 
      
   end function

   !===================================================
   ! The functions that actually pass around
   !===================================================
   function func(self, x) result(y)
      class(func_type), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp)             :: y

      real(dp) :: x1

      if (.not. self%is_associated) then
         if (associated(self%err)) then
            call self%err%push_err(ERR_INPUT, 'Trying to evaluate unassociated function')
         endif
         y = 0.0_dp
         return
      endif

      if (associated(self%mapping)) then 
           x1 = self%mapping%yx(x) 
        else
           x1 = x
      endif
     
      select case(self%n_par)
         case (1)
            y = self%p1(x1)
         case (2) 
            y = self%p2(x1, self%i)
         case (3)
            y = self%p3(x1, self%i, self%j)
         case (4)
            y = self%p4(x1, self%i, self%j, self%k)
         case default
            if (associated(self%err)) then
               call self%err%push_err(ERR_INTERNAL, 'No function is associated with func')
            endif
            y = 0.0_dp
            return
      end select
   end function func

   !================================================================
   ! Evaluate real derivative subroutine and return vector df/dp.
   !================================================================
   subroutine func_vec(self, x, yv)
      class(func_type), intent(in) :: self
      real(dp), intent(in)         :: x
      real(dp), intent(out)        :: yv(:)

      real(dp) :: x1

      if (.not. self%is_associated) then
         if (associated(self%err)) then
            call self%err%push_err(ERR_INPUT, 'Trying to evaluate unassociated function')
         endif
         yv = 0.0_dp
         return
      endif

      if (associated(self%mapping)) then 
           x1 = self%mapping%yx(x) 
        else
           x1 = x
      endif

      if (size(yv) < self%n_vec) then
         if (associated(self%err)) call self%err%push_err(ERR_INPUT, 'Derivative vector is too small')
         yv = 0.0_dp
         return
      endif

      if (associated(self%v1)) then
         call self%v1(x1, yv(1:self%n_vec))
      elseif (associated(self%v2)) then
         call self%v2(x1, self%i, yv(1:self%n_vec))
      elseif (associated(self%v3)) then
         call self%v3(x1, self%i, self%j, yv(1:self%n_vec))
      else
         if (associated(self%err)) call self%err%push_err(ERR_INTERNAL, 'No vector function is associated with func')
         yv = 0.0_dp
      endif
   end subroutine func_vec

   !================================================================
   ! Evaluate real function at complex x via its real part.
   !================================================================
   function func_c(self, x) result(y)
      class(func_type), intent(in) :: self
      complex(dp), intent(in) :: x
      real(dp)                :: y

      y = self%f(real(x, dp))

   end function func_c

   !================================================================
   ! Query complex-function association state.
   !================================================================
   function is_cfunc(self)
      class(cfunc_type), intent(in) :: self
      logical :: is_cfunc
      
      is_cfunc = self%is_associated 
      
   end function

   !================================================================
   ! Evaluate associated complex function.
   !================================================================
   function cfunc(self, x) result(y)
      class(cfunc_type), intent(in) :: self
      complex(dp), intent(in) :: x
      complex(dp)             :: y

      if (.not. self%is_associated) then
         if (associated(self%err)) then
            call self%err%push_err(ERR_INPUT, 'Trying to evaluate unassociated function')
         endif
         y = cmplx(0.0_dp, 0.0_dp, dp)
         return
      endif

      select case(self%n_par)
         case (1)
            y = self%p1(x)
         case (2) 
            y = self%p2(x, self%i)
         case (3)
            y = self%p3(x, self%i, self%j)
         case (4)
            y = self%p4(x, self%i, self%j, self%k)
         case default
            if (associated(self%err)) then
               call self%err%push_err(ERR_INTERNAL, 'No function is associated with cfunc')
            endif
            y = cmplx(0.0_dp, 0.0_dp, dp)
            return
      end select

   end function cfunc

   !================================================================
   ! Evaluate complex derivative subroutine and return vector df/dp.
   !================================================================
   subroutine cfunc_vec(self, x, yv)
      class(cfunc_type), intent(in) :: self
      complex(dp), intent(in)       :: x
      complex(dp), intent(out)      :: yv(:)

      if (.not. self%is_associated) then
         if (associated(self%err)) then
            call self%err%push_err(ERR_INPUT, 'Trying to evaluate unassociated function')
         endif
         yv = cmplx(0.0_dp, 0.0_dp, dp)
         return
      endif

      if (size(yv) < self%n_vec) then
         if (associated(self%err)) call self%err%push_err(ERR_INPUT, 'Derivative vector is too small')
         yv = cmplx(0.0_dp, 0.0_dp, dp)
         return
      endif

      if (associated(self%v1)) then
         call self%v1(x, yv(1:self%n_vec))
      elseif (associated(self%v2)) then
         call self%v2(x, self%i, yv(1:self%n_vec))
      elseif (associated(self%v3)) then
         call self%v3(x, self%i, self%j, yv(1:self%n_vec))
      else
         if (associated(self%err)) call self%err%push_err(ERR_INTERNAL, 'No vector function is associated with cfunc')
         yv = cmplx(0.0_dp, 0.0_dp, dp)
      endif
   end subroutine cfunc_vec

   !================================================================
   ! Evaluate complex function at real x (cast to complex).
   !================================================================
   function cfunc_r(self, x) result(y)
      class(cfunc_type), intent(in) :: self
      real(dp),    intent(in) :: x
      complex(dp)             :: y

      y = self%f(cmplx(x, 0.0_dp, dp))

   end function cfunc_r

   !================================================================
   ! Build a spectral grid of the given kind: fill nodes, quadrature
   ! weights, barycentric weights and the active-node range, by
   ! dispatching to the matching get_Gauss_*_grid kernel.
   !================================================================
   subroutine grid_build(self, N, kind, err)
      class(grid_type), intent(out)   :: self
      integer,          intent(in)    :: N, kind
      class(err_type),  intent(inout) :: err

      integer  :: istat, ept
      real(dp) :: xwork(N), qwork(N), lambda(N)   ! Radau work arrays

      if (err%has_error()) return

      self%kind = kind

      select case (kind)
      case (GRID_LOBATTO)
         allocate(self%x(0:N), self%qw(0:N), self%bary(0:N), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in grid_build (Lobatto)', err)) return
         call get_Gauss_Lobatto_grid(N, self%x, self%qw, self%bary, err)
         self%lo = 0
         self%hi = N

      case (GRID_LEGENDRE)
         allocate(self%x(0:N), self%qw(0:N), self%bary(0:N), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in grid_build (Legendre)', err)) return
         self%x(0)    = -1.0_dp;  self%x(N)    =  1.0_dp   ! boundary nodes (not collocation)
         self%qw(0)   =  0.0_dp;  self%qw(N)   =  0.0_dp
         self%bary(0) =  0.0_dp;  self%bary(N) =  0.0_dp
         call get_Gauss_Legendre_grid(N - 1, self%x(1:N-1), self%qw(1:N-1), self%bary(1:N-1), err)
         self%lo = 1
         self%hi = N - 1

      case (GRID_RADAU_M1, GRID_RADAU_P1)
         allocate(self%x(0:N), self%qw(0:N), self%bary(0:N), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in grid_build (Radau)', err)) return
         self%bary = 0.0_dp                          ! ghost entry stays 0
         ept = merge(-1, 1, kind == GRID_RADAU_M1)
         call get_Gauss_Radau_grid(N, ept, xwork, qwork, lambda, err)
         if (err%has_error()) return
         if (kind == GRID_RADAU_M1) then             ! endpoint -1 included; ghost at N
            self%x(0:N-1)    = xwork;  self%x(N)  =  1.0_dp
            self%qw(0:N-1)   = qwork;  self%qw(N) =  0.0_dp
            self%bary(0:N-1) = lambda
            self%lo = 0;  self%hi = N - 1
         else                                        ! endpoint +1 included; ghost at 0
            self%x(1:N)    = xwork;  self%x(0)  = -1.0_dp
            self%qw(1:N)   = qwork;  self%qw(0) =  0.0_dp
            self%bary(1:N) = lambda
            self%lo = 1;  self%hi = N
         endif

      case (GRID_PLAIN)
         !=========================================================
         ! Plain mesh: only allocate the node/weight storage. The nodes,
         ! weights and active range are filled in by the caller afterwards
         ! (e.g. copied from already-scaled spectral grids), so nothing is
         ! initialized here and the routine returns early.
         !=========================================================
         allocate(self%x(0:N), self%qw(0:N), stat=istat)
         if (istat /= 0) &
            call err%push_err(ERR_ALLOC, 'Allocation failed in grid_build (plain)')
         return

      case default
         call err%push_err(ERR_INTERNAL, 'grid_build: grid kind not implemented yet')
         return
      end select

      ! Spectral grids are constructed on the reference interval [-1,1].
      self%xmin   = -1.0_dp
      self%xmax   =  1.0_dp
      self%length =  2.0_dp
   end subroutine grid_build

   !================================================================
   ! Affine-map the nodes from the grid's current interval [self%xmin,
   ! self%xmax] onto [xmin, xmax], and rescale the quadrature weights by
   ! the length ratio. Repeatable: a later call maps from wherever the
   ! grid currently sits, so set_scale(2,3) then set_scale(-5,7) works.
   ! Barycentric weights are scale-invariant and stay untouched.
   !================================================================
   subroutine grid_set_scale(self, xmin, xmax, err)
      class(grid_type), intent(inout) :: self
      real(dp),         intent(in)    :: xmin, xmax
      class(err_type),  intent(inout) :: err

      real(dp) :: s

      if (err%has_error()) return

      s = (xmax - xmin) / self%length        ! length ratio (new / current)

      self%x      = xmin + (self%x - self%xmin) * s
      self%qw     = self%qw * s
      self%xmin   = xmin
      self%xmax   = xmax
      self%length = xmax - xmin
   end subroutine grid_set_scale

   !================================================================
   ! Build the differentiation operator on an existing grid: attach to
   ! the grid and fill D1, D2 by the barycentric + NST formulas. The
   ! matrices span only the active node range grid%lo:grid%hi, so they
   ! carry no formal-boundary (ghost) rows or columns.
   !================================================================
   subroutine diffmat_build(self, grid, err)
      class(diffmat_type),     intent(out)   :: self
      type(grid_type), target, intent(in)    :: grid
      class(err_type),         intent(inout) :: err

      integer  :: i, j, lo, hi, istat
      real(dp) :: t, t1, t2

      if (err%has_error()) return

      self%grid => grid
      lo = grid%lo
      hi = grid%hi

      allocate(self%D1(lo:hi, lo:hi), self%D2(lo:hi, lo:hi), stat=istat)
      if (alloc_failed(istat, 'Allocation failed in diffmat_build', err)) return

      associate (xv => grid%x, bary => grid%bary, D1 => self%D1, D2 => self%D2)

      ! First derivative off-diagonal entries (barycentric).
      do i = lo, hi
         do j = lo, i-1
            t = xv(i) - xv(j)
            D1(i,j) = bary(j)/(bary(i)*t)
         enddo
         do j = i+1, hi
            t = xv(i) - xv(j)
            D1(i,j) = bary(j)/(bary(i)*t)
         enddo
      enddo

      ! NST for D1 diagonal.
      do i = lo, hi
         t1 = 0.0_dp
         do j = lo, i-1
            t1 = t1 + D1(i,j)
         enddo
         do j = hi, i+1, -1
            t1 = t1 + D1(i,j)
         enddo
         D1(i,i) = -t1
      enddo

      ! Second derivative off-diagonal entries.
      do i = lo, hi
         do j = lo, i-1
            t = xv(i) - xv(j)
            D2(i,j) = 2*D1(i,j)*(D1(i,i) - 1/t)
         enddo
         do j = i+1, hi
            t = xv(i) - xv(j)
            D2(i,j) = 2*D1(i,j)*(D1(i,i) - 1/t)
         enddo
      enddo

      ! NST for D2 diagonal.
      do i = lo, hi
         t2 = 0.0_dp
         do j = lo, i-1
            t2 = t2 + D2(i,j)
         enddo
         do j = hi, i+1, -1
            t2 = t2 + D2(i,j)
         enddo
         D2(i,i) = -t2
      enddo

      end associate
   end subroutine diffmat_build

   !================================================================
   ! Build the symmetric-form differentiation operator on a Gauss-Lobatto
   ! grid: exploit the node symmetry (x_i = -x_{n-i}) for the off-diagonals
   ! and a bary-weighted NST for the diagonals. Requires a Lobatto grid.
   !================================================================
   subroutine diffmat_build_sym(self, grid, err)
      class(diffmat_type),     intent(out)   :: self
      type(grid_type), target, intent(in)    :: grid
      class(err_type),         intent(inout) :: err

      integer  :: i, j, n, istat
      real(dp) :: t, t1, t2

      if (err%has_error()) return

      if (grid%kind /= GRID_LOBATTO) then
         call err%push_err(ERR_INTERNAL, &
            & 'diffmat%build_sym requires a Gauss-Lobatto grid')
         return
      endif

      self%grid => grid
      n = ubound(grid%x, 1)

      allocate(self%D1(0:n, 0:n), self%D2(0:n, 0:n), stat=istat)
      if (alloc_failed(istat, 'Allocation failed in diffmat_build_sym', err)) return

      associate (xv => grid%x, bary => grid%bary, D1 => self%D1, D2 => self%D2)

      ! Off-diagonals from one half, reflected via x_i = -x_{n-i}.
      do i = 1, (n-1)/2
         do j = i+1, n-i
            t = xv(i) - xv(j)
            D1(i,j) =  1/t
            D1(j,i) = -D1(i,j)
            D2(i,j) = -2/(t*t)
            D2(j,i) =  D2(i,j)
         enddo
      enddo

      D1(0,0) = -dble(n)*(n+1)/4.0_dp
      D1(n,n) = -D1(0,0)

      do j = 1, n-1
         t = xv(0) - xv(j)
         D1(0,j) =  1/t
         D1(j,0) = -D1(0,j)
         D2(0,j) =  2*D1(0,j)*(D1(0,0) - 1.0_dp/t)
         D2(j,0) = -2/(t*t)
      enddo

      t = xv(0) - xv(n)

      D1(0,n) =  1/t
      D1(n,0) = -D1(0,n)
      D2(0,n) =  2*D1(0,n)*(D1(0,0) - 1.0_dp/t)
      D2(n,0) =  2*D1(n,0)*(D1(n,n) + 1.0_dp/t)

      do i = 0, n-1
         do j = i+1, n-i-1
            D1(n-i,n-j) = -D1(i,j)
            D1(n-j,n-i) = -D1(j,i)
            D2(n-i,n-j) =  D2(i,j)
            D2(n-j,n-i) =  D2(j,i)
         enddo
      enddo

      ! NST diagonals, symmetry-weighted by bary (= 1/P_N).
      do i = 0, n
         t1 = 0.0_dp
         t2 = 0.0_dp
         do j = 0, i-1
            t1 = t1 + D1(i,j)*bary(j)
            t2 = t2 + D2(i,j)*bary(j)
         enddo
         do j = n, i+1, -1
            t1 = t1 + D1(i,j)*bary(j)
            t2 = t2 + D2(i,j)*bary(j)
         enddo
         D1(i,i) = -t1/bary(i)
         D2(i,i) = -t2/bary(i)
      enddo

      end associate
   end subroutine diffmat_build_sym

   !================================================================
   ! Transform this operator's D1/D2 into its own Hermite-interpolation
   ! matrices DH1/DH2 (shape 0:n x 0:n+1, with the extra boundary-
   ! coupling column), stored on the operator itself.
   !
   ! Grid-agnostic: works for both the Gauss-Lobatto grid (bc_left =
   ! BC_DIRICHLET in the ECS split) and the Gauss-Radau endpoint = +1
   ! grid (bc_left = BC_REGULAR). Everything is indexed over the active
   ! range grid%lo:grid%hi; the free endpoint (subinterval junction) is
   ! the last active node grid%hi, and column grid%hi+1 is the extra
   ! boundary-coupling column. No formal-boundary (ghost) rows/columns.
   !================================================================
   subroutine diffmat_build_hermite(self, err)
      class(diffmat_type), intent(inout) :: self
      class(err_type),     intent(inout) :: err

      integer  :: i, j, lo, hi, istat
      real(dp) :: t1, t2

      if (err%has_error()) return

      if (.not. associated(self%grid)) then
         call err%push_err(ERR_INTERNAL, 'diffmat%build_hermite called on an unbuilt operator')
         return
      endif

      lo = self%grid%lo
      hi = self%grid%hi

      allocate(self%DH1(lo:hi, lo:hi+1), self%DH2(lo:hi, lo:hi+1), stat=istat)
      if (alloc_failed(istat, 'Allocation failed in diffmat%build_hermite', err)) return

      associate (D1m => self%D1, D2m => self%D2, xv => self%grid%x, &
                 DmH1 => self%DH1, DmH2 => self%DH2)
         do i = lo, hi
            do j = lo, hi-1
               DmH1(i,j) =  D1m(i,j)*(xv(i)-xv(hi))/(xv(j)-xv(hi))
               DmH2(i,j) = (D2m(i,j)*(xv(i)-xv(hi)) + 2*D1m(i,j))/(xv(j)-xv(hi))
            enddo
         enddo

         do i = lo, hi-1
            DmH1(i,hi)   = D1m(i,hi)*(1-(xv(i)-xv(hi))*D1m(hi,hi))
            DmH2(i,hi)   = D2m(i,hi)*(1-(xv(i)-xv(hi))*D1m(hi,hi)) - 2*D1m(i,hi)*D1m(hi,hi)
            DmH1(i,hi+1) = D1m(i,hi)*(xv(i)-xv(hi))
            DmH2(i,hi+1) = D2m(i,hi)*(xv(i)-xv(hi)) + 2*D1m(i,hi)
         enddo

         DmH1(hi,hi+1) = 1.0_dp
         DmH2(hi,hi+1) = 2*D1m(hi,hi)

         !===================================================
         ! Negative sum trick (NST) for the diagonal entries.
         ! See R. Baltensperger, M. R. Trummer, Spectral Differencing
         ! with a Twist, SIAM J. Sci. Comput. 24 (5) (2003) 1465-1487.
         !===================================================
         do i = lo, hi
            t1 = 0.0_dp
            t2 = 0.0_dp
            do j = lo, i-1
               t1 = t1 + DmH1(i,j)
               t2 = t2 + DmH2(i,j)
            enddo
            do j = hi, i+1, -1
               t1 = t1 + DmH1(i,j)
               t2 = t2 + DmH2(i,j)
            enddo
            DmH1(i,i) = -t1
            DmH2(i,i) = -t2
         enddo
      end associate
   end subroutine diffmat_build_hermite

   !================================================================
   ! Move both the attached grid and the differentiation matrices onto
   ! [xmin, xmax]. The grid is rescaled via grid%set_scale; the operators,
   ! built on the reference interval, pick up the chain-rule factor
   ! ds = length/(xmax-xmin): D1 by ds, D2 by ds**2.
   !
   ! The Hermite matrices need care: the boundary-coupling column hi+1
   ! carries one power of length less than the rest, because it multiplies
   ! a boundary value rather than a derivative. By dimensions,
   !   DH1(:,lo:hi) ~ 1/L  -> ds,   DH1(:,hi+1) dimensionless -> invariant;
   !   DH2(:,lo:hi) ~ 1/L^2 -> ds^2, DH2(:,hi+1) ~ 1/L         -> ds.
   !================================================================
   subroutine diffmat_set_scale(self, xmin, xmax, err)
      class(diffmat_type), intent(inout) :: self
      real(dp),            intent(in)    :: xmin, xmax
      class(err_type),     intent(inout) :: err

      real(dp) :: ds

      if (err%has_error()) return

      if (.not. associated(self%grid)) then
         call err%push_err(ERR_INTERNAL, 'diffmat%set_scale called on an unbuilt operator')
         return
      endif

      ds = self%grid%length / (xmax - xmin)   ! = 1/s, captured before the grid is rescaled

      call self%grid%set_scale(xmin, xmax, err)
      if (err%has_error()) return

      self%D1 = self%D1 * ds
      self%D2 = self%D2 * ds**2

      associate (lo => self%grid%lo, hi => self%grid%hi)
         if (allocated(self%DH1)) then
            self%DH1(:, lo:hi) = self%DH1(:, lo:hi) * ds       ! coupling column hi+1 is invariant
         endif
         if (allocated(self%DH2)) then
            self%DH2(:, lo:hi) = self%DH2(:, lo:hi) * ds**2
            self%DH2(:, hi+1)  = self%DH2(:, hi+1)  * ds
         endif
      end associate
   end subroutine diffmat_set_scale

   !================================================================
   ! LAPACK wrapper for real symmetric eigenproblems (DSYEVR).
   !================================================================
   subroutine lapack_solve_real_sym(A, n, lda, E, imin, imax, el, eu, vec, err)
      integer,  intent(in)               :: n, lda
      real(dp), intent(inout)            :: A(lda,n)
      real(dp), intent(out), allocatable :: E(:)
      integer,  intent(in),  optional    :: imin, imax
      real(dp), intent(in),  optional    :: el, eu
      real(dp), intent(out), allocatable, optional :: vec(:,:)
      class(err_type),     intent(inout) :: err

      integer  :: nev, ldz, lwork, liwork
      integer  :: l_il, l_iu, info
      real(dp) :: l_el, l_eu
      integer,  allocatable :: isuppz(:), iwork(:)
      real(dp), allocatable :: work(:), l_E(:)
      real(dp), allocatable :: z(:,:)
      integer   :: iworksize(1)
      real(dp)   :: worksize(1)
      real(dp)   :: abstol
      character :: range, jobz, uplo
      character(len=128) :: msg
      integer :: istat

      if (lda < n) then
         call err%push_err(ERR_INPUT, 'lda must be >= n')
         return
      endif

      if ((present(imin) .or. present(imax)) .and.  & 
      &   (present(el)   .or. present(eu))) then 
         call err%push_err(ERR_INPUT, 'Cannot decide which range of eigenvalues to use')
         return
      endif

      if (present(imin)) then
         l_il = imin
         if (present(imax)) then
            l_iu = imax
         else
            l_iu = n
         endif  
      elseif (present(imax)) then
         l_iu = imax
         l_il = 1
      elseif (present(el)) then
         l_el = el
         if (present(eu)) then
            l_eu = eu         
         else
            l_eu = huge(1.0_dp)
         endif
      elseif (present(eu)) then
         l_eu = eu
         l_el = -huge(1.0_dp)
      endif

      if (present(imin) .or. present(imax)) then 
         range = 'I'
         if (l_il < 1 .or. l_il > n) then
            call err%push_err(ERR_INPUT, 'imin is out of range')
            return
         endif
         if (l_iu < 1 .or. l_iu > n) then
            call err%push_err(ERR_INPUT, 'imax is out of range')
            return
         endif
         if (l_iu < l_il) then
            call err%push_err(ERR_INPUT, 'imax must be >= imin')
            return
         endif
         l_el = -huge(1.0_dp)
         l_eu =  huge(1.0_dp)
      elseif (present(el) .or. present(eu)) then
         range = 'V'
         if (l_eu <= l_el) then
            call err%push_err(ERR_INPUT, 'eu must be > el')
            return
         endif
         l_il = 1
         l_iu = n
      else
         range = 'A'
         l_il = 1
         l_iu = n
         l_el = -huge(1.0_dp)
         l_eu =  huge(1.0_dp)
      endif

   if(present(vec)) then
      jobz = 'V'
      if (range == 'I') then
         allocate(z(n, l_iu - l_il + 1), stat=istat)
      else 
         allocate(z(n, n), stat=istat)
      endif 
      if (alloc_failed(istat, 'Allocation failed in lapack_solve_real_sym', err)) return
      ldz = n
   else
      jobz = 'N'
      allocate(z(1,1), stat=istat)
      if (alloc_failed(istat, 'Allocation failed in lapack_solve_real_sym', err)) return
      ldz = 1
   endif

   uplo = 'U'

   allocate(l_E(n), isuppz(2*n), stat=istat)
   if (alloc_failed(istat, 'Allocation failed in lapack_solve_real_sym', err)) return

   abstol = dlamch('S')
   lwork = -1; liwork = -1

   call dsyevr(jobz, range, uplo, n, A, lda, l_el, l_eu, l_il, l_iu, abstol, nev, & 
   &    l_E, z, ldz, isuppz, worksize, lwork, iworksize, liwork, info)

   if (info /= 0) then
      write(msg, '(A,I0)') 'dsyevr query failed, info = ', info
      call err%push_err(ERR_LAPACK, msg)
      return
   endif

   lwork  = max(1, int(worksize(1), kind=kind(lwork)))
   liwork = max(1, iworksize(1))
   allocate(work(lwork), iwork(liwork), stat=istat)
   if (alloc_failed(istat, 'Allocation failed in lapack_solve_real_sym', err)) return

   call dsyevr(jobz, range, uplo, n, A, lda, l_el, l_eu, l_il, l_iu, abstol, nev, & 
   &    l_E, z, ldz, isuppz, work, lwork, iwork, liwork, info)

   if (info /= 0) then
      write(msg, '(A,I0)') 'dsyevr failed, info = ', info
      call err%push_err(ERR_LAPACK, msg)
      return
   endif
   allocate(E(nev), stat=istat)
   if (alloc_failed(istat, 'Allocation failed in lapack_solve_real_sym', err)) return

   if (present(vec)) then
      allocate(vec(n,nev), stat=istat)
      if (alloc_failed(istat, 'Allocation failed in lapack_solve_real_sym', err)) return
      vec(1:n, 1:nev) = z(1:n, 1:nev)
   endif

   E(1:nev) = l_E(1:nev)

   end subroutine

   !================================================================
   ! LAPACK wrapper for complex Hermitian eigenproblems (ZHEEVR).
   !================================================================
   subroutine lapack_solve_cmplx_herm(A, n, lda, E, imin, imax, el, eu, vec, err)
      integer,     intent(in)    :: n, lda
      complex(dp), intent(inout) :: A(lda,n)
      real(dp),    intent(out), allocatable :: E(:)
      integer,    intent(in),  optional :: imin, imax
      real(dp),    intent(in),  optional :: el, eu
      complex(dp), intent(out), allocatable, optional :: vec(:,:)
      class(err_type), intent(inout) :: err

      ! for zheevr
      integer  :: ldz, lwork, lrwork, liwork
      integer  :: l_il, l_iu, info, nev
      real(dp)  :: l_el, l_eu
      integer,    allocatable :: isuppz(:), iwork(:)
      complex(dp), allocatable :: work(:), z(:,:)
      real(dp),    allocatable :: rwork(:), l_E(:)
      integer    :: iworksize(1)
      complex(dp) :: worksize(1)
      real(dp)    :: rworksize(1)
      real(dp)    :: abstol
      character  :: range, jobz, uplo
      character(len=128) :: msg
      integer :: istat

      if (lda < n) then
         call err%push_err(ERR_INPUT, 'lda must be >= n')
         return
      endif

      if ((present(imin) .or. present(imax)) .and.  &
      &   (present(el)   .or. present(eu))) then
         call err%push_err(ERR_INPUT, 'Cannot decide which range of eigenvalues to use')
         return
      endif

      if (present(imin)) then
         l_il = imin
         if (present(imax)) then
            l_iu = imax
         else
            l_iu = n
         endif
      elseif (present(imax)) then
         l_iu = imax
         l_il = 1
      elseif (present(el)) then
         l_el = el
         if (present(eu)) then
            l_eu = eu
         else
            l_eu = huge(1.0_dp)
         endif
      elseif (present(eu)) then
         l_eu = eu
         l_el = -huge(1.0_dp)
      endif

      if (present(imin) .or. present(imax)) then
         range = 'I'
         if (l_il < 1 .or. l_il > n) then
            call err%push_err(ERR_INPUT, 'imin is out of range')
            return
         endif
         if (l_iu < 1 .or. l_iu > n) then
            call err%push_err(ERR_INPUT, 'imax is out of range')
            return
         endif
         if (l_iu < l_il) then
            call err%push_err(ERR_INPUT, 'imax must be >= imin')
            return
         endif
         l_el = -huge(1.0_dp)
         l_eu =  huge(1.0_dp)
      elseif (present(el) .or. present(eu)) then
         range = 'V'
         if (l_eu <= l_el) then
            call err%push_err(ERR_INPUT, 'eu must be > el')
            return
         endif
         l_il = 1
         l_iu = n
      else
         range = 'A'
         l_il = 1
         l_iu = n
         l_el = -huge(1.0_dp)
         l_eu =  huge(1.0_dp)
      endif

   if(present(vec)) then
      jobz = 'V'
      if (range == 'I') then
         allocate(z(n, l_iu - l_il + 1), stat=istat)
      else
         allocate(z(n, n), stat=istat)
      endif
      if (alloc_failed(istat, 'Allocation failed in lapack_solve_cmplx_herm', err)) return
      ldz = n
   else
      jobz = 'N'
      allocate(z(1,1), stat=istat)
      if (alloc_failed(istat, 'Allocation failed in lapack_solve_cmplx_herm', err)) return
      ldz = 1
   endif

   uplo = 'U'

   allocate(l_E(n), isuppz(2*n), stat=istat)
   if (alloc_failed(istat, 'Allocation failed in lapack_solve_cmplx_herm', err)) return

   abstol = dlamch('S')
   lwork  = -1
   lrwork = -1
   liwork = -1

   call zheevr(jobz, range, uplo, n, A, lda, l_el, l_eu, l_il, l_iu, &
   & abstol, nev, l_E, z, ldz, isuppz, worksize, lwork, &
   & rworksize, lrwork, iworksize, liwork, info)

   if (info /= 0) then
      write(msg, '(A,I0)') 'zheevr query failed, info = ', info
      call err%push_err(ERR_LAPACK, msg)
      return
   endif

   lwork  = max(1, int(dble(worksize(1)), kind=kind(lwork)))
   lrwork = max(1, int(rworksize(1), kind=kind(lrwork)))
   liwork = max(1, iworksize(1))
   allocate(work(lwork), rwork(lrwork), iwork(liwork), stat=istat)
   if (alloc_failed(istat, 'Allocation failed in lapack_solve_cmplx_herm', err)) return

   call zheevr(jobz, range, uplo, n, A, lda, l_el, l_eu, l_il, l_iu, abstol, nev, &
   &    l_E, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)

   if (info /= 0) then
      write(msg, '(A,I0)') 'zheevr failed, info = ', info
      call err%push_err(ERR_LAPACK, msg)
      return
   endif
   allocate(E(nev), stat=istat)
   if (alloc_failed(istat, 'Allocation failed in lapack_solve_cmplx_herm', err)) return
   E(1:nev) = l_E(1:nev)

   if (present(vec)) then
      allocate(vec(n, nev), stat=istat)
      if (alloc_failed(istat, 'Allocation failed in lapack_solve_cmplx_herm', err)) return
      vec(1:n, 1:nev) = z(1:n, 1:nev)
   endif

   end subroutine

   !================================================================
   ! LAPACK wrapper for complex non-Hermitian eigenproblems (ZGEEV).
   !================================================================
   subroutine lapack_solve_cmplx_nonherm(A, n, lda, Ec, vecl, vecr, err)
      integer,     intent(in)    :: n, lda
      complex(dp), intent(inout) :: A(lda,n)
      complex(dp), intent(out), allocatable :: Ec(:)
      complex(dp), intent(out), allocatable, optional :: vecl(:,:)
      complex(dp), intent(out), allocatable, optional :: vecr(:,:)
      class(err_type), intent(inout) :: err

      ! for zgeev
      integer  :: ldvl, ldvr, lwork, info
      complex(dp), allocatable :: work(:)
      real(dp),    allocatable :: rwork(:)
      complex(dp) :: worksize(1), dummy_vl(1,1), dummy_vr(1,1)
      character  :: jobvl, jobvr
      character(len=128) :: msg
      integer :: istat

      if (lda < n) then
         call err%push_err(ERR_INPUT, 'lda must be >= n')
         return
      endif

      if(present(vecl)) then
         jobvl = 'V'
         allocate(vecl(n,n), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in lapack_solve_cmplx_nonherm', err)) return
         ldvl = n
      else
         jobvl = 'N'
         ldvl = 1
      endif

      if(present(vecr)) then
         jobvr = 'V'
         allocate(vecr(n,n), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in lapack_solve_cmplx_nonherm', err)) return
         ldvr = n
      else
         jobvr = 'N'
         ldvr = 1
      endif
      allocate(rwork(2*n), Ec(n), stat=istat)
      if (alloc_failed(istat, 'Allocation failed in lapack_solve_cmplx_nonherm', err)) return

      lwork  = -1

      if (present(vecl) .and. present(vecr)) then
         call zgeev(jobvl, jobvr, n, A, lda, Ec, vecl, ldvl, vecr, ldvr, &
         &          worksize, lwork, rwork, info)
      elseif (present(vecl)) then
         call zgeev(jobvl, jobvr, n, A, lda, Ec, vecl, ldvl, dummy_vr, ldvr, &
         &          worksize, lwork, rwork, info)
      elseif (present(vecr)) then
         call zgeev(jobvl, jobvr, n, A, lda, Ec, dummy_vl, ldvl, vecr, ldvr, &
         &          worksize, lwork, rwork, info)
      else
         call zgeev(jobvl, jobvr, n, A, lda, Ec, dummy_vl, ldvl, dummy_vr, ldvr, &
         &          worksize, lwork, rwork, info)
      endif

      if (info /= 0) then
         write(msg, '(A,I0)') 'zgeev query failed, info = ', info
         call err%push_err(ERR_LAPACK, msg)
         return
      endif

      lwork = max(1, int(dble(worksize(1)), kind=kind(lwork)))
      allocate(work(lwork), stat=istat)
      if (alloc_failed(istat, 'Allocation failed in lapack_solve_cmplx_nonherm', err)) return

      if (present(vecl) .and. present(vecr)) then
         call zgeev(jobvl, jobvr, n, A, lda, Ec, vecl, ldvl, vecr, ldvr, &
         &          work, lwork, rwork, info)
      elseif (present(vecl)) then
         call zgeev(jobvl, jobvr, n, A, lda, Ec, vecl, ldvl, dummy_vr, ldvr, &
         &          work, lwork, rwork, info)
      elseif (present(vecr)) then
         call zgeev(jobvl, jobvr, n, A, lda, Ec, dummy_vl, ldvl, vecr, ldvr, &
         &          work, lwork, rwork, info)
      else
         call zgeev(jobvl, jobvr, n, A, lda, Ec, dummy_vl, ldvl, dummy_vr, ldvr, &
         &          work, lwork, rwork, info)
      endif

      if (info /= 0) then
         write(msg, '(A,I0)') 'zgeev failed, info = ', info
         call err%push_err(ERR_LAPACK, msg)
         return
      endif

   end subroutine

   !================================================================
   ! LAPACK wrapper for real non-symmetric eigenproblems (DGEEV).
   !================================================================
   subroutine lapack_solve_real_nonsym(A, n, lda, Ec, vecl, vecr, err)
      integer,     intent(in)    :: n, lda
      real(dp),    intent(inout) :: A(lda,n)
      complex(dp), intent(out), allocatable :: Ec(:)
      complex(dp), intent(out), allocatable, optional :: vecl(:,:)
      complex(dp), intent(out), allocatable, optional :: vecr(:,:)
      class(err_type), intent(inout) :: err

      ! for dgeev
      integer :: ldvl, ldvr, lwork, info, j
      real(dp), allocatable :: work(:)
      real(dp) :: worksize(1)

      real(dp), allocatable :: wr(:), wi(:)
      real(dp), allocatable :: vl(:,:), vr(:,:)

      character  :: jobvl, jobvr
      character(len=128) :: msg
      integer :: istat

      if (lda < n) then
         call err%push_err(ERR_INPUT, 'lda must be >= n')
         return
      endif
      allocate(Ec(n), wr(n), wi(n), stat=istat)
      if (alloc_failed(istat, 'Allocation failed in lapack_solve_real_nonsym', err)) return

      if(present(vecl)) then
         jobvl = 'V'
         allocate(vl(n,n), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in lapack_solve_real_nonsym', err)) return
         ldvl = n
      else
         jobvl = 'N'
         ldvl = 1
         allocate(vl(1,1), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in lapack_solve_real_nonsym', err)) return
      endif

      if(present(vecr)) then
         jobvr = 'V'
         allocate(vr(n,n), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in lapack_solve_real_nonsym', err)) return
         ldvr = n
      else
         jobvr = 'N'
         ldvr = 1
         allocate(vr(1,1), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in lapack_solve_real_nonsym', err)) return
      endif

      lwork  = -1

      call dgeev(jobvl, jobvr, n, A, lda, wr, wi, vl, ldvl, vr, ldvr, &
      &          worksize, lwork, info)

      if (info /= 0) then
         write(msg, '(A,I0)') 'dgeev query failed, info = ', info
         call err%push_err(ERR_LAPACK, msg)
         return
      endif

      lwork = max(1, int(worksize(1), kind=kind(lwork)))
      allocate(work(lwork), stat=istat)
      if (alloc_failed(istat, 'Allocation failed in lapack_solve_real_nonsym', err)) return

      call dgeev(jobvl, jobvr, n, A, lda, wr, wi, vl, ldvl, vr, ldvr, &
      &          work, lwork, info)

      if (info /= 0) then
         write(msg, '(A,I0)') 'dgeev failed, info = ', info
         call err%push_err(ERR_LAPACK, msg)
         return
      endif

      Ec(:) = cmplx(wr(:), wi(:), kind=dp)

      if (present(vecl)) then
         allocate(vecl(n,n), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in lapack_solve_real_nonsym', err)) return
         j = 1
         do while (j <= n)
            if (wi(j) > 0.0_dp) then
               if (j == n) then
                  call err%push_err(ERR_LAPACK, 'Complex eigenvalue pair truncated by DGEEV')
                  return
               endif
               vecl(1:n, j)   = cmplx(vl(1:n, j),  vl(1:n, j+1), kind=dp)
               vecl(1:n, j+1) = cmplx(vl(1:n, j), -vl(1:n, j+1), kind=dp)
               j = j + 2
            else if (wi(j) < 0.0_dp) then
               call err%push_err(ERR_LAPACK, 'Unexpected DGEEV conjugate-pair ordering')
               return
            else
               vecl(1:n, j) = cmplx(vl(1:n, j), 0.0_dp, kind=dp)
               j = j + 1
            endif
         enddo
      endif

      if (present(vecr)) then
         allocate(vecr(n,n), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in lapack_solve_real_nonsym', err)) return
         j = 1
         do while (j <= n)
            if (wi(j) > 0.0_dp) then
               if (j == n) then
                  call err%push_err(ERR_LAPACK, 'Complex eigenvalue pair truncated by DGEEV')
                  return
               endif
               vecr(1:n, j)   = cmplx(vr(1:n, j),  vr(1:n, j+1), kind=dp)
               vecr(1:n, j+1) = cmplx(vr(1:n, j), -vr(1:n, j+1), kind=dp)
               j = j + 2
            else if (wi(j) < 0.0_dp) then
               call err%push_err(ERR_LAPACK, 'Unexpected DGEEV conjugate-pair ordering')
               return
            else
               vecr(1:n, j) = cmplx(vr(1:n, j), 0.0_dp, kind=dp)
               j = j + 1
            endif
         enddo
      endif

   end subroutine

   !================================================================
   ! ARPACK shift-invert solver for real symmetric matrices.
   !================================================================
   subroutine arpack_solve_real_sym(A, n, lda, E0, E, nev, vec, err)
      integer,         intent(in)    :: n, lda
      real(dp),        intent(inout) :: A(lda,n)
      real(dp),        intent(in)    :: E0
      integer,         intent(in)    :: nev
      real(dp),        intent(out), allocatable :: E(:)
      real(dp),        intent(out), allocatable, optional :: vec(:,:)
      class(err_type), intent(inout) :: err

      logical, allocatable :: select(:)
      real(dp), dimension(:), allocatable :: d, resid, workd, workl
      real(dp), allocatable :: v(:,:)

      character(1) :: bmat
      character(2) :: which

      integer :: iparam(11), ipntr(14), l_nev
      integer :: ncv, info, lworkl, nconv, maxitr
      integer :: ishfts, mode, ido, ierr, ldv, j
      real(dp) :: tol
      real(dp) :: sigma
      logical :: rvec
      ! paramerers for lapack linear solver
      integer :: lwork
      real(dp) :: twork(1)
      integer, allocatable :: ipiv(:)
      real(dp), allocatable :: work(:)
      character(len=128) :: msg
      integer :: istat

      if (lda < n) then
         call err%push_err(ERR_INPUT, 'lda must be >= n')
         return
      endif

      if (nev < 1) then
         call err%push_err(ERR_INPUT, 'nE must be positive')
         return
      endif

      if (nev >= n) then
         call err%push_err(ERR_INPUT, 'nE must be less than N')
         return
      endif

      ncv = min(n, max(2*nev + 2, 20)) 

      lworkl = ncv*(ncv + 8)
      ldv    = n

      allocate(resid(n), &
      &  v(n, ncv),      &
      &  d(ncv),         &
      &  workd(3*n),     &
      &  workl(lworkl),  &
      &  ipiv(n),        &
      &  select(ncv), stat=istat)
      if (alloc_failed(istat, 'Allocation failed in arpack_solve_real_sym', err)) return

      sigma = E0

      do j = 1,n
         A(j,j) = A(j,j) - sigma
      enddo

      lwork = -1
      call dsytrf('U', n, A, lda, ipiv, twork, lwork, info)
      if (info /= 0) then
         write(msg, '(A,I0)') 'dsytrf workspace query failed, info = ', info
         call err%push_err(ERR_LAPACK, msg)
         return
      endif
      lwork = max(1, int(twork(1), kind=kind(lwork)))

      allocate(work(lwork), stat=istat)
      if (alloc_failed(istat, 'Allocation failed in arpack_solve_real_sym', err)) return
      call dsytrf('U', n, A, lda, ipiv, work, lwork, info)
      if (info /= 0) then
         write(msg, '(A,I0)') 'dsytrf factorization failed, info = ', info
         call err%push_err(ERR_LAPACK, msg)
         return
      endif

      bmat  = 'I'
      which = 'LM' 

      tol    = -1.0_dp
      ido    = 0
      info   = 0

      ishfts = 1
      maxitr = 300
      mode   = 3

      iparam = 0
      ipntr  = 0

      iparam(1) = ishfts 
      iparam(3) = maxitr 
      iparam(7) = mode 

      l_nev = nev

      do
         call dsaupd(ido, bmat, n, which, l_nev, tol, resid,   &
         &  ncv, v, ldv, iparam, ipntr, workd, workl, &
         &  lworkl, info)

         if (ido .eq. -1 .or. ido .eq. 1 ) then

            call dcopy (n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
            call dsytrs('U', n, 1, A, lda, ipiv, workd(ipntr(2)), n, ierr)
            if (ierr.ne. 0) then 
               write(msg, '(A,I0)') 'Error with linear solver dsytrs in arpack diagonalization, info = ', ierr
               call err%push_err(ERR_LAPACK, msg)
               return
            endif
         else
            exit
         endif
      enddo

      if (info == 1) then
         call err%push_warning('dsaupd: maximum number of iterations reached')
      elseif (info /= 0) then
         write(msg, '(A,I0)') 'Error with dsaupd, info = ', info
         call err%push_err(ERR_ARPACK, msg)
         return
      endif

      if (present(vec)) then 
         rvec = .true.
      else 
         rvec = .false.
      endif

      call dseupd (rvec, 'A', select, d, v, ldv, sigma,    &
      & bmat, n, which, l_nev, tol, resid, ncv, &
      & v, ldv, iparam, ipntr, workd, workl,    & 
      & lworkl, ierr)

      if (ierr .ne. 0) then
         write(msg, '(A,I0)') 'Error with dseupd, info = ', ierr
         call err%push_err(ERR_ARPACK, msg)
         nconv = 0
         allocate(E(0), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in arpack_solve_real_sym', err)) return

         if (present(vec)) then
            allocate(vec(n, 0), stat=istat)
            if (alloc_failed(istat, 'Allocation failed in arpack_solve_real_sym', err)) return
         endif
         return
      else
         nconv = iparam(5) 
         allocate(E(nconv), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in arpack_solve_real_sym', err)) return

         E(1:nconv) = d(1:nconv)

         if (present(vec)) then
            allocate(vec(n, nconv), stat=istat)
            if (alloc_failed(istat, 'Allocation failed in arpack_solve_real_sym', err)) return

            vec(1:n, 1:nconv) = v(1:n, 1:nconv)
         endif
      endif

   end subroutine

   !================================================================
   ! ARPACK shift-invert solver for real non-symmetric matrices.
   !================================================================
   subroutine arpack_solve_real_nonsym(A, n, lda, E0, Ec, nev, vec, err)
      integer,         intent(in)    :: n, lda
      real(dp),        intent(inout) :: A(lda,n)
      real(dp),        intent(in)    :: E0
      integer,         intent(in)    :: nev
      complex(dp),     intent(out), allocatable :: Ec(:)
      complex(dp),     intent(out), allocatable, optional :: vec(:,:)
      class(err_type), intent(inout) :: err

      logical, allocatable :: select(:)
      real(dp), dimension(:), allocatable :: resid, workd, workev, &
      & workl, dr, di
      real(dp), dimension(:,:), allocatable :: v

      character(1) :: bmat
      character(2) :: which

      integer :: iparam(11), ipntr(14), l_nev
      integer :: ncv, info, lworkl, nconv, maxitr
      integer :: ishfts, mode, ido, ierr, ldv, j
      real(dp) :: tol
      real(dp) :: sigmar, sigmai
      logical :: rvec
      ! paramerers for lapack linear solver
      integer, allocatable :: ipiv(:)
      character(len=128) :: msg
      integer :: istat
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (lda < n) then
         call err%push_err(ERR_INPUT, 'lda must be >= n')
         return
      endif

      if (nev < 1) then
         call err%push_err(ERR_INPUT, 'nE must be positive')
         return
      endif

      if (nev >= n-1) then
         call err%push_err(ERR_INPUT, 'nE must be less than N-1')
         return
      endif

      ncv = min(n, max(2*nev + 2, 20)) 

      lworkl = 3*ncv**2 + 6*ncv
      ldv    = n

      allocate(resid(n), &
      &  v(ldv, ncv),    &
      &  dr(ncv),        &
      &  di(ncv),        &
      &  workd(3*n),     &
      &  workev(3*ncv),  &
      &  workl(lworkl),  &
      &  ipiv(n),        &
      &  select(ncv), stat=istat)
      if (alloc_failed(istat, 'Allocation failed in arpack_solve_real_nonsym', err)) return

      sigmar = E0
      sigmai = 0.0_dp

      do j = 1,n
         A(j,j) = A(j,j) - sigmar
      enddo

      call dgetrf(n, n, A, lda, ipiv, info)
      if (info /= 0) then
         write(msg, '(A,I0)') 'dgetrf factorization failed, info = ', info
         call err%push_err(ERR_LAPACK, msg)
         return
      endif

      bmat  = 'I'
      which = 'LM' 

      tol    = -1.0_dp
      ido    = 0
      info   = 0

      ishfts = 1
      maxitr = 300
      mode   = 3

      iparam = 0
      ipntr  = 0

      iparam(1) = ishfts 
      iparam(3) = maxitr 
      iparam(7) = mode 

      l_nev = nev

      do
         call dnaupd(ido, bmat, n, which, l_nev, tol, resid, ncv, v,   & 
         &   ldv, iparam, ipntr, workd, workl, lworkl, info) 

         if (ido .eq. -1 .or. ido .eq. 1) then

            call dcopy (n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
            call dgetrs('N', n, 1, A, lda, ipiv, workd(ipntr(2)), n, ierr)
            if (ierr.ne. 0) then 
               write(msg, '(A,I0)') 'Error with linear solver dgetrs in arpack diagonalization, info = ', ierr
               call err%push_err(ERR_LAPACK, msg)
               return
            endif
         else
            exit
         endif
      enddo

      if (info == 1) then
         call err%push_warning('dnaupd: maximum number of iterations reached')
      elseif (info /= 0) then
         write(msg, '(A,I0)') 'Error with dnaupd, info = ', info
         call err%push_err(ERR_ARPACK, msg)
         return
      endif

      if (present(vec)) then 
         rvec = .true.
      else 
         rvec = .false.
      endif

      call dneupd(rvec, 'A', select, dr, di, v, ldv, sigmar, sigmai, &
      & workev, bmat, n, which, l_nev, tol, resid, ncv, v,  &
      & ldv, iparam, ipntr, workd, workl, lworkl, ierr)

      if (ierr .ne. 0) then
         write(msg, '(A,I0)') 'Error with dneupd, info = ', ierr
         call err%push_err(ERR_ARPACK, msg)
         nconv = 0
         allocate(Ec(0), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in arpack_solve_real_nonsym', err)) return

         if (present(vec)) then
            allocate(vec(n, 0), stat=istat)
            if (alloc_failed(istat, 'Allocation failed in arpack_solve_real_nonsym', err)) return
         endif
         return
      else
         nconv = iparam(5)
         allocate(Ec(nconv), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in arpack_solve_real_nonsym', err)) return
         Ec(1:nconv) = cmplx(dr(1:nconv), di(1:nconv), kind=dp)

         if (present(vec)) then
            allocate(vec(n, nconv), stat=istat)
            if (alloc_failed(istat, 'Allocation failed in arpack_solve_real_nonsym', err)) return

            j = 1
            do while (j <= nconv)

               if (di(j) > 0.0_dp) then
                  ! Complex conjugate pair: (j, j+1)
                  if (j == nconv) then
                     call err%push_err(ERR_ARPACK, 'Complex pair truncated')
                     return
                  endif

                  ! eigenvector for lambda = dr(j) + i*di(j)
                  vec(1:n, j)   = cmplx(v(1:n, j),  v(1:n, j+1), kind=dp)

                  ! eigenvector for conj(lambda) = dr(j) - i*di(j)
                  vec(1:n, j+1) = cmplx(v(1:n, j), -v(1:n, j+1), kind=dp)

                  j = j + 2
               else if (di(j) < 0.0_dp) then
                  call err%push_err(ERR_ARPACK, 'Unexpected ARPACK conjugate-pair ordering')
                  return
               else
                  ! Real eigenvalue => real eigenvector
                  vec(1:n, j) = cmplx(v(1:n, j), 0.0_dp, kind=dp)
                  j = j + 1
               endif
            enddo
         endif
      endif

   end subroutine

   !================================================================
   ! ARPACK shift-invert solver for complex non-Hermitian matrices.
   !================================================================
   subroutine arpack_solve_cmplx_nonherm(A, n, lda, Ec0, Ec, nev, vec, err)
      integer,     intent(in)    :: n, lda
      complex(dp), intent(inout) :: A(lda,n)
      complex(dp), intent(in)    :: Ec0
      integer,    intent(in)    :: nev
      complex(dp), intent(out), allocatable :: Ec(:)
      complex(dp), intent(out), allocatable, optional :: vec(:,:)
      class(err_type), intent(inout) :: err

      integer, allocatable :: ipiv(:)
      logical, allocatable :: select(:)
      real(dp),    dimension(:),   allocatable :: rwork
      complex(dp), dimension(:),   allocatable :: d, resid, workd, &
      & workev, workl
      complex(dp), dimension(:,:), allocatable :: v

      character(1) :: bmat
      character(2) :: which

      integer    :: iparam(11), ipntr(14), l_nev
      integer    :: ncv, info, lworkl, nconv, maxitr
      integer    :: ishfts, mode, ido, ierr, ldv, j
      real(dp)    :: tol
      complex(dp) :: sigma
      logical    :: rvec
      character(len=128) :: msg
      integer :: istat

      if (lda < n) then
         call err%push_err(ERR_INPUT, 'lda must be >= n')
         return
      endif

      if (nev < 1) then
         call err%push_err(ERR_INPUT, 'nE must be positive')
         return
      endif

      if (nev >= n-1) then
         call err%push_err(ERR_INPUT, 'nE must be less than N - 1')
         return
      endif

      ncv = min(n, max(2*nev + 2, 20))

      lworkl = 3*ncv**2 + 5*ncv
      ldv    = n

      allocate(resid(n), &
      &  v(n, ncv),      &
      &  d(ncv),         &
      &  workd(3*n),     &
      &  workev(2*ncv),  &
      &  workl(lworkl),  &
      &  rwork(ncv),     &
      &  ipiv(n),        &
      &  select(ncv), stat=istat)
      if (alloc_failed(istat, 'Allocation failed in arpack_solve_cmplx_nonherm', err)) return

      sigma = Ec0

      do j = 1,n
         A(j,j) = A(j,j) - sigma
      enddo

      call zgetrf(n, n, A, lda, ipiv, info)
      if (info /= 0) then
         write(msg, '(A,I0)') 'zgetrf failed, info = ', info
         call err%push_err(ERR_LAPACK, msg)
         return
      endif

      bmat  = 'I'
      which = 'LM'

      tol    = -1.0_dp
      ido    = 0
      info   = 0

      ishfts = 1
      maxitr = 300
      mode   = 3

      iparam = 0
      ipntr  = 0

      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode

      l_nev = nev

      do
         call znaupd(ido, bmat, n, which, l_nev, tol, resid, ncv, &
         &  v, ldv, iparam, ipntr, workd, workl,         &
         &  lworkl, rwork, info)

         if (ido .eq. -1 .or. ido .eq. 1) then

            call zcopy(n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
            call zgetrs('N', n, 1, A, lda, ipiv, workd(ipntr(2)), n, ierr)

            if (ierr .ne. 0 ) then
               write(msg, '(A,I0)') 'Error with linear solver in arpack diagonalization, info = ', ierr
               call err%push_err(ERR_LAPACK, msg)
               return
            end if
         else
            exit
         endif
      enddo

      if (info == 1) then
         call err%push_warning('znaupd: maximum number of iterations reached')
      elseif (info /= 0) then
         write(msg, '(A,I0)') 'Error with znaupd, info = ', info
         call err%push_err(ERR_ARPACK, msg)
         return
      endif

      if (present(vec)) then
         rvec = .true.
      else
         rvec = .false.
      endif

      call zneupd (rvec, 'A', select, d, v, ldv, sigma, &
      &     workev, bmat, n, which, l_nev, tol,  &
      &     resid, ncv, v, ldv, iparam, ipntr,   &
      workd, workl, lworkl, rwork, ierr)

      if (ierr .ne. 0) then
         write(msg, '(A,I0)') 'Error with zneupd, info = ', ierr
         call err%push_err(ERR_ARPACK, msg)
         nconv = 0
         allocate(Ec(0), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in arpack_solve_cmplx_nonherm', err)) return

         if (present(vec)) then
            allocate(vec(n, 0), stat=istat)
            if (alloc_failed(istat, 'Allocation failed in arpack_solve_cmplx_nonherm', err)) return
         endif
         return
      else
         nconv = iparam(5)
         allocate(Ec(nconv), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in arpack_solve_cmplx_nonherm', err)) return
         Ec(1:nconv) = d(1:nconv)

         if (present(vec)) then
            allocate(vec(n, nconv), stat=istat)
            if (alloc_failed(istat, 'Allocation failed in arpack_solve_cmplx_nonherm', err)) return
            vec(1:n, 1:nconv) = v(1:n, 1:nconv)
         endif
      endif

   end subroutine

   !================================================================
   ! ARPACK shift-invert solver for complex Hermitian matrices.
   !================================================================
   subroutine arpack_solve_cmplx_herm(A, n, lda, E0, E, nev, vec, err)
      integer,         intent(in)    :: n, lda
      complex(dp),     intent(inout) :: A(lda,n)
      real(dp),        intent(in)    :: E0
      integer,         intent(in)    :: nev
      real(dp),        intent(out), allocatable :: E(:)
      complex(dp),     intent(out), allocatable, optional :: vec(:,:)
      class(err_type), intent(inout) :: err

      integer, allocatable :: ipiv(:)
      complex(dp)          :: twork(1)
      integer              :: lwork

      logical, allocatable :: select(:)
      real(dp),    dimension(:),   allocatable :: rwork
      complex(dp), dimension(:),   allocatable :: d, resid, workd, &
                                                & workev, workl, work
      complex(dp), dimension(:,:), allocatable :: v

      character(1) :: bmat
      character(2) :: which

      integer     :: iparam(11), ipntr(14), l_nev
      integer     :: ncv, info, lworkl, nconv, maxitr
      integer     :: ishfts, mode, ido, ierr, ldv, j
      real(dp)    :: tol
      complex(dp) :: sigma
      logical    :: rvec
      character(len=128) :: msg
      integer :: istat

      if (lda < n) then
         call err%push_err(ERR_INPUT, 'lda must be >= n')
         return
      endif

      if (nev < 1) then
         call err%push_err(ERR_INPUT, 'nE must be positive')
         return
      endif

      if (nev >= n-1) then
         call err%push_err(ERR_INPUT, 'nE must be less than N - 1')
         return
      endif

      ncv = min(n, max(2*nev + 2, 20))

      lworkl = 3*ncv**2 + 5*ncv
      ldv    = n

      allocate(resid(n), &
      &  v(n, ncv),      &
      &  d(ncv),         &
      &  workd(3*n),     &
      &  workev(2*ncv),  &
      &  workl(lworkl),  &
      &  rwork(ncv),     &
      &  ipiv(n),        &
      &  select(ncv), stat=istat)
      if (alloc_failed(istat, 'Allocation failed in arpack_solve_cmplx_herm', err)) return

      sigma = cmplx(E0, 0.0_dp, kind=dp)

      do j = 1,n
         A(j,j) = A(j,j) - sigma
      enddo

      lwork = -1
      call zhetrf('U', n, A, lda, ipiv, twork, lwork, info)
      if (info /= 0) then
         write(msg, '(A,I0)') 'zhetrf workspace query failed, info = ', info
         call err%push_err(ERR_LAPACK, msg)
         return
      endif

      lwork = max(1, int(dble(twork(1)), kind=kind(lwork)))
      allocate(work(lwork), stat=istat)
      if (alloc_failed(istat, 'Allocation failed in arpack_solve_cmplx_herm', err)) return

      call zhetrf('U', n, A, lda, ipiv, work, lwork, info)
      if (info /= 0) then
         write(msg, '(A,I0)') 'zhetrf failed, info = ', info
         call err%push_err(ERR_LAPACK, msg)
         return
      endif

      bmat  = 'I'
      which = 'LM'

      tol    = -1.0_dp
      ido    = 0
      info   = 0

      ishfts = 1
      maxitr = 300
      mode   = 3

      iparam = 0
      ipntr  = 0

      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode

      l_nev = nev

      do
         call znaupd(ido, bmat, n, which, l_nev, tol, resid, ncv, &
         &  v, ldv, iparam, ipntr, workd, workl,         &
         &  lworkl, rwork, info)

         if (ido .eq. -1 .or. ido .eq. 1) then

            call zcopy (n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
            call zhetrs('U', n, 1, A, lda, ipiv, workd(ipntr(2)), n, ierr)

            if (ierr .ne. 0 ) then
               write(msg, '(A,I0)') 'Error with linear solver in arpack diagonalization, info = ', ierr
               call err%push_err(ERR_LAPACK, msg)
               return
            end if
         else
            exit
         endif
      enddo

      if (info == 1) then
         call err%push_warning('znaupd: maximum number of iterations reached')
      elseif (info /= 0) then
         write(msg, '(A,I0)') 'Error with znaupd, info = ', info
         call err%push_err(ERR_ARPACK, msg)
         return
      endif

      if (present(vec)) then
         rvec = .true.
      else
         rvec = .false.
      endif

      call zneupd (rvec, 'A', select, d, v, ldv, sigma, &
      &     workev, bmat, n, which, l_nev, tol,  &
      &     resid, ncv, v, ldv, iparam, ipntr,   &
      workd, workl, lworkl, rwork, ierr)

      if (ierr .ne. 0) then
         write(msg, '(A,I0)') 'Error with zneupd, info = ', ierr
         call err%push_err(ERR_ARPACK, msg)
         nconv = 0
         allocate(E(0), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in arpack_solve_cmplx_herm', err)) return

         if (present(vec)) then
            allocate(vec(n, 0), stat=istat)
            if (alloc_failed(istat, 'Allocation failed in arpack_solve_cmplx_herm', err)) return
         endif
         return
      else
         nconv = iparam(5)
         allocate(E(nconv), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in arpack_solve_cmplx_herm', err)) return
         E(1:nconv) = dble(d(1:nconv))

         if (present(vec)) then
            allocate(vec(n, nconv), stat=istat)
            if (alloc_failed(istat, 'Allocation failed in arpack_solve_cmplx_herm', err)) return
            vec(1:n, 1:nconv) = v(1:n, 1:nconv)
         endif
      endif

   end subroutine
   
   !================================================================
   ! ARPACK shift-invert solver for complex symmetric non-Hermitian matrices.
   !================================================================
   subroutine arpack_solve_cmplx_symm_nonherm(A, n, lda, Ec0, Ec, nev, vec, err)
      integer,     intent(in)    :: n, lda
      complex(dp), intent(inout) :: A(lda,n)
      complex(dp), intent(in)    :: Ec0
      integer,     intent(in)    :: nev
      complex(dp), intent(out), allocatable :: Ec(:)
      complex(dp), intent(out), allocatable, optional :: vec(:,:)
      class(err_type), intent(inout) :: err
   
      integer,  allocatable :: ipiv(:)
      logical,  allocatable :: select(:)
      real(dp), allocatable :: rwork(:)
      complex(dp), allocatable :: d(:), resid(:), workd(:), workev(:), workl(:)
      complex(dp), allocatable :: v(:,:), fwork(:)
      complex(dp) :: work_query(1)
   
      character(1) :: bmat, uplo
      character(2) :: which
   
      integer :: iparam(11), ipntr(14), l_nev
      integer :: ncv, info, lworkl, nconv, maxitr
      integer :: ishfts, mode, ido, ierr, ldv, j
      integer :: lwork_fact
      integer :: istat
      real(dp) :: tol
      complex(dp) :: sigma
      logical :: rvec
      character(len=128) :: msg
   
      if (lda < n) then
         call err%push_err(ERR_INPUT, 'lda must be >= n')
         return
      endif
   
      if (nev < 1) then
         call err%push_err(ERR_INPUT, 'nev must be positive')
         return
      end if
   
      if (nev >= n-1) then
         call err%push_err(ERR_INPUT, 'nev must be less than N - 1')
         return
      end if
   
      ncv = min(n, max(2*nev + 2, 20))
      ldv = n
      lworkl = 3*ncv*ncv + 5*ncv
   
      allocate(resid(n),      &
               v(n,ncv),      &
               d(ncv),        &
               workd(3*n),    &
               workev(2*ncv), &
               workl(lworkl), &
               rwork(ncv),    &
               ipiv(n),       &
               select(ncv),   &
               stat=istat)
      if (alloc_failed(istat, 'Allocation failed in arpack_solve_cmplx_symm_nonherm', err)) return
   
      sigma = Ec0
   
      !---------------------------------------------------------------
      ! Shift: A <- A - sigma*I
      !---------------------------------------------------------------
      do j = 1, n
         A(j,j) = A(j,j) - sigma
      end do
   
      !---------------------------------------------------------------
      ! Complex symmetric factorization of (A - sigma I)
      !
      ! Use only one triangle of A. Here we choose the lower triangle.
      ! zsytrf/zsytrs are for complex symmetric (not Hermitian) matrices.
      !---------------------------------------------------------------
      uplo = 'L'
   
      ! Workspace query
      call zsytrf(uplo, n, A, lda, ipiv, work_query, -1, info)
      if (info /= 0) then
         write(msg, '(A,I0)') 'zsytrf workspace query failed, info = ', info
         call err%push_err(ERR_LAPACK, msg)
         return
      end if
   
      lwork_fact = max(1, nint(real(work_query(1), kind=dp)))
      allocate(fwork(lwork_fact), stat=istat)
      if (alloc_failed(istat, 'Allocation failed in arpack_solve_cmplx_symm_nonherm', err)) return
   
      call zsytrf(uplo, n, A, lda, ipiv, fwork, lwork_fact, info)
      if (info /= 0) then
         write(msg, '(A,I0)') 'zsytrf failed, info = ', info
         call err%push_err(ERR_LAPACK, msg)
         return
      end if
   
      !---------------------------------------------------------------
      ! ARPACK setup
      ! mode = 3  => shift-invert mode for general complex problems
      !---------------------------------------------------------------
      bmat  = 'I'
      which = 'LM'
   
      tol    = -1.0_dp
      ido    = 0
      info   = 0
   
      ishfts = 1
      maxitr = 300
      mode   = 3
   
      iparam = 0
      ipntr  = 0
   
      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode
   
      l_nev = nev
   
      do
         call znaupd(ido, bmat, n, which, l_nev, tol, resid, ncv, &
                     v, ldv, iparam, ipntr, workd, workl,         &
                     lworkl, rwork, info)
   
         if (ido == -1 .or. ido == 1) then
            ! y = OP*x = inv(A - sigma I) * x
            call zcopy(n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
   
            call zsytrs(uplo, n, 1, A, lda, ipiv, workd(ipntr(2)), n, ierr)
            if (ierr /= 0) then
               write(msg, '(A,I0)') 'zsytrs failed in shift-invert solve, info = ', ierr
               call err%push_err(ERR_LAPACK, msg)
               return
            end if
         else
            exit
         end if
      end do
   
      if (info == 1) then
         call err%push_warning('znaupd: maximum number of iterations reached')
      else if (info /= 0) then
         write(msg, '(A,I0)') 'Error with znaupd, info = ', info
         call err%push_err(ERR_ARPACK, msg)
         return
      end if
   
      rvec = present(vec)
   
      call zneupd(rvec, 'A', select, d, v, ldv, sigma,   &
                  workev, bmat, n, which, l_nev, tol,    &
                  resid, ncv, v, ldv, iparam, ipntr,     &
                  workd, workl, lworkl, rwork, ierr)
   
      if (ierr /= 0) then
         write(msg, '(A,I0)') 'Error with zneupd, info = ', ierr
         call err%push_err(ERR_ARPACK, msg)
   
         nconv = 0
         allocate(Ec(0), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in arpack_solve_cmplx_symm_nonherm', err)) return
   
         if (present(vec)) then
            allocate(vec(n,0), stat=istat)
            if (alloc_failed(istat, 'Allocation failed in arpack_solve_cmplx_symm_nonherm', err)) return
         end if
         return
      end if
   
      nconv = iparam(5)
   
      allocate(Ec(nconv), stat=istat)
      if (alloc_failed(istat, 'Allocation failed in arpack_solve_cmplx_symm_nonherm', err)) return
      Ec(1:nconv) = d(1:nconv)
   
      if (present(vec)) then
         allocate(vec(n,nconv), stat=istat)
         if (alloc_failed(istat, 'Allocation failed in arpack_solve_cmplx_symm_nonherm', err)) return
         vec(:,1:nconv) = v(:,1:nconv)
      end if
   
   end subroutine arpack_solve_cmplx_symm_nonherm

   !===================================================
   ! Extrapolate a real-axis function to complex points.
   ! The real sampling nodes are the Gauss-Legendre nodes,
   ! so neither endpoint is evaluated and the function is
   ! never sampled at an infinite endpoint, whatever the
   ! integration limits. Coefficient chopping uses the
   ! Legendre expansion, and the value is obtained by the
   ! ordinary barycentric interpolation formula.
   !===================================================
   subroutine extrapol_cmplx(xcv, fcv, a, b, func, err)
      implicit none
      complex(dp), intent(in)  :: xcv(:)
      complex(dp), intent(out) :: fcv(:)
      real(dp), intent(in)     :: a, b

      type(func_type) :: func
      class(err_type), intent(inout) :: err

      integer, parameter :: nmax = 500

      real(dp), allocatable :: yv(:), qv(:), wv(:), x0v(:)
      real(dp), allocatable :: fv(:), cv(:), acv(:)
      real(dp), allocatable :: Pm(:,:)

      complex(dp) :: cnum, cden
      integer     :: nx, i, j, k, n, istat

      !===================================================
      ! Variables to find cutoff of the Legendre series
      ! for details see
      ! L. Aurentz, L. N. Trefethen, Chopping a Chebyshev
      ! Series, ACM Trans. Math. Softw. 43 (4) (Jan. 2017).
      !===================================================
      integer  :: j2, j3, plateauPoint, cutoff
      real(dp) :: e1, e2, r, tol, tol76
      logical  :: plateau

      n = 50

      nx = size(xcv)
      if (size(fcv) /= nx) then
         call err%push_err(ERR_INPUT, 'fcv has wrong size')
         fcv = (0.0_dp, 0.0_dp)
         return
      endif

      tol = epsilon(1.0_dp)
      plateau = .false.

      allocate(yv(nmax), qv(nmax), wv(nmax), x0v(nmax), &
             & fv(nmax), cv(nmax), acv(0:nmax-1), stat=istat)
      if (istat /= 0) then
         call err%push_err(ERR_ALLOC, 'Allocation failed in extrapol_cmplx')
         fcv = (0.0_dp, 0.0_dp)
         return
      endif

      do while (.true.)

         if (n > nmax) then
            call err%push_err(ERR_INTERNAL, 'extrapol_cmplx: n exceeds nmax buffer size')
            fcv = (0.0_dp, 0.0_dp)
            return
         endif

         !===================================================
         ! Gauss-Legendre nodes, quadrature weights and
         ! barycentric weights on [-1,1], scaled to (a;b).
         !===================================================
         call get_Gauss_Legendre_grid(n, yv(1:n), qv(1:n), wv(1:n), err)
         if (err%has_error()) return

         x0v(1:n) = ((b - a)*yv(1:n) + a + b)/2

         do i = 1, n
            fv(i) = func%f(x0v(i))
         enddo

         if (plateau) exit

         !===================================================
         ! Legendre coefficients used only for chopping.
         !===================================================
         call create_legendre_matrix(yv(1:n), Pm, err)
         if (err%has_error()) return

         call compute_real_legendre_coefficients(fv(1:n), Pm, qv(1:n), 0, acv(0:n-1), err)
         if (err%has_error()) return

         cv(1:n) = abs(acv(0:n-1))

         call normalize_decay_envelope(cv(1:n))

         do j = 2, n
            j2 = nint(1.25_dp*j + 5.0_dp)
            if (j2 > n) exit ! there is no plateau

            e1 = cv(j); e2 = cv(j2)

            if (e1 <= 0.0_dp) then
               plateau = .true.
            else
               r = 3.0_dp*(1.0_dp - log(e1) / log(tol))
               plateau = (e2 / e1 > r)
            endif

            if (plateau) then
               plateauPoint = j - 1
               exit
            endif
         enddo

         ! Setting cutoff point
         if (plateau) then
            if (cv(plateauPoint) <= tiny(1.0_dp)) then
               cutoff = plateauPoint
            else

               tol76 = tol**(7.0_dp/6.0_dp)
               j3    = count(cv(1:n) >= tol76)

               if (j3 < j2) then
                  j2    = j3 + 1
                  cv(j2) = tol76
               end if

               cv(1:j2) = log10(max(cv(1:j2), tiny(1.0_dp)))

               do j = 1, j2
                  cv(j) = cv(j) - ((j-1)*log10(tol))/(3*(j2-1))
               end do

               k      = minloc(cv(1:j2), 1)
               cutoff = max(k - 1, 1)
            endif

            n = cutoff
         else
            if (n >= 300) then
               call err%push_err(ERR_INTERNAL, 'Approximation of V, B, G, or their derivatives did not converge &
          & with 300 collocation points. The result of extrapolation to the complex plane &
          & will be bad. Use complex valued functions instead.')
               fcv = (0.0_dp, 0.0_dp)
               return
            endif

            n = n + 30
         endif

      enddo

      !===================================================
      ! Ordinary barycentric interpolation at complex points xcv.
      ! Barycentric weights wv are supplied by the grid routine.
      !===================================================
      do j = 1, nx
         cnum = sum(wv(1:n)*fv(1:n)/(xcv(j) - x0v(1:n)))
         cden = sum(wv(1:n)/(xcv(j) - x0v(1:n)))
         fcv(j) = cnum/cden
      enddo

   end subroutine extrapol_cmplx

   !===================================================
   ! Extrapolate a real-axis vector function to complex
   ! points on the Gauss-Legendre nodes. Used for derivative
   ! callbacks that return all d/dp components in one call.
   !===================================================
   subroutine extrapol_cmplx_vec(xcv, fcv, a, b, func, err)
      implicit none
      complex(dp), intent(in)  :: xcv(:)
      complex(dp), intent(out) :: fcv(:,:)
      real(dp), intent(in)     :: a, b

      type(func_type) :: func
      class(err_type), intent(inout) :: err

      integer, parameter :: nmax = 500

      real(dp), allocatable :: yv(:), qv(:), wv(:), x0v(:)
      real(dp), allocatable :: fv(:,:)
      real(dp), allocatable :: acv(:), cvnorm(:), cwork(:), tmpv(:)
      real(dp), allocatable :: Pm(:,:)

      complex(dp) :: cnum, cden
      integer     :: nx, nv, i, j, k, n, p, istat

      !===================================================
      ! Variables to find cutoff of the Legendre series
      ! for details see
      ! L. Aurentz, L. N. Trefethen, Chopping a Chebyshev
      ! Series, ACM Trans. Math. Softw. 43 (4) (Jan. 2017).
      !===================================================
      integer  :: j2, j3, plateauPoint, cutoff
      real(dp) :: e1, e2, r, tol, tol76
      logical  :: plateau

      n = 50

      nx = size(xcv)
      nv = size(fcv, 2)
      if (size(fcv, 1) /= nx) then
         call err%push_err(ERR_INPUT, 'fcv has wrong size')
         fcv = (0.0_dp, 0.0_dp)
         return
      endif

      if (nv /= func%n_vec) then
         call err%push_err(ERR_INPUT, 'fcv has wrong derivative-vector size')
         fcv = (0.0_dp, 0.0_dp)
         return
      endif

      allocate(yv(nmax), qv(nmax), wv(nmax), x0v(nmax),           &
         &     fv(nv,nmax), acv(0:nmax-1), cvnorm(nmax),          &
         &     cwork(nmax), tmpv(nv), stat=istat)
      if (istat /= 0) then
         call err%push_err(ERR_ALLOC, 'Allocation failed in extrapol_cmplx_vec')
         fcv = (0.0_dp, 0.0_dp)
         return
      endif

      tol = epsilon(1.0_dp)
      plateau = .false.

      do while (.true.)

         if (n > nmax) then
            call err%push_err(ERR_INTERNAL, 'extrapol_cmplx_vec: n exceeds nmax buffer size')
            fcv = (0.0_dp, 0.0_dp)
            return
         endif

         !===================================================
         ! Gauss-Legendre nodes, quadrature weights and
         ! barycentric weights on [-1,1], scaled to (a;b).
         !===================================================
         call get_Gauss_Legendre_grid(n, yv(1:n), qv(1:n), wv(1:n), err)
         if (err%has_error()) return

         x0v(1:n) = ((b - a)*yv(1:n) + a + b)/2

         do i = 1, n
            call func%fv(x0v(i), tmpv)
            fv(1:nv, i) = tmpv(1:nv)
         enddo

         if (plateau) exit

         !===================================================
         ! Legendre coefficients used only for chopping.
         !===================================================
         call create_legendre_matrix(yv(1:n), Pm, err)
         if (err%has_error()) return

         cvnorm(1:n) = 0.0_dp
         do p = 1, nv
            call compute_real_legendre_coefficients(fv(p, 1:n), Pm, qv(1:n), 0, acv(0:n-1), err)
            if (err%has_error()) return
            cwork(1:n) = abs(acv(0:n-1))
            call normalize_decay_envelope(cwork(1:n))
            cvnorm(1:n) = max(cvnorm(1:n), cwork(1:n))
         enddo

         do j = 2, n
            j2 = nint(1.25_dp*j + 5.0_dp)
            if (j2 > n) exit ! there is no plateau

            e1 = cvnorm(j); e2 = cvnorm(j2)

            if (e1 <= 0.0_dp) then
               plateau = .true.
            else
               r = 3*(1 - log(e1) / log(tol))
               plateau = (e2 / e1 > r)
            endif

            if (plateau) then
               plateauPoint = j - 1
               exit
            endif
         enddo

         ! Setting cutoff point
         if (plateau) then
            if (cvnorm(plateauPoint) <= tiny(1.0_dp)) then
               cutoff = plateauPoint
            else

               tol76 = tol**(7.0_dp/6.0_dp)
               j3    = count(cvnorm(1:n) >= tol76)

               if (j3 < j2) then
                  j2    = j3 + 1
                  cvnorm(j2) = tol76
               end if

               cvnorm(1:j2) = log10(max(cvnorm(1:j2), tiny(1.0_dp)))

               do j = 1, j2
                  cvnorm(j) = cvnorm(j) - ((j-1)*log10(tol))/(3*(j2-1))
               end do

               k      = minloc(cvnorm(1:j2), 1)
               cutoff = max(k - 1, 1)
            endif

            n = cutoff
         else
            if (n >= 300) then
               call err%push_err(ERR_INTERNAL, 'Approximation of V, B, G, or their derivatives did not converge &
          & with 300 collocation points. The result of extrapolation to the complex plane &
          & will be bad. Use complex valued functions instead.')
               fcv = (0.0_dp, 0.0_dp)
               return
            endif

            n = n + 30
         endif

      enddo

      !===================================================
      ! Ordinary barycentric interpolation at complex points xcv.
      ! Barycentric weights wv are supplied by the grid routine.
      !===================================================
      do j = 1, nx
         cden = sum(wv(1:n)/(xcv(j) - x0v(1:n)))
         do p = 1, nv
            cnum = sum(wv(1:n)*fv(p, 1:n)/(xcv(j) - x0v(1:n)))
            fcv(j, p) = cnum/cden
         enddo
      enddo

   end subroutine extrapol_cmplx_vec

   !================================================================
   ! Subroutine for sorting real energies.
   ! sort_mode: 0 no sort, 1 ascending, 2 descending.
   ! If E0 is present, sorting key is |E - E0|.
   !================================================================
   subroutine sort_real_energies(Ev, Eindv, E0, sort_mode, err)
      implicit none
      real(dp),             intent(in)              :: Ev(:)
      integer, allocatable, intent(out)             :: Eindv(:)
      real(dp),             intent(in),    optional :: E0
      integer,              intent(in),    optional :: sort_mode
      class(err_type),      intent(inout), optional :: err

      integer :: nE, i, j, temp_ind, istat
      integer :: l_sort_mode
      real(dp), allocatable :: dev(:)
      real(dp) :: temp_dev

      nE = size(Ev)
      allocate(Eindv(nE), dev(nE), stat=istat)
      if (istat /= 0) then
         if (present(err)) call err%push_err(ERR_ALLOC, 'Allocation failed in sort_real_energies')
         return
      endif

      if (present(sort_mode)) then
         l_sort_mode = sort_mode
      else
         l_sort_mode = 1
      endif

      if (l_sort_mode < 0 .or. l_sort_mode > 2) then
         if (present(err)) call err%push_err(ERR_INPUT, 'sortE must be 0, 1, or 2')
         return
      endif

      !================================================================
      ! initialize indices + precompute sorting keys
      !================================================================
      do i = 1, nE
         Eindv(i) = i
         if (present(E0)) then
            dev(i) = abs(Ev(i) - E0)
         else
            dev(i) = Ev(i)
         endif
      end do

      if (l_sort_mode == 0) return

      !================================================================
      ! insertion sort indices by selected order of sorting key
      !================================================================
      do i = 2, nE
         temp_ind = Eindv(i)
         temp_dev = dev(temp_ind)

         j = i
         do while (j > 1)
            if (l_sort_mode == 1) then
               if (dev(Eindv(j-1)) <= temp_dev) exit
            else
               if (dev(Eindv(j-1)) >= temp_dev) exit
            endif
            Eindv(j) = Eindv(j-1)
            j = j - 1
         end do

         Eindv(j) = temp_ind
      end do
   end subroutine sort_real_energies

   !================================================================
   ! Sort complex energies by target distance, real energy, or width.
   ! sortE/sortW: 0 no sort, 1 ascending, 2 descending.
   !================================================================
   subroutine sort_complex_energies(Ecv, Eindv, Ec0, sortE, sortW, rtol, atol, err)
      complex(dp),          intent(in)              :: Ecv(:)
      integer, allocatable, intent(out)             :: Eindv(:)
      complex(dp),          intent(in),    optional :: Ec0
      integer,              intent(in),    optional :: sortE, sortW
      real(dp),             intent(in),    optional :: rtol, atol
      class(err_type),      intent(inout), optional :: err
   
      integer  :: nE, i, j, temp_ind, istat
      integer  :: l_sortE_mode, l_sortW_mode
      real(dp) :: rtol_, atol_, t1, t2
      real(dp), allocatable :: key1(:), key2(:)
      logical :: use_target, is_desc
   
      nE = size(Ecv)
      allocate(Eindv(nE), key1(nE), key2(nE), stat = istat)
      if (istat /= 0) then
         if (present(err)) call err%push_err(ERR_ALLOC, 'Allocation failed in sort_complex_energies')
         return
      endif
   
      if (present(rtol)) then
         rtol_ = rtol
      else
         rtol_ = 100.0_dp*epsilon(1.0_dp)
      endif

      if (present(atol)) then
         atol_ = atol
      else
         atol_ = 0.0_dp
      endif

      if (present(sortE)) then
         l_sortE_mode = sortE
      else
         l_sortE_mode = 0
      endif

      if (present(sortW)) then
         l_sortW_mode = sortW
      else
         l_sortW_mode = 0
      endif

      if (l_sortE_mode < 0 .or. l_sortE_mode > 2) then
         if (present(err)) call err%push_err(ERR_INPUT, 'sortE must be 0, 1, or 2')
         return
      endif

      if (l_sortW_mode < 0 .or. l_sortW_mode > 2) then
         if (present(err)) call err%push_err(ERR_INPUT, 'sortW must be 0, 1, or 2')
         return
      endif

      if (l_sortE_mode /= 0 .and. l_sortW_mode /= 0) then
         if (present(err)) call err%push_err(ERR_INPUT, 'sortE and sortW cannot both be nonzero')
         return
      endif

      use_target = present(Ec0) .and. (.not. present(sortE)) .and. (.not. present(sortW))

      do i = 1, nE
         Eindv(i) = i
      enddo

      if (.not. use_target .and. l_sortE_mode == 0 .and. l_sortW_mode == 0) then
         return
      endif

      is_desc = (l_sortE_mode == 2) .or. (l_sortW_mode == 2)
   
      ! Precompute sort keys
      do i = 1, nE
         if (use_target) then
            key1(i) = abs(Ecv(i) - Ec0)
            key2(i) = real(Ecv(i), kind = dp)
         elseif (l_sortW_mode /= 0) then
            key1(i) = -2.0_dp*aimag(Ecv(i))
            key2(i) = real(Ecv(i), kind = dp)
         else
            key1(i) = real(Ecv(i), kind = dp)
            key2(i) = -2.0_dp*aimag(Ecv(i))
         endif
      end do
   
      ! Insertion sort of indices by (key1, then key2)
      do i = 2, nE
         temp_ind = Eindv(i)
         t1 = key1(temp_ind)
         t2 = key2(temp_ind)
         j = i
         do while (j > 1)
            if (.not. comes_after(Eindv(j-1), t1, t2)) exit
            Eindv(j) = Eindv(j-1)
            j = j - 1
         end do
         Eindv(j) = temp_ind
      end do

   contains

      pure logical function approx_equal(a, b) result(eq)
         real(dp), intent(in) :: a, b
         real(dp) :: tol
         tol = max(atol_, rtol_ * max(abs(a), abs(b)))
         eq = (abs(a - b) <= tol)
      end function approx_equal

      pure logical function comes_after(idx, t1, t2) result(after)
         integer,  intent(in) :: idx
         real(dp), intent(in) :: t1, t2
         if (.not. approx_equal(key1(idx), t1)) then
            if (is_desc) then
               after = (key1(idx) < t1)
            else
               after = (key1(idx) > t1)
            endif
         else if (.not. approx_equal(key2(idx), t2)) then
            if (is_desc) then
               after = (key2(idx) < t2)
            else
               after = (key2(idx) > t2)
            endif
         else
            after = .false.  
         end if
      end function comes_after
   
   end subroutine sort_complex_energies
   
   !================================================================
   ! Inverse tangent function for complex argument
   !================================================================
   pure elemental function atan_cmplx(z)
      complex(dp), intent(in) :: z
      complex(dp)             :: atan_cmplx

      complex(dp) :: i_u, t0, t1

      i_u = (0.0_dp, 1.0_dp)
      t0  = i_u + z
      t1  = i_u - z
      atan_cmplx = 0.5_dp*i_u*log(t0/t1)

   end function atan_cmplx

   !================================================================
   ! Inverse hyperbolic tangent for real argument
   !================================================================
   pure elemental function atanh_real(x) 
      real(dp), intent(in) :: x
      real(dp)             :: atanh_real

      atanh_real = 0.5_dp*(log(1.0_dp + x) - log(1.0_dp - x))

   end function atanh_real

   !================================================================
   ! Inverse hyperbolic sine for real argument
   !================================================================
   pure elemental function asinh_real(x) 
      real(dp), intent(in) :: x
      real(dp)             :: asinh_real

      asinh_real = sign(1.0_dp, x)*log(abs(x) + sqrt(1.0_dp + x*x))

   end function asinh_real

   !================================================================
   ! Hyperbolic tangent function for complex argument
   !================================================================
   pure elemental function tanh_cmplx(z) 
      complex(dp), intent(in) :: z
      complex(dp)             :: tanh_cmplx

      complex(dp) :: tz

      tz = exp(2*z)
      tanh_cmplx = (tz - 1.0_dp)/(tz + 1.0_dp)

   end function tanh_cmplx

   !================================================================
   ! Hyperbolic sin function for complex argument
   !================================================================
   pure elemental function sinh_cmplx(z) result(sinh_z)
      implicit none
      complex(dp), intent(in) :: z
      complex(dp) :: sinh_z

      sinh_z = 0.5_dp*(exp(z) - exp(-z))

   end function sinh_cmplx

   !================================================================
   ! Mapping function xy(y).
   !================================================================
   pure elemental function xy(self, y)
      class(mapping_type), intent(in) :: self
      real(dp), intent(in) :: y

      real(dp) :: xy, t0

      associate(mn => self%mn, beta => self%beta, ye => self%ye)
         select case (mn)
            case (0)
               xy = y
            case (1)
               xy = 2/pi*atan(beta*(y - ye))
            case (2)
               t0 = (y/ye)**beta
               xy = (t0 - 1)/(t0 + 1)
            case (3)
               xy = tanh(beta*(y - ye))
            case (4)
               xy = tanh(sinh(beta*(y - ye)))
            case (5)
               xy = 1 - exp(-(y/ye)**beta)
         end select
      end associate

   end function

   !================================================================
   ! Mapping function x(y) for complex coordinate input.
   !================================================================
   pure elemental function xyc(self, y)
      class(mapping_type), intent(in) :: self
      complex(dp), intent(in) :: y

      complex(dp) :: xyc, t0

      associate(mn => self%mn, beta => self%beta, ye => self%ye)
         select case (mn)
            case (0)
               xyc = y
            case (1)
               xyc = 2/pi*atan_cmplx(beta*(y - ye))
            case (2)
               t0  = (y/ye)**beta
               xyc = (t0 - 1)/(t0 + 1)
            case (3)
               xyc = tanh_cmplx(beta*(y - ye))
            case (4)
               xyc = tanh_cmplx(sinh_cmplx(beta*(y - ye)))
            case (5)
               xyc = 1 - exp(-(y/ye)**beta)
         end select
      end associate

   end function

   !================================================================
   ! Upper mapped boundary x(+infinity) for active mapping.
   !================================================================
   function xyInf(self)
      class(mapping_type), intent(in) :: self
      real(dp) :: xyInf

      xyInf = 0.0_dp

      associate(mn => self%mn, beta => self%beta, ye => self%ye)
         select case (mn)
            case (0)
               if (associated(self%err)) then
                  call self%err%push_err(ERR_INPUT, 'Mapping must be specified if rmax = +infinity')
               endif
               return
            case (1, 3, 4)
               if (beta > 0.0_dp) then 
                  xyInf = 1.0_dp
               else        
                  if (associated(self%err)) then
                     call self%err%push_err(ERR_INPUT, 'rmax = +infinity is incompatible with current mapping values')
                  endif
                  return
               endif
            case (2, 5)
               if (beta > 0.0_dp .and. ye > 0.0_dp) then 
                  xyInf = 1.0_dp
               else        
                  if (associated(self%err)) then
                     call self%err%push_err(ERR_INPUT, 'rmax = +infinity is incompatible with current mapping values')
                  endif
                  return
               endif
         end select
      end associate

   end function

   !================================================================
   ! Lower mapped boundary x(-infinity) for active mapping.
   !================================================================
   function xyMInf(self)
      class(mapping_type), intent(in) :: self
      real(dp) :: xyMInf

      xyMInf = 0.0_dp

      associate(mn => self%mn, beta => self%beta, ye => self%ye)
         select case (mn)
            case (0)
               if (associated(self%err)) then
                  call self%err%push_err(ERR_INPUT, 'Mapping must be specified if rmin = -infinity')
               endif
               return
            case (1, 3, 4)
               if (beta > 0.0_dp) then 
                  xyMInf = -1.0_dp
               else        
                  if (associated(self%err)) then
                     call self%err%push_err(ERR_INPUT, 'rmin = -infinity is incompatible with current mapping values')
                  endif
                  return
               endif
            case (2)
               if (beta < 0.0_dp .and. ye < 0.0_dp) then 
                  xyMInf = -1.0_dp
               else        
                  if (associated(self%err)) then
                     call self%err%push_err(ERR_INPUT, 'rmin = -infinity is incompatible with current mapping values')
                  endif
                  return
               endif
            case (5)
               if (beta < 0.0_dp .and. ye < 0.0_dp) then 
                  xyMInf = 0.0_dp
               else        
                  if (associated(self%err)) then
                     call self%err%push_err(ERR_INPUT, 'rmin = -infinity is incompatible with current mapping values')
                  endif
                  return
               endif
         end select
      end associate

   end function

   !================================================================
   ! Inverse mapping y(x) for real coordinate input.
   !================================================================
   function yx(self, x)
      class(mapping_type), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: yx

      associate(mn => self%mn, beta => self%beta, ye => self%ye)
         select case (mn)
            case (0)
               yx = x
            case (1)
               yx = ye + tan(x*pi/2)/beta
            case (2)
               yx = ye*((1 + x)/(1 - x))**(1/beta)
            case (3)
               yx = ye + atanh_real(x)/beta
            case (4)
               yx = ye + asinh_real(atanh_real(x))/beta
            case (5)
               yx = ye*(-log(1 - x))**(1/beta)
         end select
      end associate

   end function

   !================================================================
   ! Derivative dy/dx of the selected mapping.
   !================================================================
   pure elemental function dydx(self, x) 
      class(mapping_type), intent(in) :: self
      real(dp), intent(in):: x
      real(dp) :: dydx

      associate(mn => self%mn, beta => self%beta, ye => self%ye)
         select case (mn)
            case (0)
               dydx = 1.0_dp
            case (1)
               dydx = pi/(2*beta*cos(pi/2*x)**2)
            case (2)
               dydx = 2*ye*((1 + x)/(1 - x))**(1/beta - 1)/ &
               &      (beta*(1 - x)**2)
            case (3)
               dydx = 1/(beta*(1 - x*x)) 
            case (4)
               dydx = 1/(beta*(1 - x*x)*sqrt(atanh_real(x)**2+1)) 
            case (5)
               dydx = ye*(-log(1 - x))**(1/beta - 1)/(beta*(1 - x))
         end select
      end associate

   end function

   !================================================================
   ! Derivative dx/dy, inverse of dydx where defined.
   !================================================================
   pure elemental function dxdyx(self, x) 
      class(mapping_type), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: dxdyx

      associate(mn => self%mn, beta => self%beta, ye => self%ye)
         select case (mn)
            case (0)
               dxdyx = 1.0_dp
            case (1)
               dxdyx = (2*beta*cos(0.5_dp*pi*x)**2)/pi
            case (2)
               dxdyx = (beta*(1 - x)**2)/                     & 
               &       (2*ye*((1 + x)/(1 - x))**(1/beta - 1))
            case (3)
               dxdyx = beta*(1.0_dp - x*x)
            case (4)
               dxdyx = beta*(1.0_dp - x*x)*sqrt(atanh_real(x)**2 + 1.0_dp)
            case (5)
               dxdyx = (beta*(1 - x))/(ye*(-log(1 - x))**(1/beta - 1))
         end select
      end associate

   end function dxdyx

   !================================================================
   ! Mapping-related scalar function Q(x) used in transformed operator.
   !================================================================
   pure elemental function Qx(self, x) 
      class(mapping_type), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: Qx
      real(dp) :: t0, t1

      associate(mn => self%mn, beta => self%beta, ye => self%ye)
         select case (mn)
            case (0)
               Qx = 0.0_dp
            case (1)
               Qx = tan(x*pi/2)*pi/2
            case (2)
               Qx = (beta*x + 1)/(beta*(1 - x*x))
            case (3)
               Qx = x/(1 - x*x)
            case (4)
               t0 = atanh_real(x)
               t1 = 2*t0**2 + 2
               Qx = (x*t1 - t0)/(t1*(1 - x**2))
            case (5)
               t0 = -log(1 - x)
               Qx = (1 + (1/beta - 1)/t0)/(2*(1 - x))
         end select
      end associate

   end function

   !================================================================
   ! Mapping-related scalar function F(x) used in transformed operator.
   !================================================================
   pure elemental function Fx(self, x)
      class(mapping_type), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: Fx
      real(dp) :: t0

      associate(mn => self%mn, beta => self%beta, ye => self%ye)
         select case (mn)
            case (0)
               Fx = 0.0_dp
            case (1)
               Fx = (pi/2)**2
            case (2)
               Fx = (1 - 1/(beta*beta))/(1 - x*x)**2
            case (3)
               Fx = 1/(x*x - 1)**2
            case (4)
               t0 = atanh_real(x)**2 + 1
               Fx = (4*t0**2 + t0 - 3)/(4*t0**2*(x*x - 1)**2)
            case (5)
               t0 = log(1 - x)
               Fx = (1 + (1 - 1/(beta*beta))/t0**2)/(4*(1 - x)**2)
         end select
      end associate

   end function

   !================================================================
   ! Compute Gauss-Legendre nodes, quadrature weights, and
   ! barycentric weights.
   !
   ! For N nodes, the quadrature is exact for polynomials of degree
   ! up to 2*N - 1. Nodes are computed by Newton iteration in theta,
   ! x = cos(theta), using symmetry of the roots.
   !================================================================
   subroutine get_Gauss_Legendre_grid(N, xv, qw, baryv, err)
      implicit none
      integer, intent(in) :: N
      real(dp), dimension(1:N), intent(out) :: xv, qw, baryv
      class(err_type), intent(inout) :: err

      integer :: i, k, iter, maxiter
      integer :: Nh, mid, psign
      real(dp), dimension(1:N) :: thetav, toldv
      real(dp), dimension(1:N) :: xhv, sinv
      real(dp), dimension(1:N) :: pkm2v, pkm1v, pkv
      real(dp) :: eps, maxdiff
      real(dp) :: t0, t1, t2

      if (N < 1) then
         call err%push_err(ERR_INPUT, 'Number of Gauss-Legendre nodes must be positive')
         xv = 0.0_dp
         qw = 0.0_dp
         baryv = 0.0_dp
         return
      endif

      if (N == 1) then
         xv(1) = 0.0_dp
         qw(1) = 2.0_dp
         baryv(1) = 1.0_dp
         return
      endif

      Nh = N/2
      if (mod(N - 1, 2) == 0) then
         psign = 1
      else
         psign = -1
      endif

      ! For large N, Newton corrections stall at roundoff of the
      ! Legendre recurrence; a stricter tolerance gives false failures.
      eps = 512*epsilon(1.0_dp)
      maxiter = 50

      do i = 1, Nh
         thetav(i) = (dble(i) - 0.25_dp)*pi/(dble(N) + 0.5_dp)
      enddo

      maxdiff = huge(1.0_dp)
      iter = 0
      do while (maxdiff > eps)
         iter = iter + 1
         if (iter > maxiter) then
            call err%push_err(ERR_INTERNAL, 'Gauss_Legendre: Newton iteration did not converge')
            xv = 0.0_dp
            qw = 0.0_dp
            baryv = 0.0_dp
            return
         endif

         toldv(1:Nh) = thetav(1:Nh)
         xhv(1:Nh) = cos(toldv(1:Nh))
         sinv(1:Nh) = sin(toldv(1:Nh))

         pkm2v(1:Nh) = 1.0_dp
         pkm1v(1:Nh) = xhv(1:Nh)
         do k = 2, N
            pkv(1:Nh)   = ((2*k - 1)*xhv(1:Nh)*pkm1v(1:Nh) - (k - 1)*pkm2v(1:Nh))/k
            pkm2v(1:Nh) = pkm1v(1:Nh)
            pkm1v(1:Nh) = pkv(1:Nh)
         enddo

         thetav(1:Nh) = toldv(1:Nh) + pkm1v(1:Nh)*sinv(1:Nh) &
                      & /(dble(N)*(pkm2v(1:Nh) - xhv(1:Nh)*pkm1v(1:Nh)))

         if (minval(thetav(1:Nh)) <= 0.0_dp .or. maxval(thetav(1:Nh)) >= pi) then
            call err%push_err(ERR_INTERNAL, 'Gauss_Legendre: Newton step left (0, pi)')
            xv = 0.0_dp
            qw = 0.0_dp
            baryv = 0.0_dp
            return
         endif

         maxdiff = maxval(abs(thetav(1:Nh) - toldv(1:Nh)))
      enddo

      !---------------------------------------------------------------
      ! Recompute P_{N-1} at the final roots and build nodes, weights,
      ! and barycentric weights using symmetry.
      !---------------------------------------------------------------
      xhv(1:Nh)  = cos(thetav(1:Nh))
      sinv(1:Nh) = sin(thetav(1:Nh))

      pkm2v(1:Nh) = 1.0_dp
      pkm1v(1:Nh) = xhv(1:Nh)
      do k = 2, N
         pkv(1:Nh)   = ((2*k - 1)*xhv(1:Nh)*pkm1v(1:Nh) - (k - 1)*pkm2v(1:Nh))/k
         pkm2v(1:Nh) = pkm1v(1:Nh)
         pkm1v(1:Nh) = pkv(1:Nh)
      enddo

      do i = 1, Nh
         k = N + 1 - i
         t0 = sinv(i)**2
         t1 = t0/pkm2v(i)

         xv(i) = -xhv(i)
         xv(k) =  xhv(i)
         qw(i) = 2*t0/(dble(N)*N*pkm2v(i)**2)
         qw(k) = qw(i)
         baryv(i) = psign*t1
         baryv(k) = t1
      enddo

      if (mod(N, 2) == 1) then
         mid = Nh + 1
         xv(mid) = 0.0_dp

         t1 = 1.0_dp
         t2 = 0.0_dp
         do k = 2, N - 1
            t0 = -dble(k - 1)*t1/k
            t1 = t2
            t2 = t0
         enddo

         qw(mid) = 2/(dble(N)*N*t2**2)
         baryv(mid) = 1/t2
      endif

   end subroutine get_Gauss_Legendre_grid

   !================================================================
   ! Compute Gauss-Lobatto nodes and optional P_N values.
   !================================================================
   subroutine get_Gauss_Lobatto_grid(N, xv, qw, bary, err)
      implicit none
      integer, intent(in) :: N
      real(dp), dimension(0:N), intent(out) :: xv, qw
      real(dp), dimension(0:N), optional, intent(out) :: bary   ! barycentric weights = 1/P_N
      class(err_type), intent(inout) :: err

      integer :: k, i, iter, maxiter
      integer :: N1, Nh, Ncalc, psign
      real(dp), dimension(0:N) :: xoldv
      real(dp), dimension(0:N) :: pkm2v, pkm1v, pkv, pnv
      real(dp) :: eps, maxdiff

      ! ---- guard against N=0 ----
      if (N == 0) then
         xv(0) = 0.0_dp
         qw(0) = 2.0_dp
         if (present(bary)) bary(0) = 1.0_dp
         return
      end if

      N1 = N + 1
      Nh = (N - 1)/2
      Ncalc = N/2
      if (mod(N, 2) == 0) then
         psign = 1
      else
         psign = -1
      endif

      ! Tolerance
      eps = 10.0_dp*epsilon(1.0_dp)

      maxiter = 50

      ! ---- initial guess ----
      do i = 0, N
         xv(i) = -cos((pi*i)/N)
      end do

      xv(0) = -1.0_dp
      xv(N) =  1.0_dp

      if (Nh > 0) then
         maxdiff = huge(1.0_dp)
         iter = 0
         do while (maxdiff > eps)
            iter = iter + 1
            if (iter > maxiter) then
               call err%push_err(ERR_INTERNAL, 'Gauss_Lobatto: Newton iteration did not converge')
               xv = 0.0_dp
               qw = 0.0_dp
               if (present(bary)) bary = 0.0_dp
               return
            end if

            xoldv(1:Nh) = xv(1:Nh)

            ! Compute only P_N and P_{N-1} on one half of the grid.
            ! Gauss-Lobatto nodes are symmetric, x_i = -x_{N-i}, so the
            ! other half is obtained by reflection. This avoids both the
            ! scalar recurrence over nodes and unnecessary work on paired nodes.
            pkm2v(1:Nh) = 1.0_dp
            pkm1v(1:Nh) = xoldv(1:Nh)
            do k = 2, N
               pkv(1:Nh)   = ((2*k-1)*xoldv(1:Nh)*pkm1v(1:Nh) - (k-1)*pkm2v(1:Nh))/k
               pkm2v(1:Nh) = pkm1v(1:Nh)
               pkm1v(1:Nh) = pkv(1:Nh)
            end do

            ! Update only one half of the interior nodes
            xv(1:Nh) = xoldv(1:Nh) - (xoldv(1:Nh)*pkm1v(1:Nh) - pkm2v(1:Nh)) &
            & / (N1*pkm1v(1:Nh))

            do i = 1, Nh
               xv(N-i) = -xv(i)
            enddo

            maxdiff = maxval(abs(xv(1:Nh) - xoldv(1:Nh)))
         enddo
      endif

      if (mod(N, 2) == 0) xv(N/2) = 0.0_dp

      ! enforce endpoints exactly
      xv(0) = -1.0_dp
      xv(N) =  1.0_dp

      ! ---- recompute P_N at final xv (needed for the weights, and bary if asked) ----
      ! Compute P_N only on one half and restore the other half
      ! using P_N(-x) = (-1)^N P_N(x).
      pkm2v(0:Ncalc) = 1.0_dp
      pkm1v(0:Ncalc) = xv(0:Ncalc)

      do k = 2, N
         pkv(0:Ncalc)   = ((2*k-1)*xv(0:Ncalc)*pkm1v(0:Ncalc) - (k-1)*pkm2v(0:Ncalc))/k
         pkm2v(0:Ncalc) = pkm1v(0:Ncalc)
         pkm1v(0:Ncalc) = pkv(0:Ncalc)
      end do

      pnv(0:Ncalc) = pkm1v(0:Ncalc)
      do i = 0, Ncalc
         pnv(N-i) = psign*pnv(i)
      enddo

      qw(0:N) = 2.0_dp/(dble(N)*(N + 1)*pnv(0:N)**2)
      if (present(bary)) bary(0:N) = 1.0_dp/pnv(0:N)

   end subroutine get_Gauss_Lobatto_grid

   !================================================================
   ! Compute Gauss-Radau-Legendre nodes and quadrature weights.
   !
   ! endpoint = -1 gives nodes including x = -1.
   ! endpoint =  1 gives nodes including x =  1.
   !
   ! For N nodes, the quadrature is exact for polynomials of degree
   ! up to 2*N - 2. If baryv is present, it contains barycentric
   ! weights up to a common factor.
   !================================================================
   subroutine get_Gauss_Radau_grid(N, endpoint, xv, qw, baryv, err)
      implicit none
      integer, intent(in) :: N, endpoint
      real(dp), dimension(1:N), intent(out) :: xv, qw
      real(dp), dimension(1:N), optional, intent(out) :: baryv
      class(err_type), intent(inout) :: err

      integer :: i, iter, maxiter
      real(dp) :: eps, theta, theta_new, step, x, r, dr, dfdtheta
      real(dp) :: xm, xp
      real(dp), dimension(1:N) :: PNm1v, qfacv
      logical :: converged

      if (N < 1) then
         call err%push_err(ERR_INPUT, 'Number of Gauss-Radau nodes must be positive')
         xv = 0.0_dp
         qw = 0.0_dp
         if (present(baryv)) baryv = 0.0_dp
         return
      endif

      if (endpoint /= -1 .and. endpoint /= 1) then
         call err%push_err(ERR_INPUT, 'Gauss-Radau endpoint must be -1 or 1')
         xv = 0.0_dp
         qw = 0.0_dp
         if (present(baryv)) baryv = 0.0_dp
         return
      endif

      eps = 10.0_dp*epsilon(1.0_dp)
      maxiter = 50

      if (endpoint == -1) then
         xv(1) = -1.0_dp
         qw(1) = 2.0_dp/(dble(N)*N)
         qfacv(1) = 2.0_dp
         PNm1v(1) = merge(1.0_dp, -1.0_dp, mod(N - 1, 2) == 0)
         do i = 2, N
            theta = pi - 2.0_dp*pi*dble(i - 1)/dble(2*N - 1)
            converged = .false.

            do iter = 1, maxiter
               x = cos(theta)
               call jacobi_01(N - 1, x, r, dr)
               dfdtheta = -sin(theta)*dr
               if (abs(dfdtheta) <= tiny(1.0_dp)) then
                  call err%push_err(ERR_INTERNAL, 'Gauss_Radau: zero theta derivative')
                  xv = 0.0_dp
                  qw = 0.0_dp
                  if (present(baryv)) baryv = 0.0_dp
                  return
               endif

               step = r/dfdtheta
               theta_new = theta - step

               if (theta_new <= 0.0_dp .or. theta_new >= pi) then
                  call err%push_err(ERR_INTERNAL, 'Gauss_Radau: Newton step left (0, pi)')
                  xv = 0.0_dp
                  qw = 0.0_dp
                  if (present(baryv)) baryv = 0.0_dp
                  return
               endif

               if (abs(theta_new - theta) <= 64.0_dp*eps*max(1.0_dp, abs(theta_new))) then
                  theta = theta_new
                  converged = .true.
                  exit
               endif

               theta = theta_new
            enddo

            if (.not. converged) then
               call err%push_err(ERR_INTERNAL, 'Gauss_Radau: Newton iteration did not converge')
               xv = 0.0_dp
               qw = 0.0_dp
               if (present(baryv)) baryv = 0.0_dp
               return
            endif

            x = cos(theta)
            call jacobi_01(N - 1, x, r, dr)
            xm = 2.0_dp*sin(0.5_dp*theta)**2
            xp = 2.0_dp*cos(0.5_dp*theta)**2

            xv(i) = x
            qw(i) = 4.0_dp/(xp*xp*xm*dr*dr)
            qfacv(i) = xm
            PNm1v(i) = xm*xp*dr/(2.0_dp*dble(N))
         enddo
      else
         do i = 1, N - 1
            theta = 2.0_dp*pi*dble(N - i)/dble(2*N - 1)
            converged = .false.

            do iter = 1, maxiter
               x = cos(theta)
               call jacobi_10(N - 1, x, r, dr)
               dfdtheta = -sin(theta)*dr
               if (abs(dfdtheta) <= tiny(1.0_dp)) then
                  call err%push_err(ERR_INTERNAL, 'Gauss_Radau: zero theta derivative')
                  xv = 0.0_dp
                  qw = 0.0_dp
                  if (present(baryv)) baryv = 0.0_dp
                  return
               endif

               step = r/dfdtheta
               theta_new = theta - step

               if (theta_new <= 0.0_dp .or. theta_new >= pi) then
                  call err%push_err(ERR_INTERNAL, 'Gauss_Radau: Newton step left (0, pi)')
                  xv = 0.0_dp
                  qw = 0.0_dp
                  if (present(baryv)) baryv = 0.0_dp
                  return
               endif

               if (abs(theta_new - theta) <= 64.0_dp*eps*max(1.0_dp, abs(theta_new))) then
                  theta = theta_new
                  converged = .true.
                  exit
               endif

               theta = theta_new
            enddo

            if (.not. converged) then
               call err%push_err(ERR_INTERNAL, 'Gauss_Radau: Newton iteration did not converge')
               xv = 0.0_dp
               qw = 0.0_dp
               if (present(baryv)) baryv = 0.0_dp
               return
            endif

            x = cos(theta)
            call jacobi_10(N - 1, x, r, dr)
            xm = 2.0_dp*sin(0.5_dp*theta)**2
            xp = 2.0_dp*cos(0.5_dp*theta)**2

            xv(i) = x
            qw(i) = 4.0_dp/(xm*xm*xp*dr*dr)
            qfacv(i) = xp
            PNm1v(i) = -xm*xp*dr/(2.0_dp*dble(N))
         enddo
         xv(N) = 1.0_dp
         qw(N) = 2.0_dp/(dble(N)*N)
         qfacv(N) = 2.0_dp
         PNm1v(N) = 1.0_dp
      endif

      if (endpoint == -1) then
         xv(1) = -1.0_dp
         qw(1) = 2.0_dp/(dble(N)*N)
         qfacv(1) = 2.0_dp
         PNm1v(1) = merge(1.0_dp, -1.0_dp, mod(N - 1, 2) == 0)
      else
         xv(N) = 1.0_dp
         qw(N) = 2.0_dp/(dble(N)*N)
         qfacv(N) = 2.0_dp
         PNm1v(N) = 1.0_dp
      endif

      if (present(baryv)) baryv = qfacv/PNm1v

   end subroutine get_Gauss_Radau_grid

   !================================================================
   ! Jacobi polynomial P_n^(1,0)(x) and its derivative.
   ! This special form is used for Radau nodes with x = 1 included.
   !================================================================
   pure subroutine jacobi_10(n, x, p, dpdx)
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: x
      real(dp), intent(out) :: p, dpdx

      integer :: k
      real(dp) :: p0, p1, p2, dp0, dp1, dp2
      real(dp) :: a0, a1, a2, a3, den

      if (n < 0) then
         p    = 0.0_dp
         dpdx = 0.0_dp
         return
      endif

      if (n == 0) then
         p    = 1.0_dp
         dpdx = 0.0_dp
         return
      endif

      p0  = 1.0_dp
      dp0 = 0.0_dp
      p1  = 0.5_dp*(1 + 3*x)
      dp1 = 1.5_dp

      if (n == 1) then
         p    = p1
         dpdx = dp1
         return
      endif

      do k = 1, n - 1
         a0  = 2*k + 2
         a1  = (2*k + 1.0_dp)*(2*k + 3)
         a2  = 1
         a3  = a0*k*(2*k + 3)
         den = a0*(k + 2)*(2*k + 1)

         p2  = (a0*((a1*x + a2)*p1) - a3*p0)/den
         dp2 = (a0*(a1*p1 + (a1*x + a2)*dp1) - a3*dp0)/den

         p0  = p1
         p1  = p2
         dp0 = dp1
         dp1 = dp2
      enddo

      p    = p1
      dpdx = dp1
   end subroutine jacobi_10

   !================================================================
   ! Jacobi polynomial P_n^(0,1)(x) and its derivative.
   ! This special form is used for Radau nodes with x = -1 included.
   !================================================================
   pure subroutine jacobi_01(n, x, p, dpdx)
      implicit none
      integer, intent(in)   :: n
      real(dp), intent(in)  :: x
      real(dp), intent(out) :: p, dpdx

      integer  :: k
      real(dp) :: p0, p1, p2, dp0, dp1, dp2
      real(dp) :: a0, a1, a2, a3, den

      if (n < 0) then
         p    = 0.0_dp
         dpdx = 0.0_dp
         return
      endif

      if (n == 0) then
         p    = 1.0_dp
         dpdx = 0.0_dp
         return
      endif

      p0  = 1.0_dp
      dp0 = 0.0_dp
      p1  = 0.5_dp*(-1 + 3*x)
      dp1 = 1.5_dp

      if (n == 1) then
         p = p1
         dpdx = dp1
         return
      endif

      do k = 1, n - 1
         a0  = 2*k + 2
         a1  = (2*k + 1.0_dp)*(2*k + 3)
         a2  = -1
         a3  = a0*k*(2*k + 3)
         den = a0*(k + 2.0_dp)*(2*k + 1)

         p2  = (a0*((a1*x + a2)*p1) - a3*p0)/den
         dp2 = (a0*(a1*p1 + (a1*x + a2)*dp1) - a3*dp0)/den

         p0  = p1
         p1  = p2
         dp0 = dp1
         dp1 = dp2
      enddo

      p    = p1
      dpdx = dp1
   end subroutine jacobi_01

   !================================================================
   ! Build Legendre polynomial table P_0,...,P_kmax on a given grid.
   ! Used by the eigenvalue filter and the extrapolation series chopping.
   !================================================================
   subroutine create_legendre_matrix(xv, Pm, err)
      implicit none
      real(dp), dimension(:), intent(in) :: xv
      real(dp), allocatable, intent(out) :: Pm(:,:)
      class(err_type), intent(inout) :: err

      integer :: k, kmax, istat

      kmax = size(xv) - 1         ! full transform: one Legendre coeff per node
      if (kmax < 0) then
         call err%push_err(ERR_INTERNAL, 'Wrong size in create_legendre_matrix')
         return
      endif

      allocate(Pm(0:kmax, 0:kmax), stat=istat)
      if (alloc_failed(istat, 'Allocation failed in create_legendre_matrix', err)) return

      Pm(:,0) = 1
      if (kmax >= 1) Pm(:,1) = xv
      do k = 2, kmax
         Pm(:,k) = ((2*k - 1)*xv*Pm(:,k-1) - (k - 1)*Pm(:,k-2))/k
      enddo

   end subroutine create_legendre_matrix

   !================================================================
   ! Compute real Legendre expansion coefficients.
   ! qv must contain quadrature weights matching fv.
   ! i0 is the first row of Pm corresponding to fv(1).
   ! length is the interval the nodes live on (default 2 = reference [-1,1]);
   ! it sets the (2k+1)/length normalization, so the same kernel serves a
   ! reference grid and any stretched/shifted one.
   !================================================================
   subroutine compute_real_legendre_coefficients(fv, Pm, qv, i0, av, err, length)
      real(dp), dimension(:),      intent(in)  :: fv
      real(dp), dimension(0:, 0:), intent(in)  :: Pm
      real(dp), dimension(:),      intent(in)  :: qv
      integer,                    intent(in)  :: i0
      real(dp), dimension(0:),     intent(out) :: av
      class(err_type), intent(inout) :: err
      real(dp),          optional, intent(in)  :: length

      integer  :: i, kmax, nv, q
      real(dp) :: xlen

      nv   = size(fv)
      kmax = ubound(av, 1)

      if (size(qv) /= nv .or. i0 < 0 .or. i0 + nv - 1 > ubound(Pm, 1) .or. kmax > ubound(Pm, 2)) then
         call err%push_err(ERR_INPUT, 'Sizes of array do not match in Legendre expansion coefficients computation')
         av = 0.0_dp
         return
      endif

      xlen = 2.0_dp
      if (present(length)) xlen = length

      av  = 0.0_dp

      do q = 1, nv
         i = i0 + q - 1
         av(0:kmax) = av(0:kmax) + qv(q)*fv(q)*Pm(i,0:kmax)
      enddo

      do i = 0, kmax
         av(i) = (2*i + 1)/xlen*av(i)
      enddo
   end subroutine compute_real_legendre_coefficients

   !================================================================
   ! Grid-aware Legendre projection (real values).
   !
   ! Given a grid (collocation nodes + quadrature weights) sitting on its
   ! own interval [grid%xmin, grid%xmax] -- which may be the reference
   ! [-1,1] or any stretched/shifted interval -- and the function values
   ! fv at the active nodes grid%x(i0:i1), return the coefficients
   ! av(0:i1-i0) of the Legendre expansion on that interval.
   !
   ! The nodes are mapped back onto [-1,1] internally, and the
   ! normalization (2k+1)/(xmax-xmin) makes the coefficients independent
   ! of the interval length. The caller therefore never has to keep a
   ! separate reference grid around: a single (possibly scaled) grid with
   ! its weights is enough.
   !================================================================
   subroutine legendre_coeffs_real(grid, i0, i1, fv, av, err)
      type(grid_type),         intent(in)    :: grid
      integer,                 intent(in)    :: i0, i1
      real(dp), dimension(:),  intent(in)    :: fv
      real(dp), dimension(0:), intent(out)   :: av
      class(err_type),         intent(inout) :: err

      real(dp), allocatable :: Pm(:,:), t(:)
      integer  :: kmax, nv
      real(dp) :: center, half

      if (err%has_error()) return

      kmax = i1 - i0
      nv   = kmax + 1
      if (kmax < 0 .or. size(fv) /= nv .or. ubound(av, 1) < kmax) then
         call err%push_err(ERR_INPUT, 'legendre_coeffs_real: size mismatch')
         av = 0.0_dp
         return
      endif

      center = (grid%xmin + grid%xmax)/2.0_dp
      half   = (grid%xmax - grid%xmin)/2.0_dp

      allocate(t(nv))
      t = (grid%x(i0:i1) - center)/half        ! affine map onto [-1,1]

      call create_legendre_matrix(t, Pm, err)
      if (err%has_error()) return

      ! reuse the shared projection kernel; the grid spans [xmin, xmax], so
      ! its length sets the (2k+1)/(xmax-xmin) normalization
      call compute_real_legendre_coefficients(fv, Pm, grid%qw(i0:i1), 0, av, err, &
                                               length = grid%xmax - grid%xmin)
   end subroutine legendre_coeffs_real

   !================================================================
   ! Grid-aware Legendre projection (complex values).
   ! Complex counterpart of legendre_coeffs_real; see that routine
   ! for the conventions.
   !================================================================
   subroutine legendre_coeffs_cmplx(grid, i0, i1, fcv, acv, err)
      type(grid_type),            intent(in)    :: grid
      integer,                    intent(in)    :: i0, i1
      complex(dp), dimension(:),  intent(in)    :: fcv
      complex(dp), dimension(0:), intent(out)   :: acv
      class(err_type),            intent(inout) :: err

      real(dp), allocatable :: Pm(:,:), t(:)
      integer  :: k, kmax, q, nv
      real(dp) :: center, half

      if (err%has_error()) return

      kmax = i1 - i0
      nv   = kmax + 1
      if (kmax < 0 .or. size(fcv) /= nv .or. ubound(acv, 1) < kmax) then
         call err%push_err(ERR_INPUT, 'legendre_coeffs_cmplx: size mismatch')
         acv = (0.0_dp, 0.0_dp)
         return
      endif

      center = (grid%xmin + grid%xmax)/2.0_dp
      half   = (grid%xmax - grid%xmin)/2.0_dp

      allocate(t(nv))
      t = (grid%x(i0:i1) - center)/half        ! affine map onto [-1,1]

      call create_legendre_matrix(t, Pm, err)
      if (err%has_error()) return

      acv = (0.0_dp, 0.0_dp)
      do q = 1, nv
         acv(0:kmax) = acv(0:kmax) + grid%qw(i0+q-1)*fcv(q)*Pm(q-1, 0:kmax)
      enddo
      do k = 0, kmax
         acv(k) = (2*k + 1)/(2.0_dp*half)*acv(k)  ! 2*half = xmax - xmin
      enddo
   end subroutine legendre_coeffs_cmplx

   !================================================================
   ! Replace |coefficients| in-place with their non-increasing
   ! envelope, normalized to the leading (degree-0) coefficient.
   ! Shared spectral-decay preprocessing for the spurious-state
   ! filter and the extrapolation series chopping.
   !================================================================
   subroutine normalize_decay_envelope(c)
      real(dp), intent(inout) :: c(0:)
      integer :: k, kmax

      kmax = ubound(c, 1)
      do k = kmax - 1, 0, -1
         c(k) = max(c(k), c(k+1))
      enddo
      if (c(0) > 0.0_dp) c = c / c(0)
   end subroutine normalize_decay_envelope

   !================================================================
   ! Spectral-decay test for the spurious-eigenvalue filter.
   ! c holds |coefficients| and is overwritten with the envelope.
   ! Returns .true. if the high-order tail does not decay below thr
   ! (i.e. the state is spurious / a ghost).
   !================================================================
   logical function decay_tail_exceeds(c, thr) result(spurious)
      real(dp), intent(inout) :: c(0:)
      real(dp), intent(in)    :: thr
      integer :: kmax, ntail

      call normalize_decay_envelope(c)

      kmax = ubound(c, 1)
      if (c(0) <= 0.0_dp) then
         spurious = .true.
         return
      endif

      ntail   = min(3, kmax + 1)            ! average over last <=3 coeffs
      spurious = sum(c(kmax-ntail+1:kmax)) / ntail > thr
   end function decay_tail_exceeds

end module
