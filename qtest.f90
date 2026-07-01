!================================================================
! Test program for the qsolver subroutine, which solves coupled-channel
! radial Schrodinger equations for bound and resonant (quasibound) states.
! This program contains many examples that use almost all qsolver capabilities.
! A good practice is to start with the test example closest to the problem
! you want to solve.
!================================================================
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
!================================================================
! List of tests
! Tests 1-2   Particle-in-a-box model with an infinite potential well. Energies and wavefunctions.
! Tests 3-4   Harmonic oscillator. Energies and expectation values.
! Test  5     Eigenvalue problem for Bessel's equation. Regular left boundary.
! Tests 6-8   Eigenvalue problem for Jacobi equation. Eigenvalues, eigenfunctions, and dE/da.
! Tests 9-11  Hydrogen atom s and p states. Energies and expectation values.
! Tests 12-15 Morse potential. Energies, matrix elements, and energy derivatives.
! Test  15    Morse potential: energy derivatives with respect to mass.
! Test  16    Energies of double minimum potential from J. Chem. Phys. 67, 4086 (1977)
! Test  17    Energies of symmetric double minimum potential from J. Chem. Phys. 67, 4086 (1977)
! Test  18    Lennard-Jones 6-12 potential (15 levels). Energy of last bound state.
! Test  19    Lennard-Jones 3-6 potential (301 levels). Energy of last bound state. Extremely difficult case.
! Tests 20-23 Resonant states of 7.5*r^2*exp(-r) potential from J. Phys. B. 7(16), 2189 (1974)
!       20    Find all resonant states (for E > 0) using complex potential
!       21    The same as above but with real potential provided
!       22    Find resonant state closest to target value using complex potential
!       23    The same but with real potential.
! Tests 24-25 Lithium dimer B^1Pi_u state. Resonant levels.
!       24    Compute energies and widths with complex potential
!       25    Compute energies and widths with real potential
! Test  26    Doolen transformed Coulomb + inverse-square potential. 10 resonant states.
! Test  27    Doolen transformed Coulomb + inverse-square potential. Wavefunction comparison.
! Test  28    Doolen transformed Coulomb + inverse-square potential. dE/dgamma and dW/dgamma.
! Test  29    Two-channel Lefebvre model. Complex potential. Bound levels and resonant shifts/widths.
! Test  30    The same as Test 29 but with real potential.
! Test  31    Two-channel Sitnikov-Tolstikhin model (PRA 67, 032714). All three resonances with complex potential.
! Test  32    The same as Test 31 but with real potential.
!
!================================================================
program qtest
use qsolver_mod
use potentials
   implicit none

   !===================================================
   ! Parameters.
   !===================================================
   real(8), parameter :: pi = acos(-1.0d0)
   
   ! Maximum number of grid points
   integer, parameter :: Nmax = 1000

   ! Maximum number of eigenvalues to be found
   integer, parameter :: nEmax = 700

   !================================================================
   ! Scalars
   !================================================================

   ! number of eigenvalues requested and found
   integer :: nE, nEf

   ! number of grid points
   integer :: N 

   ! effective mass
   real(8) :: mu

   ! tolerance
   real(8) :: tol

   ! max error
   real(8) :: max_err

   ! number of passed and failed tests
   integer :: npass, nfail

   ! output unit for qtest_out.txt
   integer :: uout   

   ! I/O status code
   integer :: ios    

   ! I/O message text
   character(len=256) :: iomsg  

   ! temporary variables
   real(8) :: t0, t1, t2, t3, t4
   ! indices
   integer :: i, j

   ! number of missing levels in target-state searches
   integer :: missing_levels

   !================================================================
   ! Arrays
   !================================================================

   ! arrays with calculated energies and resonance widths
   real(8) :: E0v(nEmax), Ev(nEmax), Wv(nEmax)

   ! array with reference energies and widths
   real(8) :: Eref(nEmax), Wref(nEmax)
   real(8) :: Elit(nEmax), Wlit(nEmax)

   ! array with experimental energies and widths
   real(8) :: Eexp(nEmax), Wexp(nEmax)

   ! array with experimental rotational quantum numbers
   integer :: Jv(nEmax)

   character(len=10) :: sElit, sWlit

   ! arrays with calculated energy and width derivatives: dE(p, i), dW(p, i)
   real(8) :: dEv(10, nEmax), dWv(10, nEmax)

   ! arrays with calculated wavefunctions for the single-equation case
   real(8)    :: WF1(Nmax, nEmax)
   real(8)    :: braW1(Nmax, nEmax)
   complex(8) :: WF1c(Nmax, nEmax)

   ! array with grid points
   real(8)    :: rv(Nmax)
   complex(8) :: rvc(Nmax)

   ! array with quadrature weights
   real(8)    :: qw(Nmax)
   complex(8) :: qwc(Nmax)

   !===================================================
   ! Initialize counters of passed and failed tests
   !===================================================
   npass = 0
   nfail = 0

   uout = 17
   open(unit = uout, file = 'qtest_out.txt', status = 'replace', action = 'write', &
      & iostat = ios, iomsg = iomsg)
   if (ios /= 0) then
      write(*,*) 'ERROR: cannot open qtest_out.txt'
      write(*,*) trim(iomsg)
      stop 1
   endif

   write(*,*) ''
   write(*,*) 'qsolver version '//trim(qsolver_version)
   write(*,*) 'Detailed test results are written to qtest_out.txt'
   write(*,*) ''

   !===================================================
   ! Test 1-2
   ! Particle-in-a-box model with an infinite potential well
   ! Testing for energies and wavefunctions.
   ! Only the necessary parameters should be specified.
   ! If no potential is provided, it is assumed that V = 0.
   !===================================================
   write(uout,*) ''
   write(uout,*) 'Test 1-2. Finding energies and wavefunctions of particle in box model'
   call announce_test('Test 1. Particle in a box: energies')

   N  = 40 ! use N grid points
   nE = 10 ! request 10 lowest eigenvalues
   
   call qsolver( rmin = 0d0,   & ! the problem is solved on [0; pi]
          &      rmax = pi,    & !
          &      mu   = 0.5d0, & ! effective mass
          &      N    = N,     & ! number of grid points
          &      nE   = nE,    & ! number of eigenvalues to be found
          &      E    = Ev,    & ! output array with energies
          &      WF1  = WF1,   & ! output array with wavefunctions
          &      rv   = rv,    & ! output array with grid points
          &      qw   = qw     & ! output array with quadrature weights for L2 norm calculation
          &                )

   !===================================================
   ! Set the maximum tolerance for energies and eigenvectors
   !===================================================
   tol = 5.0d-12
   
   !===================================================
   ! Check energies by comparing with the exact values E0_i = i**2
   ! against the calculated ones
   !===================================================
   write(uout,*) ''
   write(uout,*) 'Checking energies..'
   write(uout,*) ''
   write(uout,'(A4, 2A21, A15)') '', 'Calculated energy', 'Exact energy', 'Relative error'

   max_err = 0.0d0
   do i = 1, nE
      associate( E0i      => dble(i**2),  & ! exact energy
            &    Ei       => Ev(i),       & ! calculated energy
            &    rel_err  => t1)

         rel_err = abs(Ei - E0i)/abs(E0i)
         max_err = max(max_err, rel_err)    ! max relative error
         write(uout,'(I4, 2F21.14, ES15.2)') i, Ei, E0i, rel_err
      end associate
   enddo 

   !===================================================
   ! If max_err > tol, then the test fails
   !===================================================
   call check_err(max_err, tol)

   !===================================================
   ! Check wavefunctions by comparing with the exact formula psi(i) = sqrt(2/pi)*sin(i*r)
   !===================================================
   call announce_test('Test 2. Particle in a box: wavefunctions')
   write(uout,*) 'Checking wavefunctions..'
   write(uout,*) 'Computing max and L2 - norms between exact and calculated wavefunctions.'
   write(uout,*) ''
   write(uout,'(A4, 2A15)') '', 'max-norm', 'L2-norm'

   max_err = 0.0d0
   
   do i = 1, nE
      !===================================================
      ! Compute the max norm and L2 norm of the difference between WFi and WF0i
      ! WFi - calculated eigenfunction, WF0i - exact eigenfunction
      !===================================================

      associate( WFi => WF1(1:N, i),                & ! calculated eigenfunction
                 WF0i => sqrt(2/pi)*sin(i*rv(1:N)), & ! exact eigenfunction
                 w    => qw(1:N),                   & ! quadrature weights
                 max_norm => t1,                    & 
                 L2_norm  => t2,                    &
                 s        => j)

         ! qsolver returns normalized wavefunctions,
         ! but they are defined up to a factor of +/-1
         s = merge(1, -1, WF1(1, i) > 0.0d0)          ! get the wavefunction sign
         max_norm  = maxval(abs(WFi - s*WF0i))        ! max-norm
         L2_norm   = sqrt(sum((WFi - s*WF0i)**2 * w)) ! L2 norm
         max_err   = max(max_err, max_norm, L2_norm)
         write(uout, '(I4, 2ES15.2)') i, max_norm, L2_norm
      end associate
   enddo
   call check_err(max_err, tol)

   !===================================================
   ! Tests 3-4
   ! Harmonic oscillator V = m*w^2*r^2/2 for w = 1, m = 1, hbar = 1
   ! Testing for energies and expectation values.
   ! The equation is solved on the interval (-infinity; +infinity).
   ! rmin and rmax must be omitted in this case.
   ! The interval is unbounded; therefore a mapping should be used to transform it to the interval
   ! (-1; 1). The best mapping for this case appeared to be the double-exponential transformation.
   ! The easiest way to obtain the optimal mapping parameter is by trial and error.
   !===================================================
   write(uout,*) ''
   write(uout,*) 'Tests 3-4. Harmonic oscillator'
   call announce_test('Test 3. Harmonic oscillator: energies')

   N  = 45 ! use N grid points
   nE = 10 ! request 10 lowest eigenvalues

   call qsolver( mu    = 1.0d0,      & ! effective mass
          &      N     = N,          & ! number of grid points
          &      nE    = nE,         & ! number of eigenvalues to be found
          &      V1    = V_Harm_Osc, & ! V1 = 1/2*r*r
          &      mtype = 'dexp',     & ! If the interval is unbounded, then a mapping should be used.
          &      beta  = 0.2d0,      & ! Double-exponential mapping works very well here.
          &      E     = Ev,         & ! output array with energies
          &      WF1   = WF1,        & ! output array with wavefunctions
          &      braW1 = braW1,      & ! output weighted bra vectors
          &      rv    = rv          & ! output array with grid points
          &                )

   !===================================================
   ! Set the maximum tolerance for energies and eigenvectors
   !===================================================
   tol = 1.0d-13
   
   !===================================================
   ! Check energies by comparing with the exact values E0_i = hbar*w*(n + 1/2)
   ! against the calculated ones.
   !===================================================
   write(uout,*) ''
   write(uout,*) 'Checking energies..'
   write(uout,*) ''
   write(uout,'(A4, 2A21, A15)') '', 'Calculated energy', 'Exact energy', 'Relative error'

   max_err = 0.0d0
   do i = 1, nE
      associate( E0i      => i - 0.5d0,   & ! exact energy
            &    Ei       => Ev(i),       & ! calculated energy
            &    rel_err  => t1)

         rel_err = abs(Ei - E0i)/abs(E0i)
         max_err = max(max_err, rel_err)    ! max relative error
         write(uout,'(I4, 2F21.14, ES15.2)') i, Ei, E0i, rel_err
      end associate
   enddo 

   !===================================================
   ! If max_err > tol, then the test fails
   !===================================================
   call check_err(max_err, tol)

   !===================================================
   ! Check expectation values <psi_n|r**2|psi_n> by comparing
   ! with the exact formula hbar/(m*w)*(n + 1/2)
   !===================================================
   call announce_test('Test 4. Harmonic oscillator: expectation values')
   write(uout,*) 'Checking expectation values..'
   write(uout,*) ''
   write(uout,'(A4, 2A21, A15)') '', '<psi_i|r^2|psi_i>', 'Exact value', 'Relative error'

   max_err = 0.0d0
   
   do i = 1, nE
      !===================================================
      ! Compute the r^2 expectation value
      !===================================================

      associate( braWi => braW1(1:N, i), & ! weighted bra vector
            &    ketWi => WF1(1:N, i),   & ! wavefunction (ket vector)
            &    rv   => rv(1:N),        & ! grid points
            &    xsq  => t1,             & ! <psi|r^2|psi>
            &    xsq0 => (i-1) + 0.5d0,  & ! exact value
            &    rel_err  => t2         )  ! relative error

         xsq     = dot_product(braWi, rv**2*ketWi) ! calculation of the integral
         rel_err = abs(xsq - xsq0)/abs(xsq0)
         max_err = max(max_err, rel_err)
         write(uout,'(I4, 2F21.14, ES15.2)') i, xsq, xsq0, rel_err
      end associate
   enddo
   call check_err(max_err, tol)

   !=========================================================
   ! Test 5. Eigenvalue problem for Bessel's equation:
   !
   !   d/dx ( x*y'(x) ) + lambda*x*y(x) = 0,
   !   0 <= x <= 1,
   !
   ! with boundary conditions
   !   |y(0)| < infinity,
   !    y(1) = 0.
   !
   ! Physically this is the vibrating circular membrane (drumhead) with a
   ! clamped rim: it is the radial part of the Helmholtz equation for the
   ! axially symmetric (m = 0) modes. The eigenvalues lambda_n fix the
   ! natural frequencies and the eigenfunctions y(x) = J_0(sqrt(lambda)*x)
   ! give the mode shapes.
   !
   ! For x > 0 this is equivalent to
   !
   !   -y''(x) - (1/x)*y'(x) = lambda*y(x).
   !
   ! In qsolver notation this is obtained with
   !   mu = 0.5,
   !   B(x) = -1/x,
   !   is_B_sym = .false.
   !
   ! The left endpoint is singular. The physical solution is regular
   ! there, so bc_left = BC_REGULAR is used. The right endpoint has the
   ! Dirichlet condition y(1) = 0.
   !
   ! Exact eigenvalues are
   !   lambda_n = j_{0,n}**2,
   ! where j_{0,n} is the n-th positive zero of J_0(x).
   !=========================================================
   write(uout,*) ''
   write(uout,*) 'Test 5. Eigenvalue problem for Bessel''s equation'
   call announce_test('Test 5. Eigenvalue problem for Bessel''s equation')
   write(uout,*) '         Compare first 10 eigenvalues with zeros of J0'

   N = 30
   tol = 1.0d-11
   max_err = 0.0d0

   Eref(1:10) = (/ &
              &    5.783185962946785d0, 30.471262343662087d0, 74.88700679069518d0, &
              &    139.04028442645983d0, 222.93230361763415d0, 326.5633529323284d0, &
              &    449.9335285180356d0, 593.042869655955d0, 755.891394783933d0, &
              &    938.479113475694d0 /)

   call qsolver( rmin     = 0.0d0,         & ! singular regular endpoint
          &      rmax     = 1.0d0,         & ! y(1) = 0
          &      bc_left  = BC_REGULAR,    & ! finite y(0) is allowed
          &      mu       = 0.5d0,         & ! gives kinetic term -d2/dx2
          &      N        = N,             & ! number of collocation points
          &      nE       = 10,            & ! first 10 eigenvalues
          &      nEf      = nEf,           & ! number of eigenvalues found
          &      B1       = B_Bessel,      & ! B(x) = -1/x
          &      is_B_sym = .false.,       & ! B enters as plain B*d/dx
          &      E        = Ev             & ! output eigenvalues
          &               )

   write(uout,*) ''
   write(uout,'(A4, 2A20, A12)') 'n', 'Calc. lambda', 'Exact lambda', 'Rel. diff.'

   if (nEf < 10) then
      write(uout,*) 'TEST FAILED: less than 10 eigenvalues found'
      max_err = 2*tol
   endif

   do i = 1, min(10, nEf)
      t1 = abs(Ev(i) - Eref(i))/abs(Eref(i))
      max_err = max(max_err, t1)
      write(uout,'(I4, 2F20.10, ES12.2)') i, Ev(i), Eref(i), t1
   enddo
   call check_err(max_err, tol)

   !=========================================================
   ! Tests 6-8. Eigenvalue problem for Jacobi equation at a = 0.
   !
   !   -(1 - x^2)*y''(x) + 2*(a + 1)*x*y'(x) = E*y(x),
   !   -1 < x < 1,
   !
   ! with regular boundary conditions
   !   |y(-1)| < infinity,
   !   |y( 1)| < infinity.
   !
   ! After division by (1 - x^2), the equation has qsolver form
   !
   !   -y''(x) + B(x)*y'(x) = E*G(x)*y(x),
   !
   ! where
   !   B(x) = 2*(a + 1)*x/(1 - x^2),
   !   G(x) = 1/(1 - x^2).
   !
   ! Here a = 0, so the exact eigenfunctions are Legendre polynomials
   !   y_n(x) = P_n(x),
   ! and the exact eigenvalues are
   !   E_n = n*(n + 1),  n = 0,1,2,...
   !
   ! The derivative with respect to a is also known analytically:
   !   dE_n/da = 2*n.
   !=========================================================
   write(uout,*) ''
   write(uout,*) 'Tests 6-8. Eigenvalue problem for Jacobi equation'
   call announce_test('Test 6. Jacobi equation: eigenvalues')

   N  = 10
   nE = 10
   tol = 1.0d-10
   max_err = 0.0d0

   call qsolver( rmin     = -1.0d0,      & ! regular left endpoint
          &      rmax     =  1.0d0,      & ! regular right endpoint
          &      bc_left  = BC_REGULAR,  & ! finite y(-1) is allowed
          &      bc_right = BC_REGULAR,  & ! finite y( 1) is allowed
          &      mu       = 0.5d0,       & ! gives kinetic term -d2/dx2
          &      N        = N,           & ! number of collocation points
          &      nE       = nE,          & ! first 10 eigenvalues
          &      nEf      = nEf,         & ! number of eigenvalues found
          &      B1       = B_Jacobi,    & ! B(x) = 2*x/(1 - x^2)
          &      G1       = G_Jacobi,    & ! G(x) = 1/(1 - x^2)
          &      is_B_sym = .false.,     & ! B enters as plain B*d/dx
          &      E        = Ev,          & ! output eigenvalues
          &      WF1      = WF1,         & ! output right eigenfunctions
          &      rv       = rv,          & ! output grid points
          &      qw       = qw           & ! output quadrature weights
          &               )

   write(uout,*) ''
   write(uout,'(A4, 2A20, A12)') 'n', 'Calc. E', 'Exact E', 'Abs. diff.'

   if (nEf < nE) then
      write(uout,*) 'TEST FAILED: less than 10 eigenvalues found'
      max_err = 2*tol
   endif

   do i = 1, min(nE, nEf)
      associate( nlev       => i - 1,          &
            &    E0i     => dble((i - 1)*i), &
            &    abs_err => t1 )

         abs_err = abs(Ev(i) - E0i)
         max_err = max(max_err, abs_err)
         write(uout,'(I4, 2F20.10, ES12.2)') nlev, Ev(i), E0i, abs_err
      end associate
   enddo
   call check_err(max_err, tol)

   !=========================================================
   ! Test 7. Compare calculated eigenfunctions with exact
   ! Legendre polynomials P_n(x), normalized by
   !
   !   integral_-1^1 [sqrt((2*n+1)/2)*P_n(x)]^2 dx = 1.
   !=========================================================
   call announce_test('Test 7. Jacobi equation: eigenfunctions')
   write(uout,*) ''
   write(uout,*) 'Checking eigenfunctions against normalized Legendre polynomials.'
   write(uout,*) ''
   write(uout,'(A4, 2A15)') 'n', 'max-norm', 'L2-norm'

   tol = 1.0d-8
   max_err = 0.0d0

   do i = 1, min(nE, nEf)
      associate( nlev        => i - 1,                                       &
            &    WF       => WF1(1:N, i),                                  &
            &    WF0      => sqrt((2.0d0*(i - 1) + 1.0d0)/2.0d0)           &
            &                *legendre_poly(i - 1, rv(1:N)),              &
            &    QW       => qw(1:N),                                      &
            &    max_norm => t1,                                          &
            &    L2_norm  => t2,                                          &
            &    s        => t3 )

         s = merge(1.0d0, -1.0d0, WF(1)/WF0(1) > 0.0d0)
         max_norm = maxval(abs(WF - s*WF0))
         L2_norm  = sqrt(sum((WF - s*WF0)**2*QW))
         max_err  = max(max_err, max_norm, L2_norm)
         write(uout,'(I4, 2ES15.2)') nlev, max_norm, L2_norm
      end associate
   enddo
   call check_err(max_err, tol)

   !=========================================================
   ! Test 8. Check derivatives dE/da for the same Jacobi
   ! equation. Since E_n(a)=n*(n+2*a+1),
   !
   !   dE_n/da = 2*n.
   !=========================================================
   call announce_test('Test 8. Jacobi equation: energy derivatives')
   write(uout,*) ''
   write(uout,*) 'Checking dE/da against exact values 2*n.'
   write(uout,*) ''
   write(uout,'(A4, 2A15, A12)') 'n', 'dE/da calc', 'dE/da exact', 'Abs. diff.'

   tol = 1.0d-8
   max_err = 0.0d0

   call qsolver( rmin     = -1.0d0,      & ! regular left endpoint
          &      rmax     =  1.0d0,      & ! regular right endpoint
          &      bc_left  = BC_REGULAR,  & ! finite y(-1) is allowed
          &      bc_right = BC_REGULAR,  & ! finite y( 1) is allowed
          &      mu       = 0.5d0,       & ! gives kinetic term -d2/dx2
          &      N        = N,           & ! number of collocation points
          &      nE       = nE,          & ! first 10 eigenvalues
          &      nEf      = nEf,         & ! number of eigenvalues found
          &      B1       = B_Jacobi,    & ! B(x) = 2*x/(1 - x^2)
          &      G1       = G_Jacobi,    & ! G(x) = 1/(1 - x^2)
          &      dB1      = dB_Jacobi,   & ! dB/da = 2*x/(1 - x^2)
          &      dE       = dEv,         & ! output dE/da
          &      is_B_sym = .false.,     & ! B enters as plain B*d/dx
          &      E        = Ev           & ! output eigenvalues
          &               )

   do i = 1, min(nE, nEf)
      associate( nlev       => i - 1,      &
            &    dE0     => dble(2*(i - 1)), &
            &    abs_err => t1 )

         abs_err = abs(dEv(1, i) - dE0)
         max_err = max(max_err, abs_err)
         write(uout,'(I4, 2ES15.6, ES12.2)') nlev, dEv(1, i), dE0, abs_err
      end associate
   enddo
   call check_err(max_err, tol)

   !===================================================
   ! Test 9.
   ! Hydrogen atom s states. Atomic units are used. The potential is -1/r.
   ! The equation is solved on the interval (0; +infinity).
   ! rmax must be omitted in this case.
   ! The atan mapping appears to be best here. The optimal beta was obtained by hand.
   !===================================================
   write(uout,*) ''
   write(uout,*) 'Tests 9-11. Hydrogen atom'
   call announce_test('Test 9. Hydrogen atom: s-state energies')

   N  = 55 ! use N grid points
   nE = 10 ! request 10 lowest eigenvalues

   call qsolver( rmin  = 0.0d0,      & ! the equation is solved on (0; +infinity)
          &      mu    = 1.0d0,      & ! effective mass
          &      N     = N,          & ! number of grid points
          &      nE    = nE,         & ! number of eigenvalues to be found
          &      V1    = V_Coulomb,  & ! V1 = -1/r
          &      mtype = 'atan',     & ! If the interval is unbounded, then a mapping should be used.
          &      beta  = 0.01d0,     & ! best beta value from trial and error.
          &      E     = Ev,         & ! output array with energies
          &      WF1   = WF1,        & ! output array with wavefunctions
          &      braW1 = braW1,      & ! output weighted bra vectors
          &      rv    = rv          & ! output array with grid points
          &                )

   !===================================================
   ! Set the maximum tolerance for energies and eigenvectors
   !===================================================
   tol = 1.0d-13
   
   !===================================================
   ! Checking energies by comparing exact values -1/(2n^2)
   ! against the calculated ones.
   !===================================================
   write(uout,*) ''
   write(uout,*) 'Checking energies of hydrogen s states..'
   write(uout,*) ''
   write(uout,'(A4, 2A21, A15)') '', 'Calculated energy', 'Exact energy', 'Relative error'

   max_err = 0.0d0
   do i = 1, nE
      associate( E0i      => -0.5d0/(i**2), & ! exact energy
            &    Ei       => Ev(i),         & ! calculated energy
            &    rel_err  => t1)

         rel_err = abs(Ei - E0i)/abs(E0i)
         max_err = max(max_err, rel_err)    ! max relative error
         write(uout,'(I4, 2F21.14, ES15.2)') i, Ei, E0i, rel_err
      end associate
   enddo 

   !===================================================
   ! If max_err > tol, then the test fails
   !===================================================
   call check_err(max_err, tol)

   !===================================================
   ! Test 10-11
   ! Hydrogen atom p states. Atomic units are used. The potential is -1/r.
   ! Here the keyword Jrot = 1 is used, which adds the rotational part to the
   ! potential.
   ! Testing for energies and expectation values.
   !===================================================

   N  = 60 ! use N grid points
   nE = 10 ! request 10 lowest eigenvalues
   call announce_test('Test 10. Hydrogen atom: p-state energies')

   call qsolver( rmin  = 0.0d0,      & ! the equation is solved on (0; +infinity)
          &      mu    = 1.0d0,      & ! effective mass
          &      N     = N,          & ! number of grid points
          &      nE    = nE,         & ! number of eigenvalues to be found
          &      V1    = V_Coulomb,  & ! V1 = -1/r
          &      Jrot  = 1,          & ! add J*(J+1)/r^2 to V1
          &      mtype = 'atan',     & ! If the interval is unbounded, then a mapping should be used.
          &      beta  = 0.01d0,     & ! best beta value from trial and error.
          &      E     = Ev,         & ! output array with energies
          &      WF1   = WF1,        & ! output array with wavefunctions
          &      braW1 = braW1,      & ! output weighted bra vectors
          &      rv    = rv          & ! output array with grid points
          &                )

   !===================================================
   ! Set the maximum tolerance for energies and eigenvectors
   !===================================================
   tol = 1.0d-13
   
   !===================================================
   ! Checking energies by comparing exact values -1/(2(n+1)^2)
   ! against the calculated ones.
   !===================================================
   write(uout,*) ''
   write(uout,*) 'Checking energies of hydrogen p states..'
   write(uout,*) ''
   write(uout,'(A4, 2A23, A15)') '', 'Calculated energy [au]', 'Exact energy [au]', 'Relative error'

   max_err = 0.0d0
   do i = 1, nE
      associate( E0i      => -0.5d0/((i+1)**2), & ! exact energy
            &    Ei       => Ev(i),             & ! calculated energy
            &    rel_err  => t1)

         rel_err = abs(Ei - E0i)/abs(E0i)
         max_err = max(max_err, rel_err)    ! max relative error
         write(uout,'(I4, 2F23.14, ES15.2)') i, Ei, E0i, rel_err
      end associate
   enddo 

   !===================================================
   ! If max_err > tol, then the test fails
   !===================================================
   call check_err(max_err, tol)

   !===================================================
   ! Check expectation values <psi_n|r|psi_n> (average electron-nucleus distance)
   ! by comparing with the exact formula 1/2*(3n^2 - l*(l+1)) where l = 1 (p-states)
   !===================================================
   call announce_test('Test 11. Hydrogen atom: p-state average radius')
   write(uout,*) 'Checking calculated average electron-nucleus distances for hydrogen p states..'
   write(uout,*) ''
   write(uout,'(A4, 2A23, A15)') '', '<psi_i|r|psi_i> [au]', 'Exact value [au]', 'Relative error'

   max_err = 0.0d0
   
   do i = 1, nE
      !===================================================
      ! Compute the expectation value of r
      !===================================================

      associate( Wi   => WF1(1:N, i),            & ! calculated eigenfunction
            &    braWi => braW1(1:N, i),         & ! weighted bra vector
            &    rv   => rv(1:N),                & ! grid points
            &    rav  => t1,                     & ! <psi|r|psi>
            &    rav0 => 0.5d0*(3*(i+1)**2 - 2), & ! exact value
            &    rel_err  => t2         )          ! relative error

         rav     = dot_product(braWi, rv*Wi) ! calculation of the integral
         rel_err = abs(rav - rav0)/abs(rav0)
         max_err = max(max_err, rel_err)
         write(uout,'(I4, 2F23.14, ES15.2)') i, rav, rav0, rel_err
      end associate
   enddo
   call check_err(max_err, tol)

   !===================================================
   ! Test 12-15.
   ! Morse potential with 50 bound levels.
   ! mu = 1, hbar = 1. a = 2*sqrt(2), De = 10000, re = 3.
   ! Tests of energies, matrix elements, and energy derivatives.
   ! The equation is solved on the interval (1; +infinity).
   ! In this case the number of grid points can be slightly reduced by increasing
   ! rmin from 0 because the wavefunction dies off rapidly in the repulsive region
   ! and quickly becomes negligibly small. Therefore rmin can be
   ! from 0 to about 1.5 without loss of accuracy.
   !
   ! The rmax keyword is omitted to make rmax = infinity.
   !
   ! The two-parameter atan mapping appears to be best here. The optimal beta was
   ! obtained by hand to minimize the number of grid points.
   !===================================================
   write(uout,*) ''
   write(uout,*) 'Tests 12-15. Morse potential'
   call announce_test('Test 12. Morse potential: energies')

   N  = 170 ! use N grid points
   nE = 50  ! request 50 lowest eigenvalues

   call qsolver( rmin  = 1.0d0,      & ! the equation is solved on (1; +infinity)
          &      mu    = 1.0d0,      & ! effective mass
          &      N     = N,          & ! number of grid points
          &      nE    = nE,         & ! number of eigenvalues to be found
          &      V1    = V_Morse,    & ! V1 = 50-level Morse potential
          &      mtype = 'atan',     & ! If the interval is unbounded, then a mapping should be used.
          &      beta  = 5.0d0,      & ! good beta value obtained by hand
          &      ye    = 3.0d0,      & ! best ye value is the value of r where V has its minimum
          &      E     = Ev,         & ! output array with energies
          &      WF1   = WF1,        & ! output array with wavefunctions
          &      braW1 = braW1,      & ! output weighted bra vectors
          &      rv    = rv,         & ! output array with grid points
          &      dV1   = dV_Morse,   & ! derivative subroutine for Morse potential
          &      dE    = dEv,        & ! calculated derivatives of energies
          &      ndE   = 2           & ! number of derivatives
          &                )

   !===================================================
   ! Set the maximum tolerance for energies
   !===================================================
   tol = 5.0d-12
   
   !===================================================
   ! Check energies by comparing with the exact values 10000 - (101 - 2*i)**2
   ! against the calculated ones.
   !===================================================
   write(uout,*) ''
   write(uout,*) 'Checking energies of 50-level Morse potential.'
   write(uout,*) ''
   write(uout,'(A4, 2A21, A15)') '', 'Calculated energy', 'Exact energy', 'Relative error'

   max_err = 0.0d0
   do i = 1, nE
      associate( E0i      => 10000d0 - (101 - 2*i)**2, & ! exact energy
            &    Ei       => Ev(i),                    & ! calculated energy
            &    rel_err  => t1)

         rel_err = abs(Ei - E0i)/abs(E0i)
         max_err = max(max_err, rel_err)    ! max relative error
         write(uout,'(I4, 2F21.14, ES15.2)') i, Ei, E0i, rel_err
      end associate
   enddo 

   !===================================================
   ! If max_err > tol, then the test fails
   !===================================================
   call check_err(max_err, tol)

   !===================================================
   ! Set the maximum tolerance for energies and eigenvectors
   !===================================================
   tol = 5.0d-8
    
   !===================================================
   ! Test 13
   ! Check the calculated matrix elements <psi_m|r|psi_n>
   ! by comparing with the exact formula given for this model
   !===================================================
   call announce_test('Test 13. Morse potential: <r> matrix elements')
   write(uout,*) 'Checking calculated matrix elements of r..'
   write(uout,*) ''
   write(uout,'(A4, 2A23, A15)') '', '<psi_i|r|psi_i+1>', 'Exact value', 'Relative error'

   max_err = 0.0d0
   do i = 1, nE-1
      !===================================================
      ! Compute matrix elements of the r operator
      !===================================================

      associate( WFip1 => WF1(1:N, i+1), & ! calculated eigenfunction i+1
            &    braWi => braW1(1:N, i), & ! weighted bra vector for state i
            &    rv    => rv(1:N),       & ! grid points
            &    rmn   => t1,            & ! <psi_i|r|psi_i+1>
            &    rmn0  => 1.0d0/(2*sqrt(2.0d0)*(50-i))                     &
            &             *sqrt((50.5d0 - i)*(49.5d0 - i)*i/(100.d0 - i)), & ! exact value
            &    rel_err  => t2         )           ! relative error

         rmn     = abs(dot_product(braWi, rv*WFip1)) ! calculation of the integral
         rel_err = abs(rmn - rmn0)/abs(rmn0)
         max_err = max(max_err, rel_err)
         write(uout,'(I4, 2F23.14, ES15.2)') i, rmn, rmn0, rel_err
      end associate
   enddo
   call check_err(max_err, tol)

   !===================================================
   ! Test 14
   ! Checking dE/dDe and dE/da of the 50-level Morse potential
   !
   ! Potential: V(r) = De * (1 - exp(-a*(r-re)))^2, min = 0, asymptote = De
   ! Energies (our special integer-spectrum choice): E_i = De - (101 - 2*i)^2, i = 1..50
   !
   ! Exact derivatives at this parameter point:
   !   dE_i/dDe = (2*i - 1)/100
   !   dE_i/da = (sqrt(2)/2) * (2*i - 1) * (101 - 2*i)
   !
   !===================================================
   write(uout,*) ''
   call announce_test('Test 14. Morse potential: energy derivatives')
   write(uout,*) 'Checking dE/dDe and dE/da of 50-level Morse potential.'
   write(uout,*) ''
   write(uout,'(A4, 2A11, A14, A11)') '', 'dE/dDe', 'rel_err_De', 'dE/da', 'rel_err_a'
   
   max_err = 0.0d0
   
   do i = 1, nE
      associate( dEi_De  => dEv(1, i), &
               & dEi_a   => dEv(2, i), &
               & dE0i_De => (2*i - 1) / 100.0d0, &
               & dE0i_a  => 0.5d0 * sqrt(2.0d0) * (2*i - 1) * (101 - 2*i), &
               & rel_De  => t1, &
               & rel_a   => t2 )
   
         rel_De = abs(dEi_De - dE0i_De)/abs(dE0i_De)
         rel_a  = abs(dEi_a - dE0i_a)/abs(dE0i_a)
   
         max_err = max(max_err, rel_De, rel_a)
   
         write(uout,'(I4, 1X, F10.6, 1X, ES10.2, 1X, F13.6, 1X, ES10.2)') &
              i, dEi_De, rel_De, dEi_a, rel_a
      end associate
   end do
   call check_err(max_err, tol)

   !===================================================
   ! Test 15
   ! Checking dE/dM of the 50-level Morse potential.
   !
   ! Physical equation:
   !   [-(1/(2*M)) d2/dr2 + V(r)] psi = E(M) psi.
   !
   ! To test derivatives with respect to the mass-like parameter M using
   ! the generalized matrix G, the same equation is written as
   !   [-1/2 d2/dr2 + M*V(r)] psi = E(M)*M*psi.
   ! At M = 1 this is the same problem as in Tests 12-14.
   ! Therefore dV/dM = V(r) and dG/dM = 1.
   !
   ! For the Morse potential
   !   E_n(M) = a*sqrt(2*De/M)*(n + 1/2)
   !          - a**2*(n + 1/2)**2/(2*M),  n = 0..49.
   !
   ! With De = 10000, a = 2*sqrt(2), M = 1, and i = n + 1:
   !   dE_i/dM = (2*i - 1)*(2*i - 101).
   !===================================================
   write(uout,*) ''
   write(uout,*) 'Test 15. Morse potential: energy derivatives with respect to mass.'
   call announce_test('Test 15. Morse potential: energy derivatives with respect to mass')
   write(uout,*) 'Checking dE/dM of 50-level Morse potential.'
   write(uout,*) ''
   write(uout,'(A4, 2A18, A12)') '', 'dE/dM calc', 'dE/dM exact', 'Rel. error'

   call qsolver( rmin  = 1.0d0,           & ! the equation is solved on (1; +infinity)
          &      mu    = 1.0d0,           & ! kinetic term in the mass-scaled equation
          &      N     = N,               & ! number of grid points
          &      nE    = nE,              & ! number of eigenvalues to be found
          &      nEf   = nEf,             & ! number of eigenvalues found
          &      V1    = V_Morse,         & ! V1 = M*V_Morse at M = 1
          &      G1    = G_Morse_mass,    & ! G1 = M at M = 1
          &      mtype = 'atan',          & ! same mapping as in Tests 12-14
          &      beta  = 5.0d0,           & ! same beta as in Tests 12-14
          &      ye    = 3.0d0,           & ! same mapping center as in Tests 12-14
          &      E     = Ev,              & ! output array with energies
          &      dV1   = dV_Morse_mass,   & ! d(M*V)/dM = V
          &      dG1   = dG_Morse_mass,   & ! dG/dM = 1
          &      dE    = dEv,             & ! calculated derivatives of energies
          &      ndE   = 1                & ! one derivative with respect to M
          &                )

   tol = 5.0d-8
   max_err = 0.0d0

   if (nEf < nE) then
      write(uout,*) 'TEST FAILED: less than 50 eigenvalues found'
      max_err = 2.0d0*tol
   endif

   do i = 1, min(nE, nEf)
      associate( dEi_M  => dEv(1, i), &
               & dE0i_M => dble((2*i - 1)*(2*i - 101)), &
               & rel_M  => t1 )

         rel_M = abs(dEi_M - dE0i_M)/abs(dE0i_M)
         max_err = max(max_err, rel_M)

         write(uout,'(I4, 1X, ES17.8, 1X, ES17.8, 1X, ES11.2)') &
              i, dEi_M, dE0i_M, rel_M
      end associate
   end do
   call check_err(max_err, tol)

   !===================================================
   ! Test 16.
   ! Double minimum potential from
   ! B. R. Johnson, J. Chem. Phys. 67, 4086–4093 (1977)
   ! V(r) = D * (1 - exp[-B*(r - xa)])^2 + A * exp[-C*(r - xb)^2]
   ! Comparing the first 16 energies with values from the article.
   ! The current values are about 6 orders of magnitude more accurate than the article values.
   ! The equation is solved on the interval (0.5; +infinity).
   ! Two-parameter exp mapping is used here.
   ! ye is set to the value of r at which the global minimum occurs
   ! Optimal beta was obtained by hand.
   !
   !===================================================
   write(uout,*) ''
   write(uout,*) 'Test 16. Double minimum potential from B. R. Johnson, JCP. 67, 4086–4093 (1977)'
   call announce_test('Test 16. Double-minimum potential: energies')

   N  = 130 ! use N grid points
   nE = 16  ! request 16 lowest eigenvalues

   mu = 1.5403756164035d0**2/16

   call qsolver( rmin  = 0.5d0,   & ! the equation is solved on (0.5; +infinity)
          &      mu    = mu   ,   & ! effective mass
          &      N     = N,       & ! number of grid points
          &      nE    = nE,      & ! number of eigenvalues to be found
          &      V1    = V_DMin,  & ! V1 = double minimum potential
          &      mtype = 'exp',   & ! exp mapping
          &      beta  = 3.0d0,   & ! good beta value obtained by hand
          &      ye    = 1.45d0,  & ! global minimum of potential
          &      E     = Ev       & ! output array with energies
          &                )


   !===================================================
   ! Set the maximum tolerance for energies
   !===================================================
   tol = 5.0d-7
   
   !===================================================
   ! Check energies by comparing with values from the article
   ! against the calculated ones. It should be noted that the current values are about
   ! 6 orders of magnitude more accurate.
   !===================================================
   write(uout,*) ''
   write(uout,*) 'Comparing first 16 energies of double-minimum potential..'
   write(uout,*) ''
   write(uout,'(A4, A21, A17, A15)') '', 'Calc. energy, cm-1', 'Johnson cm-1', 'Relative diff.'

   
   !===================================================
   ! Energies from Johnson's article
   !===================================================
   Eref(1:16) = (/  1302.500d0,  3205.307d0,  4227.339d0,  5144.251d0, &       
                    6064.241d0,  7092.679d0,  7614.622d0,  8911.545d0, &
                    9095.696d0, 10208.350d0, 10869.289d0, 11482.479d0, &
                   12353.799d0, 12972.473d0, 13690.455d0, 14435.350d0 /)

   max_err = 0.0d0
   do i = 1, nE
      associate( E0i      => Eref(i), & ! exact energy
            &    Ei       => Ev(i),   & ! calculated energy
            &    rel_err  => t1)

         rel_err = abs(Ei - E0i)/abs(E0i)
         max_err = max(max_err, rel_err)    ! max relative error
         write(uout,'(I4, F21.9, F17.3,  ES15.2)') i, Ei, E0i, rel_err
      end associate
   enddo 

   !===================================================
   ! If max_err > tol, then the test fails
   !===================================================
   call check_err(max_err, tol)

   !===================================================
   ! Test 17.
   ! Symmetric double minimum potential from
   ! B. R. Johnson, J. Chem. Phys. 67, 4086–4093 (1977)
   ! V(r) = (r*r - 1)**2
   ! mu = 100
   ! Comparing the first 16 energies with values from the article.
   ! The current values are about 5 orders of magnitude more accurate than the article values.
   ! The equation is solved on the interval (-infinity; +infinity).
   ! Double-exponential mapping is used here.
   ! ye is set to the value of r at which the global minimum occurs
   ! Optimal beta was obtained by hand.
   !
   !===================================================
   write(uout,*) ''
   write(uout,*) 'Test 17. Symmetric double minimum potential V = (r^2 - 1)^2'
   call announce_test('Test 17. Symmetric double-minimum potential: energies')

   N  = 70  ! use N grid points
   nE = 16  ! request 16 lowest eigenvalues

   mu = 100d0

   call qsolver(  mu    = mu,         & ! effective mass
          &       N     = N,          & ! number of grid points
          &       nE    = nE,         & ! number of eigenvalues to be found
          &       V1    = V_DMin_sym, & ! V1 = double minimum potential
          &       mtype = 'dexp',     & ! double exponential mapping
          &       beta  = 0.7d0,      & ! good beta value obtained by hand
          &       E     = Ev          & ! output array with energies
          &                )

   !===================================================
   ! Set the maximum tolerance for energies
   !===================================================
   tol = 1.0d-8
   
   !===================================================
   ! Check energies by comparing with values from the article
   ! against the calculated ones. It should be noted that the current values are about
   ! 5 orders of magnitude more accurate.
   !===================================================
   write(uout,*) ''
   write(uout,*) 'Comparing first 16 energies of potential..'
   write(uout,*) ''
   write(uout,'(A4, A21, A17, A15)') '', 'Calculated energy', 'Johnson', 'Relative diff.'

   
   !===================================================
   ! Energies from Johnson's article
   !===================================================
   Eref(1:16) = (/ 0.138811928d0, 0.138811949d0, 0.405026541d0, 0.405030240d0, &
                   0.650844055d0, 0.651100997d0, 0.864617277d0, 0.872446349d0, &
                   1.017228960d0, 1.078052090d0, 1.189379930d0, 1.301102700d0, &
                   1.425248200d0, 1.557185350d0, 1.696608050d0, 1.842778290d0 /)
      

   max_err = 0.0d0
   do i = 1, nE
      associate( E0i      => Eref(i), & ! exact energy
            &    Ei       => Ev(i),   & ! calculated energy
            &    rel_err  => t1)

         rel_err = abs(Ei - E0i)/abs(E0i)
         max_err = max(max_err, rel_err)    ! max relative error
         write(uout,'(I4, F21.13, F17.9,  ES15.2)') i, Ei, E0i, rel_err
      end associate
   enddo 

   !===================================================
   ! If max_err > tol, then the test fails
   !===================================================
   call check_err(max_err, tol)

   !===================================================
   ! Test 18.
   ! Lennard-Jones (6-12) potential from:
   !   V. Meshkov, A. Stolyarov, and R. J. Le Roy, Phys. Rev. A 78, 052510 (2008).
   ! V(r) = De*((re/r)^(2n) - (re/r)^n)
   ! n = 6, De = 3761, re = 1
   ! mu = 0.5
   !
   ! This potential supports 15 bound levels.
   ! The energy of the last level is examined. This level lies extremely close to the
   ! dissociation limit and is therefore difficult to compute accurately.
   !
   ! The energy is obtained as the eigenvalue closest to a reference target value E0.
   ! In this mode the ARPACK eigenvalue solver is used, which is the most accurate
   ! approach for such a near-threshold state.
   !
   ! If instead the energy is requested by index or over an energy range, the LAPACK solver
   ! is used. For the last bound level, LAPACK may return a less accurate value due to
   ! limitations of the dense-matrix diagonalization algorithm.
   !=================================================== 

   write(uout,*) ''
   write(uout,*) 'Test 18. Lennard-Jones 6-12 potential (15 bound levels).'
   call announce_test('Test 18. Lennard-Jones 6-12: last bound level')

   N  = 250  ! use N grid points
   nE = 1    ! request 1 eigenvalue

   mu = 0.5d0

   call qsolver( rmin  = 0.3d0,      & ! equation is solved on [0.3; +infinity)
          &      mu    = mu,         & ! effective mass
          &      N     = N,          & ! number of grid points
          &      nE    = nE,         & ! number of eigenvalues to be found
          &      E0   = -0.001d0,    & ! find eigenvalue closest to E0
          &      V1    = V_LJ_6_12,  & ! V1 = 15-level LJ 6-12 potential
          &      mtype = 'rat',      & ! rational mapping
          &      beta  = 0.7d0,      & ! good beta value obtained by hand
          &      ye    = 1d0,        & ! ye is set to the position of the potential minimum
          &      E     = Ev          & ! output array with energies
          &               )

   !===================================================
   ! Set the maximum tolerance for energies
   !===================================================

   tol = 1.0d-10
   
   !===================================================
   ! Checking energy by comparing with values from article
   ! against the calculated ones.
   !===================================================
   write(uout,*) ''
   write(uout,*) 'Comparing the last bound level'
   write(uout,*) ' which lies very close to dissociation limit ...'
   write(uout,*) ''
   write(uout,'(A4, A20, A20, A15)') '', 'Calc. energy cm^-1', 'Ref cm^-1', 'Relative diff.'

   
   !===================================================
   ! Reference energy from the article
   !===================================================
      
   max_err = 0.0d0
   do i = 1, nE
      associate( E0i     => -1.79465837533d-5, & ! exact energy
            &    Ei      => Ev(i),   & ! calculated energy
            &    rel_err => t1)

         rel_err = abs(Ei - E0i)/abs(E0i)
         max_err = max(max_err, rel_err)    ! max relative error
         write(uout,'(I4, ES20.11, ES20.11,  ES15.2)') i, Ei, E0i, rel_err
      end associate
   enddo 

   !===================================================
   ! If max_err > tol, then the test fails
   !===================================================
   call check_err(max_err, tol)


   !===================================================
   ! Test 19.
   ! Lennard-Jones 3-6 potential from
   ! V. Meshkov, A. Stolyarov, and R. Le Roy, Phys. Rev. A 78, 052510 (2008).
   ! V(r) = De*((re/r)^(2n) - (re/r)^n)
   ! n = 3, De = 107000, re = 1
   ! mu = 0.5
   ! The energy of the last level is examined. This level lies extremely
   ! close to the dissociation limit, which makes this an extremely difficult problem.
   ! The wavefunction of the last level does not decay until about 10^6 A.
   ! It is almost impossible to solve this problem without a mapping, which
   ! makes the wavefunction oscillations more uniform,
   ! so the 'rat' mapping seems to be best in this case.
   !
   !===================================================
 
   write(uout,*) ''
   write(uout,*) 'Test 19. Lennard-Jones 3-6 potential (301 bound levels). Extremly difficult case.'
   call announce_test('Test 19. Lennard-Jones 3-6: last bound level')

   N  = 1000 ! use N grid points
   nE = 1    ! request 1 eigenvalue

   mu = 0.5d0

   call qsolver( rmin  = 0.30d0,    & ! equation is solved on [0.3; +infinity)
          &      mu    = mu,        & ! effective mass
          &      N     = N,         & ! number of grid points
          &      nE    = nE,        & ! number of eigenvalues to be found
          &      E0   = -0.1d-10,   & ! find eigenvalue closest to E0
          &      V1    = V_LJ_3_6 , & ! V1 = 301-level LJ 3-6 potential
          &      mtype = 'rat',     & ! rational mapping
          &      beta  = 0.5d0,     & ! beta from article
          &      ye    = 1d0,       & ! ye is set to the position of the potential minimum
          &      E     = Ev         & ! output array with energies
          &               )

   !===================================================
   ! Set the maximum tolerance for energies
   !===================================================
   tol = 1.0d-10
   
   !===================================================
   ! Checking energies by comparing with value from article
   ! against the calculated ones.
   !===================================================
   write(uout,*) ''
   write(uout,*) 'Comparing the last bound level'
   write(uout,*) ' which lies extremely close to dissociation limit ...'
   write(uout,*) ''
   write(uout,'(A4, A20, A20, A15)') '', 'Calc. energy cm^-1', 'Ref cm^-1', 'Relative diff.'

   
   max_err = 0.0d0
   do i = 1, nE
      associate( E0i      => -1.5738671521d-13, & ! exact energy
            &    Ei       => Ev(i),             & ! calculated energy
            &    rel_err  => t1)

         rel_err = abs(Ei - E0i)/abs(E0i)
         max_err = max(max_err, rel_err)    ! max relative error
         write(uout,'(I4, ES20.10, ES20.10,  ES15.2)') i, Ei, E0i, rel_err
      end associate
   enddo 

   !===================================================
   ! If max_err > tol, then the test fails
   !===================================================
   call check_err(max_err, tol)

   !===================================================
   ! Test 20. Resonant states of the 7.5*r^2*exp(-r) potential from
   ! R. A. Bain, J. N. Burdsley, B. R. Junker and C. V. Sukumar
   ! J. Phys. B 7(16), 2189 (1974)
   ! mu = 1.0, rmin = 0, rmax = infinity
   ! This potential does not support bound levels. It has 7 resonance eigenstates
   ! with complex energies Ec_n = E_n - i*Gamma_n/2, with E_n > 0.
   ! Complex scaling is used to transform the boundary condition at infinity into a zero boundary condition.
   ! To enable complex scaling, the scaling angle (phi) must be provided. In this test
   ! phi = 1.0 rad. This relatively large value allows all 7 levels to be computed at once.
   ! In the output table there is a threshold angle, defined as theta_n =
   ! arctan(Gamma_n/(2E_n)), which estimates the minimal complex-rotation angle (phi)
   ! required to uncover the state.
   ! A spurious-eigenvalue filter is enabled automatically in this case.
   ! The program uses the complex potential 7.5*r^2*exp(-r), V_BBJS_cmplx.
   ! 
   ! The calculated energies and widths are compared with reference values that were
   ! carefully computed with qsolver; all significant digits are stable with respect to N
   ! and the model parameters.
   ! 
   ! The 'dexp' and 'tanh' mappings seem to be best here; they allow a minimal number
   ! of grid points to be used.
   !===================================================
 
   write(uout,*) ''
   write(uout,*) 'Test 20. Quasibound states of V = 7.5*r^2*exp(-r) potential'
   call announce_test('Test 20. 7.5*r^2*exp(-r): all resonances, complex potential')
   write(uout,*) ' Complex potential provided'

   N  =  80 ! use N grid points

   mu = 1.0d0

   call qsolver( rmin  = 0.0d0,        & ! equation is solved on [0.0; +infinity)
          &      mu    = mu,           & ! effective mass
          &      N     = N,            & ! number of grid points
          &      nEf   = nEf,          & ! number of eigenvalues found
          &      V1c   = V_BBJS_cmplx, & ! complex potential
          &      phi   = 1.0d0,        & ! complex-scaling angle
          &      sortW = 1,            & ! sort resonances by width (ascending)
          &      mtype = 'tanh',       & ! one-parameter tanh mapping
          &      beta  = 0.15d0,       & ! mapping parameter
          &      WF1c  = WF1c,         & !
          &      E     = Ev,           & ! output array with energies
          &      W     = Wv            & ! output array with resonance width
          &               )

   !===================================================
   ! Set the maximum tolerance for energies
   !===================================================
   tol = 7.0d-11

   !===================================================
   ! Checking energies by comparing with value from article
   ! against the calculated ones.
   !===================================================
   write(uout,*) ''
   write(uout,*) ' Comparing energy of the resonant states with reference ones..'
   write(uout,*) ''

   !===================================================
   ! Reference values of energies and resonance widths
   !===================================================
   Eref(1:7) = (/ 3.4263903101482d0, 4.8348068410965d0, 5.2772798639840d0, &
             &    5.0649296073802d0, 4.2688602992626d0, 2.9477816003290d0, &
             &    1.1471837382930d0 /)
   
   Wref(1:7) = (/ 2.55489611858d-2,  2.2357533377358d0, 6.7781065905079d0, &
             &    1.1952069575731d1, 1.7433816867896d1, 2.3061029462520d1, &
             &    2.8738014274194d1 /)

   write(uout,'(A4, A20, A12, A20, A12, A7)') '', 'Calc. energies au', 'Rel. diff.', & 
             &                         'Calc. widths au', 'Rel. diff.', 'theta'

   max_err = 0.0d0
   do i = 1, nEf
      associate( E0i       => Eref(i), & ! reference energy
            &    W0i       => Wref(i), & ! reference width
            &    Ei        => Ev(i),   & ! calculated energy
            &    Gi        => Wv(i),   & ! calculated width
            &    rel_err_E => t1,      & 
            &    rel_err_G => t2,      &
            &    theta     => t3)        ! threshold angle

         rel_err_E = abs(Ei - E0i)/abs(E0i)
         rel_err_G = abs(Gi - W0i)/abs(W0i)
         theta     = atan(Gi/(2*Ei))/2/pi*180
         max_err   = max(max_err, rel_err_E, rel_err_G)  ! max relative error
         write(uout,'(I4, F20.12, ES12.2, F20.12, ES12.2, F7.1)') i, Ei, rel_err_E, Gi, rel_err_G, theta  
      end associate
   enddo 
   
   if (nEf < 7) then
      max_err = 2*tol
      write(uout,*) 'TEST FAILED: expected 7 resonances, found ', nEf 
   endif 

   !===================================================
   ! If max_err > tol, then the test fails
   !===================================================
   call check_err(max_err, tol)

   !===================================================
   ! Test 21. The same test as Test 20, but with a real potential provided.
   ! Resonant states of the 7.5*r^2*exp(-r) potential from
   ! R. A. Bain, J. N. Burdsley, B. R. Junker and C. V. Sukumar
   ! J. Phys. B 7(16), 2189 (1974)
   ! mu = 1.0, rmin = 0, rmax = infinity
   ! This potential does not support bound levels.
   ! It has 7 resonance eigenstates with complex energies Ec_n = E_n - i*Gamma_n/2,
   ! with E_n > 0.
   ! Complex scaling is used to transform the boundary condition at infinity into a zero boundary condition.
   ! To enable complex scaling, the scaling angle (phi) must be provided. In this test
   ! phi = 1.0 rad. This relatively large value allows all 7 levels to be computed at once.
   ! In the output table there is a threshold angle, defined as theta_n =
   ! arctan(Gamma_n/(2E_n)), which estimates the minimal complex-rotation angle (phi)
   ! required to uncover the state.
   !
   ! A spurious-eigenvalue filter is enabled automatically in this case.
   !
   ! The program uses the real potential 7.5*r^2*exp(-r), V_BBJS_real, which is
   ! automatically extrapolated to the complex plane.
   ! 
   ! The calculated energies and widths are compared with reference values that were
   ! carefully computed with qsolver; all significant digits are stable with respect to N
   ! and the model parameters.
   ! 
   ! Extrapolation into the complex plane gives slightly less accurate energies and widths.
   ! The situation is complicated by the fact that the angle phi is relatively large, which means
   ! that the calculation probes relatively far into the complex plane.
   ! 
   ! If a real potential is used to calculate resonant states, automatic extrapolation
   ! to the complex plane is applied. In this case the 'atan' mapping typically should be used.
   ! Otherwise the potential V(r) in the mapped coordinate x may have a singularity at x = 1,
   ! corresponding to r = infinity, which makes it difficult to approximate with a polynomial
   ! and to extrapolate to the complex plane.
   !===================================================
 
   write(uout,*) ''
   write(uout,*) 'Test 21. Quasibound states of V = 7.5*r^2*exp(-r) potential'
   call announce_test('Test 21. 7.5*r^2*exp(-r): all resonances, real potential')
   write(uout,*) ' Real potential provided and automatically extrapolated to the complex plane'

   N  =  130 ! use N grid points

   mu = 1.0d0

   call qsolver( rmin  = 0.0d0,        & ! equation is solved on [0.0; +infinity)
          &      mu    = mu,           & ! effective mass
          &      N     = N,            & ! number of grid points
          &      nEf   = nEf,          & ! number of eigenvalues found
          &      V1    = V_BBJS_real,  & ! real potential
          &      phi   = 1.0d0,        & ! complex-scaling angle
          &      sortW = 1,            & ! sort resonances by width (ascending)
          &      mtype = 'atan',       & ! one-parameter atan mapping
          &      beta  = 0.07d0,       & ! mapping parameter
          &      E     = Ev,           & ! output array with energies
          &      W     = Wv            & ! output array with resonance width
          &               )

   !===================================================
   ! Set the maximum tolerance for energies. 
   !===================================================
   tol = 2.0d-7
   
   !===================================================
   ! Checking energies by comparing with value from article
   ! against the calculated ones.
   !===================================================
   write(uout,*) ''
   write(uout,*) ' Comparing energy of the resonant states with reference ones..'
   write(uout,*) ''
   
   write(uout,'(A4, A20, A12, A20, A12, A7)') '', 'Calc. energies au', 'Rel. diff.', & 
             &                         'Calc. widths au', 'Rel. diff.', 'theta'
   max_err = 0.0d0
   do i = 1, nEf
      associate( E0i       => Eref(i), & ! reference energy
            &    W0i       => Wref(i), & ! reference width
            &    Ei        => Ev(i),   & ! calculated energy
            &    Gi        => Wv(i),   & ! calculated width
            &    rel_err_E => t1,      & 
            &    rel_err_G => t2,      &
            &    theta     => t3)        ! threshold angle

         rel_err_E = abs(Ei - E0i)/abs(E0i)
         rel_err_G = abs(Gi - W0i)/abs(W0i)
         theta     = atan(Gi/(2*Ei))/2/pi*180
         max_err   = max(max_err, rel_err_E, rel_err_G)  ! max relative error
         write(uout,'(I4, F20.12, ES12.2, F20.12, ES12.2, F7.1)') i, Ei, rel_err_E, Gi, rel_err_G, theta  
      end associate
   enddo 

   if (nEf < 7) then
      max_err = 2*tol
      write(uout,*) 'TEST FAILED: expected 7 resonances, found ', nEf 
   endif 
   !===================================================
   ! If max_err > tol, then the test fails
   !===================================================
   call check_err(max_err, tol)

   !===================================================
   ! Test 22. Calculation of a selected resonant state of the potential
   ! V(r)=7.5*r^2*exp(-r) from
   ! R. A. Bain, J. N. Burdsley, B. R. Junker and C. V. Sukumar,
   ! J. Phys. B 7(16), 2189 (1974).
   ! mu = 1.0, rmin = 0, rmax = infinity.
   !
   ! The same problem as in Tests 20-21, but now a selected resonance is searched for
   ! by providing a target energy E0 = 3.4.
   ! The program uses the complex-rotated potential 7.5*r^2*exp(-r) (V_BBJS_cmplx).
   !
   ! Complex scaling is used to convert the outgoing-wave boundary condition at infinity
   ! into a decaying (L^2) boundary condition. To enable complex scaling, the scaling angle
   ! (phi) must be provided.
   !
   ! To transform the equation to a finite domain, the 'tanh' mapping is used.
   ! In this test the complex-scaling angle is set to phi = 0.7 rad, and the mapping
   ! parameter is beta = 0.3, which appears to be close to optimal for this case.
   ! These values allow full precision to be achieved with about 30 grid points.
   ! If the parameters are far from optimal, many more points may be needed.
   !
   ! The 'dexp' and 'tanh' mappings seem to work best here; they allow a minimal number
   ! of grid points to be used. For the 'tanh' mapping, the optimal parameter beta is
   ! approximately 0.3.
   !
   ! In the output table a threshold angle is printed, defined as
   !   theta_n = arctan(Gamma_n/(2*E_n)),
   ! which estimates the minimal complex-rotation angle (phi) required to uncover the state.
   !
   ! A spurious-eigenvalue filter is enabled automatically in this case.
   !
   ! The calculated energy and width are compared with reference values that were
   ! carefully computed with qsolver; all significant digits are stable with respect to N
   ! and the model parameters.
   !===================================================
 
   write(uout,*) ''
   write(uout,*) 'Test 22. Quasibound states of V = 7.5*r^2*exp(-r) potential'
   call announce_test('Test 22. 7.5*r^2*exp(-r): target resonance, complex potential')

   write(uout,*) ' Find resonance level with energy closest to E0 = 3.4 au'
   write(uout,*) ' Complex potential provided'

   N  =  30 ! use N grid points

   mu = 1.0d0

   call qsolver( rmin  = 0.0d0,        & ! equation is solved on [0.0; +infinity)
          &      mu    = mu,           & ! effective mass
          &      N     = N,            & ! number of grid points
          &      E0    = 3.4d0,        & ! find the eigenvalue closest to E0 = 3.4 and W0 = 0.0
          &      nE    = 1,            & ! request 1 eigenvalue
          &      nEf   = nEf,          & ! number of eigenvalues found
          &      V1c   = V_BBJS_cmplx, & ! complex potential
          &      phi   = 0.7d0,        & ! complex-scaling angle
          &      mtype = 'tanh',       & ! one-parameter tanh mapping
          &      beta  = 0.3d0,        & ! mapping parameter
          &      E     = Ev,           & ! output array with energies
          &      W     = Wv            & ! output array with resonance width
          &               )

   !===================================================
   ! Set the maximum tolerance for energies
   !===================================================
   tol = 1.0d-11
   
   !===================================================
   ! Checking energies by comparing with value from article
   ! against the calculated ones.
   !===================================================
   write(uout,*) ''
   write(uout,*) ' Comparing energy of the resonant state with reference one..'
   write(uout,*) ''

   write(uout,'(A4, A20, A12, A20, A12, A7)') '', 'Calc. energies au', 'Rel. diff.', & 
             &                         'Calc. widths au', 'Rel. diff.', 'theta'
   max_err = 0.0d0
   do i = 1, nEf
      associate( E0i       => Eref(i), & ! reference energy
            &    W0i       => Wref(i), & ! reference width
            &    Ei        => Ev(i),   & ! calculated energy
            &    Gi        => Wv(i),   & ! calculated width
            &    rel_err_E => t1,      & 
            &    rel_err_G => t2,      &
            &    theta     => t3)        ! threshold angle

         rel_err_E = abs(Ei - E0i)/abs(E0i)
         rel_err_G = abs(Gi - W0i)/abs(W0i)
         theta     = atan(Gi/(2*Ei))/2/pi*180
         max_err   = max(max_err, rel_err_E, rel_err_G)  ! max relative error
         write(uout,'(I4, F20.13, ES12.2, F20.13, ES12.2, F7.1)') i, Ei, rel_err_E, Gi, rel_err_G, theta  
      end associate
   enddo 

   if (nEf < 1) then
      max_err = 2*tol
      write(uout,*) 'TEST FAILED: no resonance found'
   endif 

   !===================================================
   ! If max_err > tol, then the test fails
   !===================================================
   call check_err(max_err, tol)

   !===================================================
   ! Test 23. Calculation of a selected resonant state of the potential
   ! V(r)=7.5*r^2*exp(-r) from
   ! R. A. Bain, J. N. Burdsley, B. R. Junker and C. V. Sukumar,
   ! J. Phys. B 7(16), 2189 (1974).
   ! mu = 1.0, rmin = 0, rmax = infinity.
   !
   ! The same problem as in Tests 20-21, but now a selected resonance is searched for
   ! by providing a target energy E0 = 3.4.
   !
   ! Complex scaling is used to convert the outgoing-wave boundary condition at infinity
   ! into a decaying (L^2) boundary condition. To enable complex scaling, the scaling angle
   ! (phi) must be provided.
   !
   ! The program uses the real, user-provided potential 7.5*r^2*exp(-r) (V_BBJS_real) and
   ! automatically extrapolates it into the complex plane.
   !
   ! To transform the equation to a finite domain, the 'atan' mapping is used.
   ! This is the recommended choice when a real potential is provided and resonant states
   ! are computed by complex scaling.
   !
   ! In this test the complex-scaling angle is set to phi = 0.4 rad, and the mapping
   ! parameter is beta = 0.1, which appears to be close to optimal for this case.
   ! These values allow full precision to be achieved with about 55 grid points.
   ! If the parameters are far from optimal, many more points may be needed.
   !
   ! In the output table a threshold angle is printed, defined as
   !   theta_n = arctan(Gamma_n/(2*E_n)),
   ! which estimates the minimal complex-rotation angle (phi) required to uncover the state.
   !
   ! A spurious-eigenvalue filter is enabled automatically in this case.
   !
   ! The calculated energy and width are compared with reference values that were
   ! carefully computed with qsolver; all significant digits are stable with respect to N
   ! and the model parameters.
   !===================================================

   write(uout,*) ''
   write(uout,*) 'Test 23. Quasibound states of V = 7.5*r^2*exp(-r) potential'
   call announce_test('Test 23. 7.5*r^2*exp(-r): target resonance, real potential')
   write(uout,*) ' Find resonance level with energy closest to E0 = 3.4 au'
   write(uout,*) ' Real potential provided'

   N  =  55 ! use N grid points

   mu = 1.0d0

   call qsolver( rmin  = 0.0d0,        & ! equation is solved on [0.0; +infinity)
          &      mu    = mu,           & ! effective mass
          &      N     = N,            & ! number of grid points
          &      E0    = 3.4d0,        & ! find the eigenvalue closest to E0 = 3.4 and W0 = 0.0
          &      nE    = 1,            & ! request 1 eigenvalue
          &      nEf   = nEf,          & ! number of eigenvalues found
          &      V1    = V_BBJS_real,  & ! real potential
          &      phi   = 0.4d0,        & ! complex-scaling angle
          &      mtype = 'atan',       & ! one-parameter atan mapping
          &      beta  = 0.1d0,        & ! mapping parameter
          &      E     = Ev,           & ! output array with energies
          &      W     = Wv            & ! output array with resonance width
          &               )

   !===================================================
   ! Set the maximum tolerance for energies
   !===================================================
   tol = 1.0d-11
   
   !===================================================
   ! Checking energies by comparing with value from article
   ! against the calculated ones.
   !===================================================
   write(uout,*) ''
   write(uout,*) ' Comparing energy of the resonant state with reference one..'
   write(uout,*) ''

   write(uout,'(A4, A20, A12, A20, A12, A7)') '', 'Calc. energies au', 'Rel. diff.', & 
             &                         'Calc. widths au', 'Rel. diff.', 'theta'

   max_err = 0.0d0
   
   do i = 1, nEf
      associate( E0i       => Eref(i), & ! reference energy
            &    W0i       => Wref(i), & ! reference width
            &    Ei        => Ev(i),   & ! calculated energy
            &    Gi        => Wv(i),   & ! calculated width
            &    rel_err_E => t1,      & 
            &    rel_err_G => t2,      &
            &    theta     => t3)        ! threshold angle

         rel_err_E = abs(Ei - E0i)/abs(E0i)
         rel_err_G = abs(Gi - W0i)/abs(W0i)
         theta     = atan(Gi/(2*Ei))/2/pi*180
         max_err   = max(max_err, rel_err_E, rel_err_G)  ! max relative error
         write(uout,'(I4, F20.13, ES12.2, F20.13, ES12.2, F7.1)') i, Ei, rel_err_E, Gi, rel_err_G, theta  
      end associate
   enddo 

   if (nEf < 1) then
      max_err = 2*tol
      write(uout,*) 'TEST FAILED: no resonance found'
   endif 

   !===================================================
   ! If max_err > tol, then the test fails
   !===================================================
   call check_err(max_err, tol)

   !===================================================
   ! Test 24. Li2 B^1Pi_u state: resonance levels.
   ! Potential and data from:
   !   Y. Huang and R. J. Le Roy, J. Chem. Phys. 119, 7398 (2003).
   !
   ! The potential has a rotationless barrier that extends above the dissociation
   ! asymptote.
   ! 
   ! In Test 24 the complex-rotated potential is used directly, whereas in Test 25 the
   ! real potential is provided and extrapolated into the complex plane. The results are
   ! almost the same.
   !
   ! Exterior complex scaling is used in this case because the wavefunction becomes
   ! too intricate in the complex-r region near equilibrium. This inner region should
   ! be avoided; otherwise many grid points may be required if complex scaling starts
   ! too early.
   !
   ! Here the exterior complex-scaling parameters are phi = 0.6 rad and rc = 5.0 A,
   ! so the complex rotation starts in the beginning of the long-range region.
   !
   ! In this test a complex potential is used. Using the complex potential 
   ! gives absolute errors in the resonance energy and width below 1.0e-10 cm^-1, 
   ! whereas using the real potential gives errors below 1.0e-8 cm^-1.
   !
   ! A two-parameter 'atan' mapping is used with ye = 2.9 (equilibrium distance) and
   ! beta = 1.0.
   !
   ! The calculated resonance energies and widths are compared with experimental and
   ! literature values, as well as with reference high-precision values computed by qsolver.
   !===================================================

   write(uout,*) ''
   write(uout,*) 'Test 24. B^1Pi_u state of lithium dimer (Li2). Resonant states calculation.'
   call announce_test('Test 24. Li2 B state: complex potential')
   write(uout,*) '         Complex potential is provided directly in the complex plane.'
   write(uout,*) '          Comparison with experiment, literature and reference values.' 
   write(uout,*) '          All values are in cm-1.' 

   !===================================================
   ! Here are the data for Tests 24 and 25.
   !===================================================

   !===================================================
   ! Rotational quantum number for each level
   !===================================================
   Jv(1:52) = (/ &
            &    1, 2, 3, 4, 5, 6, &
            &    7, 8, 9, 10, 11, 11, &
            &    12, 12, 13, 14, 15, 16, &
            &    17, 18, 19, 20, 21, 22, &
            &    23, 24, 25, 25, 26, 26, &
            &    27, 27, 28, 29, 30, 31, &
            &    32, 33, 35, 36, 37, 38, &
            &    39, 41, 43, 44, 45, 46, &
            &    49, 54, 56, 62 /)

   !===================================================
   ! Experimental energies
   !===================================================
   Eexp(1:52) = (/ &
              &    15382.107d0, 15383.136d0, 15384.675d0, 15386.719d0, 15389.261d0, 15392.293d0, &
              &    15395.807d0, 15399.787d0, 15404.223d0, 15328.254d0, 15335.007d0, 15414.423d0, &
              &    15342.335d0, 15420.163d0, 15350.225d0, 15358.664d0, 15367.637d0, 15377.126d0, &
              &    15387.112d0, 15397.57d0,  15408.473d0, 15419.786d0, 15431.47d0,  15443.477d0, &
              &    15455.757d0, 15468.262d0, 15396.036d0, 15480.961d0, 15412.235d0, 15493.856d0, &
              &    15428.859d0, 15506.834d0, 15445.871d0, 15463.229d0, 15480.878d0, 15498.752d0, &
              &    15516.766d0, 15534.822d0, 15485.197d0, 15507.709d0, 15530.467d0, 15553.403d0, &
              &    15576.423d0, 15526.33d0,  15580.997d0, 15608.588d0, 15636.24d0,  15663.833d0, &
              &    15652.779d0, 15710.885d0, 15783.338d0, 15904.75d0 /)

   !===================================================
   ! Experimental resonance widths
   !===================================================
   Wexp(1:52) = (/ &
              &    0.8d0, 0.91d0, 0.85d0, 1.1d0, 1.1d0, 1.2d0, &
              &    1.4d0, 1.9d0, 1.9d0, 0.000122d0, 0.000196d0, 3.6d0, &
              &    0.00043d0, 4.6d0, 0.00078d0, 0.00158d0, 0.00472d0, 0.0063d0, &
              &    0.00892d0, 0.0238d0, 0.0413d0, 0.0741d0, 0.18d0, 0.33d0, &
              &    0.7d0, 1.3d0, 0.00011d0, 2.6d0, 0.000414d0, 4.0d0, &
              &    0.00111d0, 8.0d0, 0.00351d0, 0.0115d0, 0.0364d0, 0.108d0, &
              &    0.51d0, 0.68d0, 0.000345d0, 0.00127d0, 0.00557d0, 0.0231d0, &
              &    0.079d0, 0.000034d0, 0.0005d0, 0.00315d0, 0.0194d0, 0.1d0, &
              &    0.000117d0, 0.000095d0, 0.0017d0, 0.009423d0 /)

   !===================================================
   ! Resonance widths reported by Y. Huang and R. J. Le Roy
   !===================================================
   Wlit(1:52) = (/ &
              &    0.692d0, 0.732d0, 0.794d0, 0.885d0, 1.01d0, 1.18d0, &
              &    1.4d0, 1.7d0, 2.09d0, 0.000137d0, 0.000234d0, 3.27d0, &
              &    0.000415d0, 4.14d0, 0.000757d0, 0.00142d0, 0.00274d0, 0.00537d0, &
              &    0.0107d0, 0.0216d0, 0.0436d0, 0.0881d0, 0.176d0, 0.347d0, &
              &    0.663d0, 1.22d0, 0.000119d0, 2.16d0, 0.000389d0, 3.65d0, &
              &    0.00125d0, 5.84d0, 0.00393d0, 0.012d0, 0.0356d0, 0.1d0, &
              &    0.267d0, 0.661d0, 0.000331d0, 0.00144d0, 0.00592d0, 0.0228d0, &
              &    0.0814d0, 0.0000266d0, 0.000871d0, 0.00441d0, 0.0204d0, 0.0855d0, &
              &    0.00038d0, 0.0000615d0, 0.00309d0, 0.00881d0 /)

   !===================================================
   ! Reference values of energies and widths obtained from
   ! a previous run of this same
   ! Test 24 with N = 250.
   !
   ! Eref is rounded to 10 digits after the decimal point.
   ! Wref is rounded to 11 digits after the decimal point.
   ! These values are stable to at least this precision.
   !===================================================
   Eref(1:52) = (/ &
           &    15382.1210546216d0, 15383.1508238583d0, 15384.6911535300d0, 15386.7368448813d0, &
           &    15389.2809665541d0, 15392.3148886103d0, 15395.8283675982d0, 15399.8097087374d0, &
           &    15404.2460296346d0, 15328.2479079407d0, 15335.0004900976d0, 15414.4285372244d0, &
           &    15342.3270022788d0, 15420.1469815977d0, 15350.2158260565d0, 15358.6537233534d0, &
           &    15367.6255595670d0, 15377.1139265261d0, 15387.0986273308d0, 15397.5559830402d0, &
           &    15408.4579452831d0, 15419.7710921126d0, 15431.4558161447d0, 15443.4664292671d0, &
           &    15455.7533433612d0, 15468.2683147723d0, 15396.0048350840d0, 15480.9721685919d0, &
           &    15412.2029689462d0, 15493.8420070228d0, 15428.8251488065d0, 15506.8743654129d0, &
           &    15445.8358204646d0, 15463.1921290686d0, 15480.8409840050d0, 15498.7155184633d0, &
           &    15516.7330769100d0, 15534.8002539405d0, 15485.1405587004d0, 15507.6502961271d0, &
           &    15530.4078517438d0, 15553.3432942357d0, 15576.3648390144d0, 15526.2502802080d0, &
           &    15580.9150992062d0, 15608.5050388847d0, 15636.1570743689d0, 15663.7530366572d0, &
           &    15652.6720632146d0, 15710.7538396856d0, 15783.2074401105d0, 15904.5975053828d0 /)

   Wref(1:52) = (/ &
           &    0.69362905981d0, 0.73322407497d0, 0.79611455197d0, 0.88681621634d0, &
           &    1.01188775816d0, 1.18039976475d0, 1.40451538983d0, 1.70015560872d0, &
           &    2.08769699611d0, 0.00013705569d0, 0.00023433687d0, 3.24604153469d0, &
           &    0.00041483509d0, 4.08496147603d0, 0.00075777071d0, 0.00142321003d0, &
           &    0.00273782925d0, 0.00537235356d0, 0.01070541774d0, 0.02155647848d0, &
           &    0.04361761523d0, 0.08811615719d0, 0.17639521539d0, 0.34687579819d0, &
           &    0.66367242597d0, 1.22373771976d0, 0.00011893410d0, 2.15743419826d0, &
           &    0.00038831910d0, 3.61941106636d0, 0.00124784524d0, 5.77316934883d0, &
           &    0.00392561057d0, 0.01201506136d0, 0.03550319770d0, 0.10027073545d0, &
           &    0.26715095582d0, 0.66093536765d0, 0.00032981102d0, 0.00143541180d0, &
           &    0.00590004743d0, 0.02273046123d0, 0.08116400123d0, 0.00002642858d0, &
           &    0.00086609544d0, 0.00439156175d0, 0.02035578590d0, 0.08517988705d0, &
           &    0.00037719843d0, 0.00006101623d0, 0.00306876132d0, 0.00874968170d0 /)

   write(uout,*) ''
   write(uout,'(A4, A4, A11, A18, A10, A10, A15, 3A10)') '', 'J', 'Eexp', 'Ecalc',  &
            &     'Wexp', 'Wlit', 'Wcalc', 'E - Eref', 'W - Wref'

   N  = 250

   ! Mass for 7Li7Li isotopomer
   mu = 7.0160034366d0/(2*5.485799090441d-4*2.1947463136314d5*0.529177210544d0**2)

   ! Tolerance for comparison with the reference values.
   tol = 5.0d-10
   max_err = 0.0d0

   do i = 1, 52
      Jrot = Jv(i)
      call qsolver( rmin  = 1.0d0,        & ! equation is solved on [1.0; +infinity)
             &      mu    = mu,           & ! effective mass
             &      N     = N,            & ! number of grid points
             &      E0    = Eexp(i),      & ! find the resonance closest to the experimental energy
             &      nE    = 1,            & ! request 1 eigenvalue
             &      nEf   = nEf,          & ! number of eigenvalues found
             &      V1c   = V_Li2B_cmplx, & ! complex Li2(B) potential
             &      phi   = 0.60d0,       & ! complex scaling angle
             &      rc    = 5.0d0,        & ! start of exterior complex scaling
             &      mtype = 'atan',       & ! mapping
             &      ye    = 2.9d0,        & ! two-parameter mapping center
             &      beta  = 1.0d0,        & ! mapping parameter
             &      E     = Ev,           & ! output array with resonance energies
             &      W     = Wv            & ! output array with resonance widths
             &               )

       if (nEf < 1) then
          max_err = 2*tol
          write(uout,*) 'TEST FAILED: no resonance found'
          exit
       endif 

       t1 = Ev(1) - Eref(i) ! absolute difference between calculated and reference E
       t2 = Wv(1) - Wref(i) ! absolute difference between calculated and reference W

       max_err = max(max_err, abs(t1), abs(t2))
       write(uout,'(I4, I4, F11.3, F18.10, ES10.2, ES10.2, F15.11, ES10.2, ES10.2, ES10.2)') &
               & i, Jv(i), Eexp(i), Ev(1), Wexp(i), Wlit(i), Wv(1), t1, t2

   enddo
   call check_err(max_err, tol)

   !===================================================
   ! Test 25. Li2 B^1Pi_u state: resonance levels.
   ! Potential and data from:
   !   Y. Huang and R. J. Le Roy, J. Chem. Phys. 119, 7398 (2003).
   !
   ! The potential has a rotationless barrier that extends above the dissociation
   ! asymptote.
   !
   ! Exterior complex scaling is used in this case because the wavefunction becomes
   ! too intricate in the complex-r region near equilibrium. This inner region should
   ! be avoided; otherwise many grid points may be required if complex scaling starts
   ! too early.
   !
   ! Here the exterior complex-scaling parameters are phi = 0.6 rad and rc = 5.0 A,
   ! so the complex rotation starts in the beginning of the long-range region.
   !
   ! In this test a real potential is used and is automatically extrapolated into the
   ! complex plane by qsolver. Using the complex potential gives absolute errors in the
   ! resonance energy and width below 1.0e-10 cm^-1, whereas using the real potential gives
   ! errors below 1.0e-8 cm^-1.
   !
   ! A two-parameter 'atan' mapping is used with ye = 2.9 (equilibrium distance) and
   ! beta = 1.0.
   !
   ! The calculated resonance energies and widths are compared with experimental and
   ! literature values, as well as with reference high-precision values computed by qsolver.
   !===================================================

   write(uout,*) ''
   write(uout,*) 'Test 25. B^1Pi_u state of lithium dimer (Li2). Resonant states calculation.'
   call announce_test('Test 25. Li2 B state: real potential')
   write(uout,*) '     Real potential is automatically extrapolated to the complex plane.'
   write(uout,*) '     Comparison with experiment, literature and reference values.' 
   write(uout,*) '     All values are in cm-1.' 


   write(uout,*) ''
   write(uout,'(A4, A4, A11, A18, A10, A10, A15, 3A10)') '', 'J', 'Eexp', 'Ecalc',  &
            &     'Wexp', 'Wlit', 'Wcalc', 'E - Eref', 'W - Wref'

   N  = 200

   ! Tolerance for comparison with the reference values.
   tol = 5.0d-8
   max_err = 0.0d0

   do i = 1, 52
      Jrot = Jv(i)
      call qsolver( rmin  = 1.0d0,        & ! equation is solved on [1.0; +infinity)
             &      mu    = mu,           & ! effective mass
             &      N     = N,            & ! number of grid points
             &      E0    = Eexp(i),      & ! find the resonance closest to the experimental energy
             &      nE    = 1,            & ! request 1 eigenvalue
             &      nEf   = nEf,          & ! number of eigenvalues found
             &      V1    = V_Li2B_real,  & ! real Li2(B) potential
             &      phi   = 0.60d0,       & ! complex scaling angle
             &      rc    = 5.0d0,        & ! start of exterior complex scaling
             &      mtype = 'atan',       & ! mapping
             &      ye    = 2.9d0,        & ! two-parameter mapping center
             &      beta  = 1.0d0,        & ! mapping parameter
             &      E     = Ev,           & ! output array with resonance energies
             &      W     = Wv            & ! output array with resonance widths
             &               )

       if (nEf < 1) then
          max_err = 2*tol
          write(uout,*) 'TEST FAILED: no resonance found'
          exit
       endif 

      t1 = Ev(1) - Eref(i) ! absolute difference between calculated and reference E
      t2 = Wv(1) - Wref(i) ! absolute difference between calculated and reference W

      max_err = max(max_err, abs(t1), abs(t2))
      write(uout,'(I4, I4, F11.3, F18.10, ES10.2, ES10.2, F15.11, ES10.2, ES10.2, ES10.2)') &
               & i, Jv(i), Eexp(i), Ev(1), Wexp(i), Wlit(i), Wv(1), t1, t2
   enddo
   call check_err(max_err, tol)

   !=========================================================
   ! Test 26. Radial equation for the Coulomb + inverse-square
   ! potential from:
   !   G. Doolen, Int. J. Quantum Chem. 14, 523-529 (1976)
   !
   !   -1/2 psi''(r) + (1/r - gamma/r**2) * psi(r) = E * psi(r)
   !
   ! This Hamiltonian roughly corresponds to the model of
   ! the interaction of hydrogen anion (H-) with an electron.
   !
   ! For gamma > 1/8 the problem has resonant states with complex
   ! eigenvalues Ec_n = E_n - i*G_n/2, 
   ! where
   !   E_n = 2 * (8*gamma - 4*n**2 - 4*n - 2) / (4*n**2 + 4*n + 8*gamma)**2,
   !   G_n = 4 * (4*n + 2)*sqrt(8*gamma - 1) / (4*n**2 + 4*n + 8*gamma)**2.
   !
   ! Near the origin the solutions behave as
   !   psi(r) ~ r**((1 + i*sqrt(8*gamma - 1))/2)
   ! so the wavefunction oscillates infinitely rapidly as r -> 0 and is
   ! difficult to approximate numerically.
   !
   ! To remove this singularity we transform the wavefunction as
   !   psi(r) = r**lambda * phi(r),
   ! where lambda satisfies
   !   lambda * (lambda - 1) = -2*gamma,
   !   lambda = (1 - i*sqrt(8*gamma - 1))/2.
   !
   ! This transformation only removes the singular branch. It does not
   ! force the transformed wavefunction to vanish at r = 0; instead,
   ! phi(r) is finite there. Therefore the left boundary condition is
   ! regular/finite rather than Dirichlet.
   !
   ! The transformed equation has the qsolver form
   ! (mu = 1, is_B_sym = .false.):
   !
   !   -1/2 phi''(r) + 1/2 B(r)*phi'(r) + V(r)*phi(r) = E*phi(r)
   !
   ! with
   !   V(r) = 1/r,
   !   B(r) = -2*lambda/r.
   !
   ! The functions V(r) and B(r) are therefore complex.
   ! The left regular boundary is handled by a Radau grid without the
   ! left endpoint r = 0.
   !
   ! Parameter choice: gamma = 401/8 = 50.125,
   !   lambda = 1/2 - 10*i.
   ! This value of gamma supports exactly 10 resonant states with E_n > 0.
   ! In this case
   !
   !   E_n = 2 * (399 - 4*n*(n + 1)) / (4*n**2 + 4*n + 401)**2,
   !   G_n = 160*(2*n + 1) / (4*n**2 + 4*n + 401)**2.
   !
   ! To reduce the semi-infinite interval to a finite domain the double
   ! exponential mapping 'dexp' is used.
   !
   ! The relative error of E and G is around 1d-8 to 1d-9, which is somewhat 
   ! larger than in the other tests, but the absolute error is still
   ! around 1d-10 atomic units. The optimal mapping parameter is very small
   ! (beta = 4d-4) because the resonance wavefunctions extend to several
   ! thousand bohr before they start to damp.
   !
   ! Also, to isolate these levels it is necessary to rotate rather far into
   ! the complex plane (phi = 1.1 rad).
   !
   ! Since the problem is non-Hermitian and a target value E0 is not selected,
   ! qsolver finds all eigenvalues of the Hamiltonian matrix. The eigenvalue filter
   ! leaves only the physical ones. The parameter Emin = 0.0 keeps only eigenvalues
   ! with positive energies. The total number of eigenvalues should be 10.
   !
   !=========================================================

   write(uout,*) ''
   write(uout,*) 'Test 26. Doolen transformed Coulomb plus inverse-square potential'
   call announce_test('Test 26. Doolen model: resonances')

   N  = 150

   tol = 5.0d-6

   write(uout,*) ''
   write(uout,*) 'Finding all 10 resonant levels with E > 0... '
   write(uout,*) ''
   write(uout,'(A4, 4A15, 2A10)') 'n', 'Calc. E [au]', 'Exact E [au]', &
                              &     'Calc. G [au]', 'Exact G [au]',  &
                              &     'Rel_err E', 'Rel_err G'

   max_err = 0.0d0

   call qsolver( rmin     = 0.0d0,           & ! solved on (0; +infinity)
          &      bc_left  = BC_REGULAR,      & ! regular finite behavior at r = 0
          &      mu       = 1.0d0,           & ! atomic units
          &      N        = N,               & ! number of collocation points
          &      nEf      = nEf,             & ! number of eigenvalues found
          &      Emin     = 0.0d0,           & ! find only positive energies 
          &      V1c      = V_Doolen_cmplx,  & ! complex potential
          &      B1c      = B_Doolen_cmplx,  & ! complex B-coupling
          &      is_B_sym = .false.,         & ! B enters as plain B*d/dr
          &      mtype    = 'dexp',          & ! double exponential mapping for semi-infinite domain
          &      beta     = 4.0d-4,          & ! mapping parameter
          &      phi      = 1.1d0,           & ! complex scaling angle
          &      sortE    = 2,               & ! sort by Re(E), descending
          &      E        = Ev,              & ! output resonance energies
          &      W     = Wv,              & ! output resonance widths
          &      rvc      = rvc,             &
          &      qwc      = qwc,             &
          &      WF1c     = WF1c             &
          &               )

   do i = 1, nEf
      associate( E0   => 2*(399 - 4.0d0*(i-1)*i)/  & ! exact energy
            &            (4*(i-1)*i + 401.0d0)**2, & 
            &    W0   => 160.0d0*(2*i - 1)/        & ! exact width
            &            (4*(i-1)*i + 401.0d0)**2, &
            &    errE => t1,                       &
            &    errG => t2 )
   
         errE = abs(Ev(i) - E0)/abs(E0)
         errG = abs(Wv(i) - W0)/abs(W0)
   
         max_err = max(max_err, errE, errG)
         write(uout,'(I4, 4ES15.7, 2ES10.2)') i-1, Ev(i), E0, &
            &                                   Wv(i), W0, errE, errG
      end associate
   enddo
   call check_err(max_err, tol)

   !=========================================================
   ! Test 27. Compare calculated and analytical wavefunctions
   ! for Doolen transformed Coulomb + inverse-square potential.
   !
   ! Exact transformed wavefunction used for comparison:
   !   phi_n(r) = exp(kappa_n*r)*1F1(-n, 2*lambda, -2*kappa_n*r),
   !   kappa_n = 1/(n + lambda),
   !   lambda  = (1 - i*sqrt(8*gamma - 1))/2.
   !
   ! Here 1F1(a,b,z) is Kummer's confluent hypergeometric function.
   ! In our case a = -n (n = 0,1,2,...), so the series terminates and
   ! 1F1(-n,b,z) reduces to a finite polynomial:
   !   1F1(-n,b,z) = sum_{k=0}^{n} [(-n)_k/(b)_k] * z^k/k!
   !
   ! Physical wavefunction:
   !   psi_n(r) = r**lambda * phi_n(r).
   !
   ! The function WF_Doolen returns the transformed wavefunction phi_n(r),
   ! normalized by analytical continuation of the c-norm:
   !
   !   integral_0^infinity phi_n(r)**2 dr = 1.
   !
   ! If
   !   1F1(-n, 2*lambda, -2*kappa*r) = sum_{k=0}^{n} c_k*r^k,
   ! where
   !   c_k = (-n)_k*(-2*kappa)^k / ((2*lambda)_k*k!),
   ! then the analytical normalization factor is
   !
   !   N_n^2 = sum_{k=0}^{n} sum_{j=0}^{n}
   !           c_k*c_j*Gamma(k+j+1)/(-2*kappa)^(k+j+1),
   !
   ! and WF_Doolen returns phi_n(r)/sqrt(N_n^2).
   !=========================================================

   tol = 1.0d-4

   write(uout,*) ''
   write(uout,*) 'Test 27. Doolen wavefunction comparison'
   call announce_test('Test 27. Doolen model: wavefunctions')
   write(uout,*) 'Checking wavefunctions..'
   write(uout,*) 'Computing max and L2 - norms between analytical and calculated wavefunctions.'
   write(uout,*) ''
   write(uout,'(A4, 2A15)') '', 'max-norm', 'L2-norm'
   max_err = 0.0d0

   do i = 1, 10

      associate( WF  => WF1c(1:N, i),             & ! calculated wavefunction
            &    WF0 => WF_Doolen(i-1, rvc(1:N)), & ! exact wavefunction
            &    QW        => qwc(1:N),           &
            &    L2_norm   => t1,                 &
            &    max_norm  => t2,                 &
            &    s         => t3   )

          s = merge(1.0d0, -1.0d0, real(WF(1)/WF0(1)) > 0.0d0) ! get the wavefunction sign
          L2_norm  = sqrt(real(dot_product(s*WF - WF0, (s*WF - WF0)*abs(QW)), kind=8))
          max_norm = maxval(abs(s*WF - WF0))
          max_err  = max(max_err, L2_norm, max_norm)

          write(uout,'(I4, 2ES15.2)') i-1, max_norm, L2_norm
      end associate
   enddo 
   
   call check_err(max_err, tol)

   !=========================================================
   ! Test 28. First derivatives dE/dgamma and dW/dgamma
   ! for Doolen transformed Coulomb + inverse-square potential.
   ! See Tests 26-27.
   !
   ! Exact formulas:
   !   E_n(gamma) = 2*(8*gamma - 4*n^2 - 4*n - 2) / D^2
   !   G_n(gamma) = 4*(4*n + 2)*sqrt(8*gamma - 1) / D^2
   ! where D = 4*n^2 + 4*n + 8*gamma
   !
   !   dE_n/dgamma = 64*(3*n^2 + 3*n + 1 - 2*gamma) / D^3
   !   dW_n/dgamma = 32*(2*n + 1)/(sqrt(8*gamma - 1)*D^2) &
   !               - 128*(2*n + 1)*sqrt(8*gamma - 1)/D^3
   !=========================================================
   write(uout,*) ''
   write(uout,*) 'Test 28. Doolen potential derivatives dE/dgamma and dW/dgamma'
   call announce_test('Test 28. Doolen model: derivatives')
   write(uout,*) ''
   write(uout,*) 'Checking calulates derivatives with analytical ones... '

   tol = 5.0d-4
   max_err = 0.0d0

   N = 150

   call qsolver( rmin     = 0.0d0,           & ! solved on (0; +infinity)
          &      bc_left  = BC_REGULAR,      & ! regular finite behavior at r = 0
          &      mu       = 1.0d0,           & ! atomic units
          &      N        = N,               & ! number of collocation points
          &      Emin     = 0.0d0,           & ! keep only positive energies
          &      V1c      = V_Doolen_cmplx,  & ! complex potential
          &      B1c      = B_Doolen_cmplx,  & ! complex B-coupling
          &      dV1c     = dV_Doolen_cmplx, & ! derivative subroutine for dV/dgamma
          &      dB1c     = dB_Doolen_cmplx, & ! derivative subroutine for dB/dgamma
          &      dE       = dEv,             & ! output dE/dgamma
          &      dW       = dWv,             & ! output dW/dgamma
          &      is_B_sym = .false.,         & ! B enters as plain B*d/dr
          &      sortE    = 2,               & ! sort by Re(E), descending
          &      mtype    = 'dexp',          & ! double exponential mapping
          &      beta     = 4.0d-4,          & ! mapping parameter
          &      phi      = 1.1d0            & ! complex scaling angle
          &               )

   write(uout,*) ''
   write(uout,'(A4, 4A15, 2A12)') 'n', 'dE/dg calc', 'dE/dg exact', &
          &                    'dW/dg calc',      'dW/dg exact', &
          &                    'rel_err E', 'rel_err G'

   do i = 1, 10
      associate( nlev         => i - 1,                   &
            &    gamma     => 401.0d0/8,               &
            &    sq        => sqrt(8*(401.0d0/8) - 1), &
            &    dE_calc   => dEv(1,i),                &
            &    dW_calc   => dWv(1,i),                &
            &    dE_exact  => t1,                      &
            &    dW_exact  => t2,                      &
            &    rel_err_E => t3,                      &
            &    rel_err_G => t4 )

         t0 = 4*nlev*nlev + 4*nlev + 8*gamma

         dE_exact = 64*(3*nlev*nlev + 3*nlev + 1 - 2*gamma)/t0**3
         dW_exact = 32*(2*nlev + 1)/(sq*t0**2) &
                   & - 128*(2*nlev + 1)*sq/t0**3

         rel_err_E = abs(dE_calc - dE_exact)/abs(dE_exact)
         rel_err_G = abs(dW_calc - dW_exact)/abs(dW_exact)

         max_err = max(max_err, rel_err_E, rel_err_G)

         write(uout,'(I4, 4ES15.6, 2ES12.2)') nlev, dE_calc, dE_exact, dW_calc,    & 
                                           &  dW_exact, rel_err_E, rel_err_G
      end associate
   enddo
   call check_err(max_err, tol)

   !=========================================================
   ! Test 29. Two-channel Lefebvre model.
   ! Model from:
   !   R. Lefebvre, J. Chem. Phys. 92, 2869 (1990).
   !   doi: 10.1063/1.457933
   !
   ! This example is a coupled-channel model that represents a
   ! situation in which an excited electronic diatomic state
   ! possesses bound energy levels when considered in isolation.
   ! However, all these levels become resonant once the
   ! interaction with a lower electronic state is included.
   !
   ! The excited electronic state is described by a Morse potential
   !
   !   V11 = D*(1 - exp(-beta*(r - re)))^2
   !
   ! with parameters
   !   D    = 15000 cm^-1
   !   re   = 1.6 A
   !   beta = 1.9685 A^-1
   !
   ! For the reduced mass mu = 8 u this potential supports
   ! 43 bound levels.
   !
   ! The lower electronic state is represented by a repulsive potential
   !
   !   V22 = A*exp(-gamma*(r - ra)) + B
   !
   ! with parameters
   !   A     = 18154.95 cm^-1
   !   B     = -8000 cm^-1
   !   ra    = 2.48 A
   !   gamma = 2.2039 A^-1
   !
   ! The interaction between the two states is taken as
   !
   !   V12 = V21 = epsilon
   !
   ! where
   !   epsilon = 100 cm^-1
   !
   ! The parameters for extension into the complex plane are
   !   rc  = 4.1 A
   !   phi = 0.55 rad
   !
   ! Step 1:
   !   Find bound levels of the isolated Morse potential V11.
   !
   ! Step 2:
   !   Solve the coupled 2x2 problem with exterior complex scaling
   !   and compute resonance shifts and widths.
   !
   ! Here the shifts dE and widths G are compared with reference
   ! values obtained in a previous calculation.
   !
   ! The atan and rat mappings appear to be best here.
   !=========================================================
   write(uout,*) ''
   write(uout,*) 'Test 29. Two-channel Lefebvre model'
   call announce_test('Test 29. Lefebvre model: complex potential')
   write(uout,*) '         Step 1: Compute bound levels of isolated V11 (Morse) potential'
   write(uout,*) '         Step 2: Compute energy shifts and widths for 2x2 equation'

   N  = 130

   ! Reduced mass mu = 8 u converted to cm-1 - angstrom units
   mu = 8.0d0/(5.485799090441d-4*2.1947463136314d5*0.529177210544d0**2)

   !---------------------------------------------------------
   ! Step 1. Bound levels of isolated V11 (Morse) potential.
   ! This step uses the dedicated single-channel function V_Lefebvre_V11.
   !---------------------------------------------------------

   call qsolver( rmin  = 0.3d0,      & ! equation is solved on [0.3; +infinity)
          &      mu    = mu,         & ! reduced mass
          &      N     = N,          & ! number of grid points
          &      Emax  = 15000.0d0,  & ! find levels up to dissociation energy
          &      nEf   = nEf,        & ! number of levels found
          &      V1    = V_Lefebvre_V11, & ! isolated Morse V11 channel
          &      mtype = 'rat',      & ! rational mapping appears to be best here
          &      ye    = 1.6d0,      & ! ye corresponds to potential minimum (optimal)
          &      beta  = 3.0d0,      & ! mapping parameter
          &      E     = E0v         & ! output bound energies
          &                )
   
   write(uout,*) ''
   if (nEf .ne. 43) then
      write(uout,*) 'Calculation of bound levels for Lefebvre potential failed.'
      write(uout,'(A,I3)') ' Expected 43 bound levels, found nEf =', nEf
      max_err = huge(1.0d0)
      tol     = 5.0d-9
   else
      write(uout,'(A,I3,A)') 'Found ', nEf, ' bound levels'
      write(uout,*) ''


      !---------------------------------------------------------
      ! Step 2. Coupled 2x2 complex-scaled calculation.
      ! For each bound level, find the closest coupled resonance.
      !---------------------------------------------------------
      Eref(1:43) = (/ &
              &   -0.087292571d0, -0.093189981d0, -0.099795518d0, -0.107237852d0, &
              &   -0.115678826d0, -0.125325047d0, -0.136444801d0, -0.149393437d0, &
              &   -0.164652731d0, -0.182894354d0, -0.205087208d0, -0.232690048d0, &
              &   -0.268025072d0, -0.315086954d0, -0.381628102d0, -0.485897541d0, &
              &   -0.677650818d0, -1.042755292d0, -1.359122980d0, -0.448301388d0, &
              &    1.249178515d0,  0.134543381d0,  0.420291461d0,  0.438361976d0, &
              &    0.500946002d0,  0.160672240d0,  0.777356369d0,  0.171026657d0, &
              &    0.337699249d0,  0.684408003d0,  0.392229115d0,  0.206299735d0, &
              &    0.394262635d0,  0.570789669d0,  0.545363360d0,  0.425973338d0, &
              &    0.345358686d0,  0.334721427d0,  0.363106776d0,  0.396741400d0, &
              &    0.419886969d0,  0.430956042d0,  0.434340653d0 /)

      Wref(1:43) = (/ &
              &    0.000000000d0, -0.000000000d0,  0.000000000d0, -0.000000000d0, &
              &   -0.000000000d0, -0.000000000d0, -0.000000000d0,  0.000000000d0, &
              &    0.000000000d0,  0.000000000d0,  0.000000000d0,  0.000000001d0, &
              &    0.000000067d0,  0.000003831d0,  0.000141992d0,  0.003360600d0, &
              &    0.049369232d0,  0.428457446d0,  2.011863886d0,  4.324079290d0, &
              &    2.584850822d0,  0.110260688d0,  2.460009749d0,  0.001074642d0, &
              &    1.842657883d0,  0.176762961d0,  0.872554566d0,  1.098180452d0, &
              &    0.021920994d0,  0.532308833d0,  0.952924079d0,  0.399930647d0, &
              &    0.003684260d0,  0.179075286d0,  0.463521438d0,  0.514805938d0, &
              &    0.371188436d0,  0.194541674d0,  0.075349961d0,  0.020167464d0, &
              &    0.002933399d0,  0.000074772d0,  0.000010726d0 /)

      ! Literature values from Lefebvre (1990).
      ! The paper tabulates Gamma/2, while qsolver uses W = Gamma,
      ! therefore Wlit is stored as 2*(Gamma/2).
      Elit = huge(1.0d0)
      Wlit = huge(1.0d0)

      Elit(1:31) = (/ &
              &   -0.0873d0, -0.0932d0, -0.0998d0, -0.1072d0, -0.1157d0, &
              &   -0.1253d0, -0.1364d0, -0.1494d0, -0.1646d0, -0.1829d0, &
              &   -0.2050d0, -0.2326d0, -0.2679d0, -0.3150d0, -0.3815d0, &
              &   -0.4857d0, -0.6772d0, -1.0422d0, -1.3587d0, -0.4517d0, &
              &    1.2490d0,  0.1368d0,  0.4184d0,  0.4405d0,  0.4988d0, &
              &    0.1622d0,  0.7768d0,  0.1700d0,  0.3390d0,  0.6843d0, &
              &    0.3910d0 /)

      Wlit(14:31) = 2.0d0 * (/ &
              &    0.000002d0, 0.000071d0, 0.001667d0, 0.02461d0, 0.2137d0, &
              &    1.0043d0,  2.1597d0,  1.2954d0, 0.5419d0, 1.2307d0, &
              &    0.000638d0, 0.9222d0, 0.8724d0, 0.4378d0, 0.5483d0, &
              &    0.01059d0, 0.2672d0, 0.4768d0 /)

      write(uout,'(A4, A12, A14, A13, A9, A14, A11, A10)') 'v', 'Ebound', 'dEcalc', &
             &               'dEcalc-dEref', 'dElit', 'Wcalc', 'Wcalc-Wref', 'Wlit'


      max_err = 0.0d0
      tol = 5.0d-9
      missing_levels = 0

      N = 350

      do i = 1, 43
         call qsolver( rmin  = 0.0d0,            & ! equation is solved on [0.0; +infinity)
                &      mu    = mu,               & ! reduced mass
                &      N     = N,                & ! number of grid points
                &      nEq   = 2,                & ! two coupled channels
                &      E0    = E0v(i),           & ! target around Morse bound level
                &      nEf   = nEf,              & ! number of levels found
                &      Vc    = V_Lefebvre_cmplx, & ! coupled 2x2 Lefebvre potential
                &      mtype = 'atan',           & ! mapping for semi-infinite interval
                &      ye    = 1.6d0,            & ! mapping center
                &      beta  = 0.6d0,            & ! mapping parameter
                &      rc    = 4.1d0,            & ! ECS start radius (A)
                &      phi   = 0.55d0,            & ! ECS rotation angle (rad)
                &      E     = Ev,               & ! output resonance energy
                &      W     = Wv                & ! output resonance width
                &                )

         if (nEf >= 1) then
            t1 = Ev(1) - E0v(i)
            t2 = t1 - Eref(i)
            t3 = Wv(1) - Wref(i)
            max_err = max(max_err, abs(t2), abs(t3))

            sElit = ''
            sWlit = ''
            if (Elit(i) < huge(1.0d0)/2.0d0) write(sElit, '(F9.4)') Elit(i)
            if (Wlit(i) < huge(1.0d0)/2.0d0) write(sWlit, '(F10.6)') Wlit(i)

            write(uout,'(I4, F12.3, F14.9, ES13.2, A9, F14.9, ES11.2, A10)') &
               & i-1, E0v(i), t1, t2, sElit, Wv(1), t3, sWlit
         else
            write(uout,'(I4, A)') i-1, ' no resonance found near target E0'
            missing_levels = missing_levels + 1
            max_err = max(max_err, 2.0d0*tol)
         endif
       enddo

      if (missing_levels > 0) then
         write(uout,'(A,I4,A)') 'ERROR: ', missing_levels, &
            ' Lefebvre resonance levels were not found.'
      endif
   endif

   !===================================================
   ! If max_err > tol, then the test fails
   !===================================================
   call check_err(max_err, tol)

   !=========================================================
   ! Test 30. The same Lefebvre model as in Test 29, but now
   ! the real potential matrix V_Lefebvre is used.
   !
   ! In this case qsolver automatically extrapolates the real
   ! potential into the complex plane. The same reference values
   ! are used as in Test 29. Since automatic extrapolation is
   ! less accurate than providing the complex potential directly,
   ! a looser tolerance is used here.
   !=========================================================
   write(uout,*) ''
   write(uout,*) 'Test 30. Two-channel Lefebvre model'
   call announce_test('Test 30. Lefebvre model: real potential')
   write(uout,*) '         Step 2: Compute energy shifts and widths for 2x2 equation'
   write(uout,*) '         Real potential automatically extrapolated to the complex plane'

   write(uout,'(A4, A12, A14, A13, A9, A14, A11, A10)') 'v', 'Ebound', 'dEcalc', &
          &               'dEcalc-dEref', 'dElit', 'Wcalc', 'Wcalc-Wref', 'Wlit'

   max_err = 0.0d0
   tol = 1.0d-6
   missing_levels = 0

   N = 350

   do i = 1, 43
      call qsolver( rmin  = 0.0d0,      & ! equation is solved on [0.0; +infinity)
             &      mu    = mu,         & ! reduced mass
             &      N     = N,          & ! number of grid points
             &      nEq   = 2,          & ! two coupled channels
             &      E0    = E0v(i),     & ! target around Morse bound level
             &      nEf   = nEf,        & ! number of levels found
             &      V     = V_Lefebvre, & ! real coupled 2x2 Lefebvre potential
             &      mtype = 'atan',     & ! mapping for semi-infinite interval
             &      ye    = 1.6d0,      & ! mapping center
             &      beta  = 0.6d0,      & ! mapping parameter
             &      rc    = 4.1d0,      & ! ECS start radius (A)
             &      phi   = 0.55d0,     & ! ECS rotation angle (rad)
             &      E     = Ev,         & ! output resonance energy
             &      W     = Wv          & ! output resonance width
             &                )

      if (nEf >= 1) then
         t1 = Ev(1) - E0v(i)
         t2 = t1 - Eref(i)
         t3 = Wv(1) - Wref(i)
         max_err = max(max_err, abs(t2), abs(t3))

         sElit = ''
         sWlit = ''
         if (Elit(i) < huge(1.0d0)/2.0d0) write(sElit, '(F9.4)') Elit(i)
         if (Wlit(i) < huge(1.0d0)/2.0d0) write(sWlit, '(F10.6)') Wlit(i)

         write(uout,'(I4, F12.3, F14.9, ES13.2, A9, F14.9, ES11.2, A10)') &
            & i-1, E0v(i), t1, t2, sElit, Wv(1), t3, sWlit
      else
         write(uout,'(I4, A)') i-1, ' no resonance found near target E0'
         missing_levels = missing_levels + 1
         max_err = max(max_err, 2.0d0*tol)
      endif
   enddo

   if (missing_levels > 0) then
      write(uout,'(A,I4,A)') 'ERROR: ', missing_levels, &
         ' Lefebvre resonance levels were not found.'
   endif

   !===================================================
   ! If max_err > tol, then the test fails
   !===================================================
   call check_err(max_err, tol)

   !=========================================================
   ! Test 31. Two-channel Sitnikov-Tolstikhin model.
   !
   !   G. V. Sitnikov and O. I. Tolstikhin,
   !   Phys. Rev. A 67, 032714 (2003), Sec. B, Eq. (92a)
   !
   ! Complex 2x2 potential:
   !   V11 = 15*exp(-0.5*r)
   !   V22 = 15*(r^2 - r - 1)*exp(-r) + 15
   !   V12 = V21 = 5*r*exp(-r)
   !   mu = 1
   !
   ! All three resonances from cited article:
   !   E = 6.9682224478, 14.4530939667, 18.2847366077
   !   W = 0.1354784309, 0.0115049984, 0.1185907300
   !
   ! The calculated energies and widths are accurate to about 12-13 digits,
   ! which is more accurate than the values reported in the article.
   !
   ! In this test the complex potential is provided directly.
   ! Here the complex-scaling angle is phi = 0.50 rad.
   ! The atan mapping is used with beta = 0.20.
   !
   ! All three resonances are computed in a single diagonalization.
   ! Nonphysical roots are removed by the built-in qsolver filter.
   ! Since no energy range or target values are provided, qsolver
   ! initially computes all (2*N) eigenvalues of the Hamiltonian matrix,
   ! after which the internal filter removes the unphysical ones,
   ! leaving only the three physical eigenstates.
   !=========================================================
   write(uout,*) ''
   write(uout,*) 'Test 31. Two-channel Sitnikov model'
   call announce_test('Test 31. Two-channel Sitnikov model: complex potential')
   write(uout,*) '         Compute three resonances in one run and compare with reference values'

   Elit(1:3) = (/ 6.9682224478d0, 14.4530939667d0, 18.2847366077d0 /)
   Wlit(1:3) = (/ 0.1354784309d0, 0.0115049984d0, 0.1185907300d0 /)

   N = 60
   mu = 1.0d0
   tol = 1.0d-9
   max_err = 0.0d0

   write(uout,*) ''
   write(uout,'(A2, A14, A18, A10, A14, A18, A10)') '', 'E lit', 'E calc', &
          &                                              '|dE|', 'W lit', 'W calc', '|dW|'

   Ev = 0.0d0
   Wv = 0.0d0

   call qsolver( rmin  = 0.0d0,            & ! equation is solved on [0.0; +infinity)
          &      mu    = mu,               & ! effective mass
          &      N     = N,                & ! number of grid points
          &      nEq   = 2,                & ! two coupled channels
          &      nEf   = nEf,              & ! number of eigenvalues found
          &      Vc    = V_Sitnikov_cmplx, & ! complex 2x2 potential matrix
          &      phi   = 0.50d0,           & ! complex-scaling angle
          &      mtype = 'atan',           & ! mapping for semi-infinite interval
          &      beta  = 0.2d0,            & ! mapping parameter
          &      E     = Ev,               & ! output resonance energies
          &      W     = Wv                & ! output resonance widths
          &                )

   if (nEf < 3) then
      write(uout, *) 'Not all resonances were found.'
      max_err = max(max_err, 2.0d0*tol)
   endif

   do i = 1, 3
      t1 = abs(Ev(i) - Elit(i))
      t2 = abs(Wv(i) - Wlit(i))
      max_err = max(max_err, t1, t2)
      write(uout,'(I2, F14.10, F18.12, ES10.2, F14.10, F18.12, ES10.2)') &
         & i, Elit(i), Ev(i), t1, Wlit(i), Wv(i), t2
   enddo
   call check_err(max_err, tol)

   !=========================================================
   ! Test 32. The same two-channel Sitnikov-Tolstikhin model as
   ! in Test 31, but now the real potential matrix is provided.
   !
   ! qsolver automatically extrapolates the real potential into the
   ! complex plane. In this case the optimal numerical parameters
   ! differ slightly from the complex-potential case.
   !
   ! Here the complex-scaling angle is phi = 0.50 rad.
   ! The atan mapping is used with beta = 0.10.
   ! 
   ! All three resonances are computed in a single diagonalization.
   ! Nonphysical roots are removed by the built-in qsolver filter.
   ! Since no energy range or target values are provided, qsolver
   ! initially computes all (2*N) eigenvalues of the Hamiltonian matrix,
   ! after which the internal filter removes the unphysical ones,
   ! leaving only the three physical eigenstates.
   !=========================================================
   write(uout,*) ''
   write(uout,*) 'Test 32. Two-channel Sitnikov model'
   call announce_test('Test 32. Two-channel Sitnikov model: real potential')
   write(uout,*) '         Compute three resonances in one run and compare with reference values'

   N = 80
   mu = 1.0d0
   tol = 1.0d-7
   max_err = 0.0d0

   write(uout,*) ''
   write(uout,'(A2, A14, A18, A10, A14, A18, A10)') '', 'E lit', 'E calc', &
          &                                              '|dE|', 'W lit', 'W calc', '|dW|'

   Ev = 0.0d0
   Wv = 0.0d0

   call qsolver( rmin  = 0.0d0,           & ! equation is solved on [0.0; +infinity)
          &      mu    = mu,              & ! effective mass
          &      N     = N,               & ! number of grid points
          &      nEq   = 2,               & ! two coupled channels
          &      nEf   = nEf,             & ! number of eigenvalues found
          &      V     = V_Sitnikov_real, & ! real 2x2 potential matrix
          &      phi   = 0.50d0,          & ! complex-scaling angle
          &      mtype = 'atan',          & ! mapping for semi-infinite interval
          &      beta  = 0.1d0,           & ! mapping parameter
          &      E     = Ev,              & ! output resonance energies
          &      W     = Wv               & ! output resonance widths
          &                )

   if (nEf < 3) then
      write(uout, *) 'Not all resonances were found.'
      max_err = max(max_err, 2.0d0*tol)
   endif

   do i = 1, 3
      t1 = abs(Ev(i) - Elit(i))
      t2 = abs(Wv(i) - Wlit(i))
      max_err = max(max_err, t1, t2)
      write(uout,'(I2, F14.10, F18.12, ES10.2, F14.10, F18.12, ES10.2)') &
         & i, Elit(i), Ev(i), t1, Wlit(i), Wv(i), t2
   enddo
   call check_err(max_err, tol)

   !===================================================
   ! Final summary
   !===================================================
   write(uout,'(A, I3)') 'Number of passed tests: ', npass 
   write(uout,'(A, I3)') 'Number of failed tests: ', nfail
   close(uout)

   write(*,'(A, I3)') 'Number of passed tests: ', npass
   write(*,'(A, I3)') 'Number of failed tests: ', nfail

   write(*,*) ''
   write(*,*) 'Detailed test results are in qtest_out.txt'
   write(*,*) ''

   if (nfail > 0) stop 1

contains
   !===================================================
   ! Print the current test name on the screen.
   !===================================================
   subroutine announce_test(title)
      character(*), intent(in) :: title
      write(*,'(A)') 'Running: ' // trim(title)
   end subroutine

   !===================================================
   ! Subroutine for comparing max_err with the tolerance.
   ! If max_err > tol, then the test fails.
   !===================================================
   subroutine check_err(max_err, tol)
      real(8), intent(in)  :: max_err
      real(8), intent(in)  :: tol

      write(uout,*) ''
      if (max_err <= tol) then
         write(uout, '(A, ES9.2, A, ES9.2)') 'TEST PASSED: max_err = ', max_err, &
                                     ' < tol = ', tol
         npass = npass + 1
         write(*,'(A)') 'PASSED'
      else
         write(uout, '(A, ES9.2, A, ES9.2)') 'TEST FAILED: max_err = ', max_err, &
                                     ' > tol = ', tol
         nfail = nfail + 1
         write(*,'(A)') 'FAILED'
      endif
      write(uout,*) ''
   end subroutine
end program
