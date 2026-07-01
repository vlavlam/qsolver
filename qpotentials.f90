!======================================================================
! qpotentials.f90 -- module of model problems for the qsolver test
! program qtest.f90. It collects the model potentials V(r) (real and
! complex, single- and multi-channel), the B-coupling and G-weight
! functions, their derivatives with respect to fitting parameters, and
! the auxiliary special functions and exact reference wavefunctions
! (Legendre polynomials, the 1F1 hypergeometric polynomial and the
! analytic Doolen wavefunction) used to verify the computed solutions.
!
! Version: 1.0
! SPDX-License-Identifier: MIT
! Copyright (c) 2026 Vladimir Meshkov, Asen Pashov, Andrey Stolyarov
! qsolver is distributed under the MIT License; see the LICENSE file
! distributed with this program for the full license text.
!======================================================================
module potentials
implicit none
   integer :: JRot = 0  ! rotational quantum number
contains
   !=========================================================
   ! Harmonic oscillator potential
   !=========================================================
   function V_Harm_Osc(r)
     real(8), intent(in) :: r
     real(8)             :: V_Harm_Osc
   
     V_Harm_Osc = r*r/2
   
   end function

   !=========================================================
   ! Coulomb potential for electron-proton attraction in atomic units
   !=========================================================
   function V_Coulomb(r)
     real(8), intent(in) :: r
     real(8)             :: V_Coulomb
   
     V_Coulomb = -1/r
   
   end function

   !=========================================================
   ! B-coupling for the singular Bessel test problem:
   !   y'' + (1/x)*y' + lambda*y = 0.
   ! In qsolver notation with mu = 0.5 and is_B_sym = .false.,
   !   B(x) = -1/x.
   !=========================================================
   function B_Bessel(x) result(B)
      real(8), intent(in) :: x
      real(8)             :: B

      B = -1.0d0/x

   end function B_Bessel

   !=========================================================
   ! B-coupling for the Jacobi/Gegenbauer test problem:
   !   -(1 - x^2)*y'' + 2*(a + 1)*x*y' = E*y.
   ! After division by (1 - x^2):
   !   -y'' + B(x)*y' = E*G(x)*y.
   ! Here a = 0.
   !=========================================================
   function B_Jacobi(x) result(B)
      real(8), intent(in) :: x
      real(8)             :: B

      real(8), parameter :: a = 0.0d0

      B = 2.0d0*(a + 1.0d0)*x/(1.0d0 - x*x)

   end function B_Jacobi

   !=========================================================
   ! Diagonal matrix G(x) for the Jacobi/Gegenbauer test.
   !=========================================================
   function G_Jacobi(x) result(G)
      real(8), intent(in) :: x
      real(8)             :: G

      G = 1.0d0/(1.0d0 - x*x)

   end function G_Jacobi

   !=========================================================
   ! Derivative dB/da for the Jacobi/Gegenbauer test.
   !=========================================================
   subroutine dB_Jacobi(x, dBdp)
      real(8), intent(in)  :: x
      real(8), intent(out) :: dBdp(:)

      dBdp = 0.0d0
      if (size(dBdp) >= 1) dBdp(1) = 2.0d0*x/(1.0d0 - x*x)

   end subroutine dB_Jacobi

   !=========================================================
   ! Legendre polynomial P_n(x).
   ! Used for the analytical eigenfunctions in the Jacobi test.
   !=========================================================
   pure elemental function legendre_poly(n, x) result(Pn)
      integer, intent(in) :: n
      real(8), intent(in) :: x
      real(8)             :: Pn

      integer :: k
      real(8) :: Pkm2, Pkm1, Pk

      if (n == 0) then
         Pn = 1.0d0
      else if (n == 1) then
         Pn = x
      else
         Pkm2 = 1.0d0
         Pkm1 = x
         do k = 2, n
            Pk   = ((2*k - 1)*x*Pkm1 - (k - 1)*Pkm2)/k
            Pkm2 = Pkm1
            Pkm1 = Pk
         enddo
         Pn = Pkm1
      endif

   end function legendre_poly

   !=========================================================
   ! Morse potential with a minimum at r = re and asymptote De as r -> +inf.
   ! V(r) = De * (1 - exp(-a*(r - re)))^2
   !=========================================================
   function V_Morse(r) 
      real(8), intent(in) :: r
      real(8)             :: V_Morse
      real(8) :: t
  
      ! Parameters:
      real(8), parameter :: De = 10000.0d0
      real(8), parameter :: a  = 2 * sqrt(2.0d0)
      real(8), parameter :: re = 3.0d0
  
      t = exp( -a * (r - re) )
      V_Morse = De * (1.0d0 - t)**2

   end function V_Morse

   !=========================================================
   ! Derivatives of Morse potential with respect to parameters De and a
   !=========================================================
   subroutine dV_Morse(r, dVdp)
      real(8), intent(in)  :: r
      real(8), intent(out) :: dVdp(:)
      real(8)              :: t, dr
     
      ! Parameters:
      real(8), parameter :: De = 10000.0d0
      real(8), parameter :: a  = 2.0d0 * sqrt(2.0d0)
      real(8), parameter :: re = 3.0d0
     
      dr = r - re
      t  = exp( -a * dr )
     
      dVdp = 0.0d0
      if (size(dVdp) >= 1) dVdp(1) = (1.0d0 - t)**2
      if (size(dVdp) >= 2) dVdp(2) = 2.0d0 * De * (1.0d0 - t) * t * dr
      
   end subroutine dV_Morse

   !=========================================================
   ! Constant generalized-metric function for the Morse mass
   ! derivative test. The physical mass parameter M is introduced as
   !   (-1/2 d2/dr2 + M*V) psi = E*M*psi.
   ! At the test point M = 1, G(r) = M = 1.
   !=========================================================
   function G_Morse_mass(r)
      real(8), intent(in) :: r
      real(8)             :: G_Morse_mass

      G_Morse_mass = 1.0d0 + 0.0d0*r

   end function G_Morse_mass

   !=========================================================
   ! d(M*V_Morse)/dM at M = 1.
   !=========================================================
   subroutine dV_Morse_mass(r, dVdp)
      real(8), intent(in)  :: r
      real(8), intent(out) :: dVdp(:)

      dVdp = 0.0d0
      if (size(dVdp) >= 1) dVdp(1) = V_Morse(r)

   end subroutine dV_Morse_mass

   !=========================================================
   ! dG/dM for G(r) = M at M = 1.
   !=========================================================
   subroutine dG_Morse_mass(r, dGdp)
      real(8), intent(in)  :: r
      real(8), intent(out) :: dGdp(:)

      dGdp = 0.0d0
      if (size(dGdp) >= 1) dGdp(1) = 1.0d0 + 0.0d0*r

   end subroutine dG_Morse_mass

   !=========================================================
   ! Double-minimum potential from B. R. Johnson, J. Chem. Phys. 67, 4086–4093 (1977)
   !=========================================================
   function V_DMin(r)
      real(8), intent(in) :: r
      real(8)             :: V_DMin
      real(8) :: t1, t2

      ! Parameters (as in the paper/table):
      ! V(r) = D * (1 - exp[-B*(r - xa)])^2 + A * exp[-C*(r - xb)^2]

      real(8), parameter :: D  = 31250.0d0          ! cm^-1
      real(8), parameter :: B  = 1.5403756164035d0  ! A^-1
      real(8), parameter :: xa = 1.5d0              ! A
     
      real(8), parameter :: A  = 10000.0d0          ! cm^-1
      real(8), parameter :: C  = 200.0d0            ! A^-2
      real(8), parameter :: xb = 1.6d0              ! A
     
      t1 = exp( -B * (r - xa) )
      t2 = exp( -C * (r - xb)**2 )
     
      V_DMin = D * (1.0d0 - t1)**2 + A * t2

   end function V_DMin

   !=========================================================
   ! Symmetric Double-minimum potential from
   ! B. R. Johnson, J. Chem. Phys. 67, 4086–4093 (1977)
   !=========================================================
   function V_DMin_sym(r)
      real(8), intent(in) :: r
      real(8)             :: V_DMin_sym

      V_DMin_sym = (r*r - 1)**2

   end function V_DMin_sym

   !=========================================================
   ! Lennard-Jones 6-12 potential from
   ! V. Meshkov, A. Stolyarov and R. J. Le Roy, Phys. Rev. A 78, 052510 (2008)
   !=========================================================
   function V_LJ_6_12(r)
      real(8), intent(in) :: r
      real(8) :: V_LJ_6_12

      real(8) :: t1

      real(8), parameter :: re  = 1.0d0    ! angstrom
      real(8), parameter :: De  = 3761.0d0 ! cm-1

      t1 = (re/r)**6
      
      V_LJ_6_12 = De*(t1*t1 - 2*t1) 

   end function V_LJ_6_12

   !=========================================================
   ! Lennard-Jones 3-6 potential from
   ! V. Meshkov, A. Stolyarov and R. J. Le Roy, Phys. Rev. A 78, 052510 (2008)
   !=========================================================
   function V_LJ_3_6(r)
      real(8), intent(in) :: r
      real(8) :: V_LJ_3_6

      real(8) :: t1

      real(8), parameter :: re  = 1.0d0    ! angstrom
      real(8), parameter :: De  = 107000d0 ! cm-1

      t1 = (re/r)**3
      
      V_LJ_3_6 = De*(t1*t1 - 2*t1) 

   end function V_LJ_3_6

   !=========================================================
   ! 7.5*r^2*exp(-r) potential from
   ! R A Bain, J N Burdsley, B R Junker and C V Sukumar
   ! J. Phys. B. 7(16), 2189 (1974)
   !=========================================================
   function V_BBJS_cmplx(r)
      complex(8), intent(in) :: r
      complex(8)             :: V_BBJS_cmplx
    
      V_BBJS_cmplx = 7.5d0*r**2*exp(-r)

   end function

   function V_BBJS_real(r)
      real(8), intent(in) :: r
      real(8)             :: V_BBJS_real
    
      V_BBJS_real = 7.5d0*r**2*exp(-r)

   end function

   !=========================================================
   ! Complex Li2 PB^1 Pi_u state potential from
   ! Y. Huang and R. J. Le Roy
   ! J. Chem. Phys. 119, 7398 (2003)
   ! r in angstrom. E in cm-1.
   !=========================================================
   function V_Li2B_cmplx(r)
      complex(8), intent(in) :: r
      complex(8)             :: V_Li2B_cmplx
   
      integer            :: i, n, ilr
      integer, parameter :: p = 3, nl = 9, ns = 5, nlambda = 2, nlr = 4
      integer, parameter :: lr_pow(1:nlr) = (/ 3, 6, 8, 10 /)
      real(8), parameter :: req = 2.93617142d0
      real(8), parameter :: A = 3853.980130d0, B = 9982.620189d0
      real(8), parameter :: Dinf = 14903.983468d0, rhod = 0.464700505d0
      real(8), parameter :: mu = 3.5080019999d0/5.4857990946d-4/2.194746313708d5/0.52917721092d0**2
      real(8), parameter :: beta(0:nl) = (/ 0.970911966d0, 0.2075358d0, 0.1751542d0, 0.188843d0, &
                                          & 0.15648d0, 0.252d0, -2.185d0, 6.91598d0, -9.6903477d0, 4.7186d0 /)
      real(8), parameter :: lam(0:nlambda) = (/ 0.00031365d0, -0.000215d0, 0.000314d0 /)
      real(8), parameter :: lr_coef(1:nlr) = (/ 1.788d5, -6.97586d6, -1.378d8, -3.445d9 /)
      complex(8)         :: y, betv, Hrot, V0, V_LR, V_BOB, V_Lam, Dm
   
      Hrot = (Jrot*(Jrot + 1.0d0) - 1.0d0)/(2.0d0*mu*r**2)
   
      if (real(r,8) >= req) then
         n = nl
      else
         n = ns
      end if
   
      y = (r**p - req**p)/(r**p + req**p)
   
      betv = (0.0d0, 0.0d0)
      do i = 0,n
         betv = betv + beta(i)*y**i
      end do
   
      V0 = A*exp(-2.0d0*betv*(r - req)) - B*exp(-betv*(r - req))
   
      V_LR = Dinf
      do ilr = 1,nlr
         i  = lr_pow(ilr)
         Dm = exp(-3.968424576d0*(rhod*r/dble(i)) - 0.3892460101d0*(rhod*r)**2/sqrt(dble(i)))
         Dm = (1.0d0 - Dm)**i
         V_LR = V_LR + Dm*lr_coef(ilr)/r**i
      end do
   
      V_Lam = (0.0d0, 0.0d0)
      do i = 0,nlambda
         V_Lam = V_Lam + lam(i)*y**i
      end do
      V_Lam = V_Lam*Jrot*(Jrot + 1.0d0)/(2.0d0*mu*r**2)**2
   
      V_BOB = (0.0d0, 0.0d0)
   
      V_Li2B_cmplx = V0 + V_LR + Hrot + V_Lam + V_BOB 
   end function V_Li2B_cmplx

   !=========================================================
   ! Real Li2 PB^1 Pi_u state potential from
   ! Y. Huang and R. J. Le Roy
   ! J. Chem. Phys. 119, 7398 (2003)
   ! r in angstrom. E in cm-1.
   !=========================================================
   function V_Li2B_real(r)
      real(8), intent(in) :: r
      real(8)             :: V_Li2B_real

      integer            :: i, n, ilr
      integer, parameter :: p = 3, nl = 9, ns = 5, nlambda = 2, nlr = 4
      integer, parameter :: lr_pow(1:nlr) = (/ 3, 6, 8, 10 /)
      real(8), parameter :: req = 2.93617142d0
      real(8), parameter :: A = 3853.980130d0, B = 9982.620189d0
      real(8), parameter :: Dinf = 14903.983468d0, rhod = 0.464700505d0
      real(8), parameter :: mu = 3.5080019999d0/5.4857990946d-4/2.194746313708d5/0.52917721092d0**2
      real(8), parameter :: beta(0:nl) = (/ 0.970911966d0, 0.2075358d0, 0.1751542d0, 0.188843d0, &
                                          & 0.15648d0, 0.252d0, -2.185d0, 6.91598d0, -9.6903477d0, 4.7186d0 /)
      real(8), parameter :: lam(0:nlambda) = (/ 0.00031365d0, -0.000215d0, 0.000314d0 /)
      real(8), parameter :: lr_coef(1:nlr) = (/ 1.788d5, -6.97586d6, -1.378d8, -3.445d9 /)
      real(8)            :: y, betv, Hrot, V0, V_LR, V_BOB, V_Lam, Dm

      Hrot = (Jrot*(Jrot + 1.0d0) - 1.0d0)/(2.0d0*mu*r**2)

      if (r >= req) then
         n = nl
      else
         n = ns
      end if

      y = (r**p - req**p)/(r**p + req**p)

      betv = 0.0d0
      do i = 0,n
         betv = betv + beta(i)*y**i
      end do

      V0 = A*exp(-2.0d0*betv*(r - req)) - B*exp(-betv*(r - req))

      V_LR = Dinf
      do ilr = 1,nlr
         i  = lr_pow(ilr)
         Dm = exp(-3.968424576d0*(rhod*r/dble(i)) - 0.3892460101d0*(rhod*r)**2/sqrt(dble(i)))
         Dm = (1.0d0 - Dm)**i
         V_LR = V_LR + Dm*lr_coef(ilr)/r**i
      end do

      V_Lam = 0.0d0
      do i = 0,nlambda
         V_Lam = V_Lam + lam(i)*y**i
      end do
      V_Lam = V_Lam*Jrot*(Jrot + 1.0d0)/(2.0d0*mu*r**2)**2

      V_BOB = 0.0d0

      V_Li2B_real = V0 + V_LR + Hrot + V_Lam + V_BOB
   end function V_Li2B_real

   !=========================================================
   ! Complex potential for the transformed Doolen equation.
   ! After psi(r) = r**lambda*phi(r),
   !   V(r) = 1/r
   ! with gamma = 401/8 = 50.125 and lambda = 1/2 - 10*i.
   !=========================================================
   function V_Doolen_cmplx(r) result(V)
      complex(8), intent(in) :: r
      complex(8)             :: V

      V = 1/r

   end function V_Doolen_cmplx

   !=========================================================
   ! Complex B-coupling for the transformed Doolen equation.
   ! After psi(r) = r**lambda*phi(r),
   !   B(r) = -2*lambda/r
   ! with gamma = 401/8 = 50.125 and lambda = 1/2 - 10*i.
   ! Used as a simple B*d/dr term: set is_B_sym = .false.
   !=========================================================
   function B_Doolen_cmplx(r) result(B)
      complex(8), intent(in) :: r
      complex(8)             :: B

      real(8),    parameter :: gamma = 401.0d0/8.0d0
      complex(8), parameter :: lam = cmplx(0.5d0, -0.5d0*sqrt(8.0d0*gamma - 1.0d0), 8)

      B = -2*lam/r

   end function B_Doolen_cmplx

   !=========================================================
   ! Derivative dV/dgamma for transformed Doolen equation.
   ! V(r) = 1/r does not depend on gamma.
   !=========================================================
   subroutine dV_Doolen_cmplx(r, dVdp)
      complex(8), intent(in)  :: r
      ! r is unused because dV/dgamma is constant for this model.
      complex(8), intent(out) :: dVdp(:)

      dVdp = (0.0d0, 0.0d0)*r

   end subroutine dV_Doolen_cmplx

   !=========================================================
   ! Derivative dB/dgamma for transformed Doolen equation.
   ! B(r) = -2*lambda/r,
   ! lambda = (1 - i*sqrt(8*gamma - 1))/2.
   !=========================================================
   subroutine dB_Doolen_cmplx(r, dBdp)
      complex(8), intent(in)  :: r
      complex(8), intent(out) :: dBdp(:)

      real(8),    parameter :: gamma = 401.0d0/8.0d0
      real(8),    parameter :: s = sqrt(8.0d0*gamma - 1.0d0)
      complex(8), parameter :: dlam_dgamma = cmplx(0.0d0, -2.0d0/s, 8)

      dBdp = (0.0d0, 0.0d0)
      if (size(dBdp) >= 1) dBdp(1) = -2.0d0*dlam_dgamma/r

   end subroutine dB_Doolen_cmplx

   !=========================================================
   ! Finite polynomial for 1F1(-n,b,z)
   !=========================================================
   pure elemental function hyp1f1_negint(n, b, z) result(val)
      integer,    intent(in) :: n
      complex(8), intent(in) :: b, z
      complex(8)             :: val

      complex(8) :: term
      integer    :: k

      val  = (1.0d0, 0.0d0)
      term = (1.0d0, 0.0d0)

      do k = 1, n
         term = term * (cmplx(dble(-n + k - 1), 0.0d0, kind=8)/(b + dble(k - 1))) * (z/dble(k))
         val  = val + term
      enddo
   end function hyp1f1_negint

   !=========================================================
   ! Exact transformed Doolen wavefunction for Test 26:
   !
   !   psi(r) = r**lambda * phi(r)
   !
   ! where phi(r) satisfies the transformed equation used by qsolver.
   !
   ! Returned value:
   !   WF_Doolen(n,r) = phi_n(r), normalized so that
   !   integral_0^infinity phi_n(r)**2 dr = 1 by analytic continuation.
   !
   ! Exact outgoing (Siegert) branch used for ECS comparison:
   !   phi_n(r) = exp(+kappa*r) * 1F1(-n, 2*lambda, -2*kappa*r)
   !
   ! with
   !   lambda = (1 - i*sqrt(8*gamma - 1))/2 = 1/2 - 10*i
   !   gamma  = 401/8
   !   kappa  = 1/(n + lambda)
   !
   ! n = 0,1,2,...
   !=========================================================
   pure elemental function WF_Doolen(n, r) result(wf)
      integer,    intent(in) :: n
      complex(8), intent(in) :: r
      complex(8)             :: wf

      real(8),    parameter :: gamma = 401.0d0/8.0d0
      complex(8), parameter :: lam = cmplx(0.5d0, -0.5d0*sqrt(8.0d0*gamma - 1.0d0), 8)

      complex(8) :: kappa, z
      complex(8) :: ck, cj, norm
      real(8)    :: gamma_int
      integer    :: j, k, m

      if (n < 0) then
         wf = (0.0d0, 0.0d0)
         return
      endif

      kappa = 1.0d0 / ( dble(n) + lam )
      z     = 2.0d0 * kappa * r

      ! Compute the analytical c-norm. Since 1F1(-n,2*lambda,-2*kappa*r)
      ! is a finite polynomial, the norm is a finite double sum.
      norm = (0.0d0, 0.0d0)
      ck   = (1.0d0, 0.0d0)

      do k = 0, n
         cj = (1.0d0, 0.0d0)

         do j = 0, n
            gamma_int = 1.0d0
            do m = 1, k + j
               gamma_int = gamma_int*dble(m)
            enddo

            norm = norm + ck*cj*gamma_int/(-2.0d0*kappa)**(k + j + 1)

            if (j < n) then
               cj = cj * (cmplx(dble(-n + j), 0.0d0, 8)/(2.0d0*lam + dble(j))) &
                  &    * (-2.0d0*kappa/dble(j + 1))
            endif
         enddo

         if (k < n) then
            ck = ck * (cmplx(dble(-n + k), 0.0d0, 8)/(2.0d0*lam + dble(k))) &
               &    * (-2.0d0*kappa/dble(k + 1))
         endif
      enddo

      wf = exp(+kappa*r) * hyp1f1_negint(n, 2.0d0*lam, -z) / sqrt(norm)
   end function WF_Doolen

   !=========================================================
   ! Two-channel model from:
   !   R. Lefebvre, J. Chem. Phys. 92, 2869 (1990).
   !   Diatomic predissociative widths and shifts from multichannel
   !   Siegert quantization and asymptotic analysis of resonance
   !   wave functions. doi: 10.1063/1.457933
   !
   ! V11 = D*(1 - exp(-beta*(r-re)))^2
   ! V22 = A*exp(-gamma*(r-ra)) + B
   ! V12 = V21 = epsilon
   !
   !=========================================================
   ! Isolated Morse channel V11 used in Step 1 of Tests 29-30.
   !=========================================================
   function V_Lefebvre_V11(r)
      real(8), intent(in) :: r
      real(8)             :: V_Lefebvre_V11

      real(8), parameter :: Dm = 15000.0d0
      real(8), parameter :: re = 1.6d0
      real(8), parameter :: bm = 1.9685d0

      V_Lefebvre_V11 = Dm*(1.0d0 - exp(-bm*(r - re)))**2
   end function V_Lefebvre_V11

   !=========================================================
   ! Complex-valued Morse channel V11.
   !=========================================================
   function V_Lefebvre_V11_cmplx(r)
      complex(8), intent(in) :: r
      complex(8)             :: V_Lefebvre_V11_cmplx

      real(8), parameter :: Dm = 15000.0d0
      real(8), parameter :: re = 1.6d0
      real(8), parameter :: bm = 1.9685d0

      V_Lefebvre_V11_cmplx = Dm*(1.0d0 - exp(-bm*(r - re)))**2
   end function V_Lefebvre_V11_cmplx

   !=========================================================
   ! Real two-channel Lefebvre 2x2 potential matrix.
   !=========================================================
   function V_Lefebvre(r, k, m)
      real(8), intent(in) :: r
      integer, intent(in) :: k, m
      real(8)             :: V_Lefebvre

      real(8), parameter :: A  =  18154.95d0
      real(8), parameter :: B0 = -8000.0d0
      real(8), parameter :: ra =  2.48d0
      real(8), parameter :: gr =  2.2039d0
      real(8), parameter :: eps = 100.0d0

      if (k == 1 .and. m == 1) then
         V_Lefebvre = V_Lefebvre_V11(r)
      else if (k == 2 .and. m == 2) then
         V_Lefebvre = A*exp(-gr*(r - ra)) + B0
      else if ((k == 1 .and. m == 2) .or. (k == 2 .and. m == 1)) then
         V_Lefebvre = eps
      else
         V_Lefebvre = 0.0d0
      endif

   end function V_Lefebvre

   !=========================================================
   ! Complex version of Lefebvre 2x2 potential matrix.
   !=========================================================
   function V_Lefebvre_cmplx(r, k, m)
      complex(8), intent(in) :: r
      integer,    intent(in) :: k, m
      complex(8)             :: V_Lefebvre_cmplx

      real(8), parameter :: A   = 18154.95d0
      real(8), parameter :: B0  = -8000.0d0
      real(8), parameter :: ra  = 2.48d0
      real(8), parameter :: gr  = 2.2039d0
      real(8), parameter :: eps = 100.0d0

      if (k == 1 .and. m == 1) then
         V_Lefebvre_cmplx = V_Lefebvre_V11_cmplx(r)
      else if (k == 2 .and. m == 2) then
         V_Lefebvre_cmplx = A*exp(-gr*(r - ra)) + B0
      else if ((k == 1 .and. m == 2) .or. (k == 2 .and. m == 1)) then
         V_Lefebvre_cmplx = cmplx(eps, 0.0d0, 8)
      else
         V_Lefebvre_cmplx = (0.0d0, 0.0d0)
      endif

   end function V_Lefebvre_cmplx

   !=========================================================
   ! Two-channel model from:
   !   G. V. Sitnikov and O. I. Tolstikhin,
   !   Phys. Rev. A 67, 032714 (2003), Sec. B, Eq. (92a)
   !
   !   V11 = 15*exp(-0.5*r)
   !   V22 = 15*(r^2 - r - 1)*exp(-r) + 15
   !   V12 = V21 = 5*r*exp(-r)
   !
   ! Asymptotic thresholds:
   !   v1 = 0
   !   v2 = 15
   !=========================================================
   function V_Sitnikov_real(r, k, m)
      real(8), intent(in) :: r
      integer, intent(in) :: k, m
      real(8)             :: V_Sitnikov_real

      real(8) :: e1, e2

      e1 = exp(-0.5d0*r)
      e2 = exp(-r)

      if (k == 1 .and. m == 1) then
         V_Sitnikov_real = 15.0d0*e1
      else if (k == 2 .and. m == 2) then
         V_Sitnikov_real = 15.0d0*(r*r - r - 1.0d0)*e2 + 15.0d0
      else if ((k == 1 .and. m == 2) .or. (k == 2 .and. m == 1)) then
         V_Sitnikov_real = 5.0d0*r*e2
      else
         V_Sitnikov_real = 0.0d0
      endif

   end function V_Sitnikov_real

   !=========================================================
   ! Complex version of the Sitnikov-Tolstikhin 2x2 potential.
   !=========================================================
   function V_Sitnikov_cmplx(r, k, m)
      complex(8), intent(in) :: r
      integer,    intent(in) :: k, m
      complex(8)             :: V_Sitnikov_cmplx

      complex(8) :: e1, e2

      e1 = exp(-0.5d0*r)
      e2 = exp(-r)

      if (k == 1 .and. m == 1) then
         V_Sitnikov_cmplx = 15.0d0*e1
      else if (k == 2 .and. m == 2) then
         V_Sitnikov_cmplx = 15.0d0*(r*r - r - 1.0d0)*e2 + 15.0d0
      else if ((k == 1 .and. m == 2) .or. (k == 2 .and. m == 1)) then
         V_Sitnikov_cmplx = 5.0d0*r*e2
      else
         V_Sitnikov_cmplx = (0.0d0, 0.0d0)
      endif

   end function V_Sitnikov_cmplx

end module potentials   
