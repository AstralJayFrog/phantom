!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module granular
!
! Routines related to the μ(I) model for granular viscosity
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - 
!
! :Dependencies: dim, eos, infile_utils, io, part, timestep
!
 implicit none

 real :: ds,mus,mu2,I0,rhos

 public :: granularshearfunc,set_defaults_granular

 private

contains

!----------------------------------------------------------------
!+
!  set default values for the viscosity coefficients
!+
!----------------------------------------------------------------
subroutine set_defaults_granular
   use units, only:in_code_units

   integer :: gerr

   ! All units that don't call 'in_code_units' are unitless.
   ds = in_code_units('0.53.*mm', gerr, 'length')
   mus = 0.38
   mu2 = 0.64 
   I0 = 0.279
   rhos = in_code_units('2.5.*g/cm^3', gerr, 'density')

   write(*,*) 'ds =', ds, 'rhos =', rhos, 'mus =', mus, 'mu2 =', mu2, 'I0 =', I0

end subroutine set_defaults_granular

!----------------------------------------------------------------
!+
!  function which returns shear parameter \nu as per
!  the mu I model as explained in (Bui 2021)
!+
!----------------------------------------------------------------
real function granularshearfunc(strain, granularpressure)

 real, intent(in) :: strain(6), granularpressure
 real :: I, muI, epsilondotcolon
 integer :: derr

 ! -- mu I model
 ! Go through how shear viscosity should be calculated in the mu I model

 epsilondotcolon = strain(1)**2 + 2 * strain(2)**2 + 2 * strain(3)**2 + strain(4)**2 + 2*strain(5)**2 + strain(6)**2

 if(granularpressure > 0) then
   I = ds * (2 * epsilondotcolon) ** (0.5) / (granularpressure/rhos)**(0.5)
 else
   I = 1.0
 endif

! Check if I is zero (if so, make sure there is no division by zero)
 if(I == 0) then
   muI = mus
 else
   muI = mus + (mu2 - mus) / ((I0/I) + 1)
 endif

 if(epsilondotcolon /= 0) then
   granularshearfunc = muI * granularpressure / (2 * epsilondotcolon) ** (0.5)
 else
   granularshearfunc = 0
 endif

end function granularshearfunc

end module granular
