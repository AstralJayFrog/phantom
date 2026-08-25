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
! :Dependencies: dim, eos, infile_utils, io, part, timestep, units
!
 implicit none

 real :: ds,mus,mu2,I0,rhos
 character(len=12) :: ds_in, rho_ref 

 public :: granularshearfunc,set_defaults_granular
 public :: write_options_granular, read_options_granular

 private

contains

!----------------------------------------------------------------
!+
!  set default values for the viscosity coefficients
!+
!----------------------------------------------------------------
subroutine set_defaults_granular
   use units, only:in_code_units

   integer :: ierr

   ! Units that don't call 'in_code_units' are unitless.
   ! The importance of these values is taken from Bui 2021 (eqn 64 and 65)
   ! The default values are for glass balls, taken from Jop et al., 2006
   
   ! Original values from Forterre 2003
   mus = 0.38 ! static friction coefficient, unitless
   mu2 = 0.64 ! critical friction angle at high I, unitless
   I0 = 0.279 ! a material constant (unnamed?), unitless

   ds_in = '0.053*cm' ! grain diameter   
   ds = in_code_units(ds_in, ierr, 'length')
   !if (ierr /= 0) then
     !call error(label,'could not convert units of length for ds')
     !nerr = nerr + 1
   !endif  
   ierr = 0
   
   rho_ref = '2.5*g/cm^3' ! solid density of the material
   rhos = in_code_units(rho_ref, ierr, 'density')
   !if (ierr /= 0) then
     !call error(label,'could not convert units of length for rhos')
     !nerr = nerr + 1
   !endif  
   ierr = 0

end subroutine set_defaults_granular

!----------------------------------------------------------------
!+
!  function which returns shear parameter \nu as per
!  the mu I model as explained in (Bui 2021)
!+
!----------------------------------------------------------------
real function granularshearfunc(strain, granularpressure, rho1i)

 real, intent(in) :: strain(6), granularpressure, rho1i
 real :: I, muI, epsilondotcolon

 ! -- mu I model
 ! Go through how shear viscosity should be calculated in the mu I model
 ! NOTE this gives eta/rho (a value called nu in Price 2018)

 epsilondotcolon = strain(1)**2 + 2 * strain(2)**2 + 2 * strain(3)**2 + strain(4)**2 + 2*strain(5)**2 + strain(6)**2

 if(granularpressure > 0) then
   I = ds * (2 * epsilondotcolon) ** (0.5) / (granularpressure/rhos)**(0.5)

   ! If I is very low or zero, then muI --> mus
   if(I == 0) then
     muI = mus
   else
     muI = mus + (mu2 - mus) / ((I0/I) + 1)
   endif
   
   ! Ensure no division by zero if the strain rate is zero
   if(epsilondotcolon /= 0) then
    granularshearfunc = muI * granularpressure / (2 * epsilondotcolon) ** (0.5) * rho1i
   else
    granularshearfunc = 0
   endif

 else
   ! If pressure is negative, then the shear stress should be zero as per Bui 2021
   granularshearfunc = 0
 endif

end function granularshearfunc

!----------------------------------------------------------------
!+
!  routine to write physical granularity options to input file
!+
!----------------------------------------------------------------
subroutine write_options_granular(iwritein)
  use infile_utils, only:write_inopt
  integer, intent(in) :: iwritein
 
  write(iwritein,"(/,a)") '# options controlling granular flow model'
  call write_inopt(ds_in,'ds','grain diameter of the material (code units or e.g. 0.053 cm)',iwritein)
  call write_inopt(mus,'mus','static friction coefficient (e.g 0.38), unitless',iwritein)
  call write_inopt(mu2,'mu2','material constant for yielding (e.g 0.64), unitless',iwritein)
  call write_inopt(I0,'I0','material constant for shear (e.g 0.279), unitless',iwritein)
  call write_inopt(rhos,'rhos','solid density of the material (code units or e.g 2.5 g/cm^3)',iwritein)
 
 end subroutine write_options_granular
 
 !----------------------------------------------------------------
 !+
 !  routine to read physical granularity options from input file
 !+
 !----------------------------------------------------------------
 subroutine read_options_granular(db,nerr)
  use io,           only:error
  use infile_utils, only:inopts,read_inopt
  use units,        only:in_code_units
  type(inopts), intent(inout) :: db(:)
  integer,      intent(inout) :: nerr
  character(len=*), parameter :: label = 'read_infile'
  integer :: ierr = 0
 
  call read_inopt(ds_in,'ds',db,errcount=nerr)
  call read_inopt(rho_ref,'rhos',db,errcount=nerr)

  call read_inopt(mu2,'mu2',db,errcount=nerr,min=0.,max=1.,default=mu2)
  call read_inopt(I0,'I0',db,errcount=nerr,min=0.,max=1.,default=I0)
  call read_inopt(mus,'mus',db,errcount=nerr,min=0.,max=1.,default=mus)
 
  ds = in_code_units(ds_in, ierr, 'length')
  if (ierr /= 0) then
     call error(label,'could not convert units of length for ds')
     nerr = nerr + 1
  endif  
  ierr = 0

  rhos = in_code_units(rho_ref, ierr, 'density')
  if (ierr /= 0) then
     call error(label,'could not convert units of density for rhos')
     nerr = nerr + 1
  endif  

 end subroutine read_options_granular


end module granular
