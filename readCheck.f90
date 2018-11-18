module rstReader

   implicit none
   
   integer :: version 
   integer :: n_r_max,n_theta_max,n_phi_tot,minc,nalias,n_r_ic_max,n_time_step
   integer :: rscheme_oc_n_max, rscheme_oc_order_boundary
   real(kind=8) :: alph1, alph2, fd_stretch, fd_ratio
   real(kind=8) :: time, dt, dtNew
   real(kind=8) :: ra, pr, raxi, sc, prmag, ek, radratio, sigma_ratio
   real(kind=8) :: lorentz_torque_ic,lorentz_torque_ma
   real(kind=8) :: omega_ic1,omegaOsz_ic1,tOmega_ic1
   real(kind=8) :: omega_ic2,omegaOsz_ic2,tOmega_ic2
   real(kind=8) :: omega_ma1,omegaOsz_ma1,tOmega_ma1
   real(kind=8) :: omega_ma2,omegaOsz_ma2,tOmega_ma2
   
   logical :: l_heat, l_chemical_conv, l_mag, l_cond_ic
   
   complex(kind=8), allocatable :: wlm(:,:), zlm(:,:), dwdt(:,:), dzdt(:,:)
   complex(kind=8), allocatable :: plm(:,:), slm(:,:), dpdt(:,:), dsdt(:,:)
   complex(kind=8), allocatable :: xilm(:,:), dxidt(:,:)
   complex(kind=8), allocatable :: blm(:,:), jlm(:,:), dbdt(:,:), djdt(:,:)
   complex(kind=8), allocatable :: blm_ic(:,:), jlm_ic(:,:)
   complex(kind=8), allocatable :: dbdt_ic(:,:),djdt_ic(:,:)
 

   character(len=4) :: rscheme_oc_version


contains

   subroutine readRst(filename, endian, l_cheb, lm_max)

      !-- Input variables
      character(len=*), intent(in) :: filename
      character(len=1), intent(in) :: endian
      integer, intent(in) :: lm_max
      logical, intent(in) :: l_cheb

      if ( endian == 'B' ) then
         open(unit=10, file=filename, form='unformatted', convert='big_endian')
      else
         open(unit=10, file=filename, form='unformatted', convert='little_endian')
      end if

      !-- Header
      read(10) version
      read(10) time,dt,n_time_step
      read(10) ra,pr,raxi,sc,prmag,ek,radratio,sigma_ratio
      read(10) n_r_max,n_theta_max,n_phi_tot,minc,nalias,n_r_ic_max
      
      if(l_cheb) then
         read(10) rscheme_oc_version, rscheme_oc_n_max,  &
         &        rscheme_oc_order_boundary, alph1, alph2
      else
         read(10) rscheme_oc_version, rscheme_oc_n_max,  &
         &        rscheme_oc_order_boundary, fd_stretch, fd_ratio
      end if

      read(10) lorentz_torque_ic,lorentz_torque_ma,         &
      &        omega_ic1,omegaOsz_ic1,tOmega_ic1,           &
      &        omega_ic2,omegaOsz_ic2,tOmega_ic2,           &
      &        omega_ma1,omegaOsz_ma1,tOmega_ma1,           &
      &        omega_ma2,omegaOsz_ma2,tOmega_ma2,           &
      &        dtNew


      !-- Logicals
      read(10) l_heat, l_chemical_conv, l_mag, l_cond_ic

      
      !-- Allocate potential arrays, take lm_max as input
      allocate(wlm(lm_max,n_r_max), zlm(lm_max,n_r_max))
      allocate(dwdt(lm_max,n_r_max), dzdt(lm_max,n_r_max))
      allocate(plm(lm_max,n_r_max), dpdt(lm_max,n_r_max))
      if (l_heat) allocate(slm(lm_max,n_r_max), dsdt(lm_max,n_r_max))
      if (l_chemical_conv) allocate(xilm(lm_max,n_r_max), dxidt(lm_max,n_r_max))
      if (l_mag) then
         allocate(blm(lm_max,n_r_max), jlm(lm_max,n_r_max))
         allocate(dbdt(lm_max,n_r_max), djdt(lm_max,n_r_max))
         if (l_cond_ic) then
            allocate(blm_ic(lm_max,n_r_ic_max), jlm_ic(lm_max,n_r_ic_max))
            allocate(dbdt_ic(lm_max,n_r_ic_max), djdt_ic(lm_max,n_r_ic_max))
         end if
      end if


      ! Read potential arrays
      read(10) wlm
      read(10) dwdt
      read(10) zlm
      read(10) dzdt
      read(10) plm
      read(10) dpdt
      if (l_heat) then
         read(10) slm
         read(10) dsdt
      end if
      if (l_chemical_conv) then
         read(10) xilm
         read(10) dxidt
      end if
      if (l_mag) then
         read(10) blm
         read(10) dbdt
         read(10) jlm
         read(10) djdt
         if (l_cond_ic) then
            read(10) blm_ic
            read(10) dbdt_ic
            read(10) jlm_ic
            read(10) djdt_ic
         end if
      end if


      close(10)


   end subroutine readRst

end module rstReader
