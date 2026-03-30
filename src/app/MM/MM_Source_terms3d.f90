module hydraMM_Source_Terms3D
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  use MOSE_Types_m
implicit none
private
public :: Assign_source

  procedure(source_if), pointer, public :: source

abstract interface
subroutine source_if(domain, i_rk)
  use MOSE_Types_m
  implicit none
  type(MOSE_domain_type) :: domain
  integer :: i_rk
  
end subroutine source_if
  end interface
contains

subroutine Assign_source(source_name)
  Implicit none
character(len=*), intent(inout) :: source_name

nullify(source)

select case (source_name)

case('FAKECHEM')
  source => Combustion3D

case default 

  source => null()
  


end select
end subroutine Assign_source


subroutine Combustion3D(domain, i_rk)
  use MOSE_Types_m
  use MOSE_Global_m
  use IO_RF,only: obj_RF
  use omp_lib
  use response_function
  implicit none
  type(MOSE_domain_type) :: domain
  integer :: i_rk, i

  call Combustion(domain)


  
end subroutine Combustion3D


subroutine Combustion (domain)
  use MOSE_Types_m
  use MOSE_Global_m
  use IO_RF,only: obj_RF
  use MOSKA_IO_Chemistry,only: obj_chem
  use FLINT_Lib_Thermodynamic
  use omp_lib
  !use response_function
  implicit none
  type(MOSE_domain_type) :: domain
  integer , parameter :: p2 = selected_real_kind(15,307) !Double Precision
  integer :: i,Bm,nbb,is,n1,n2,n3,j,k
  real(p2), parameter :: half = 0.5_p2
  real(p2), parameter :: one = 1_p2
  real(p2), parameter :: zero = 0_p2
  real(8)             :: hin_f
  real(8)             ::u,v,w,Mfvx,omegaMf,omo_ijk,omf1_ijk,omfr_ijk,omf_ijk,flux_diff,QQ
  real(8)             ::yfr,yi(nsc),om_base,omp_ijk,vel_mag,pi,x_i
  real(8)             ::rotot,temp,Rtot,mfr_build,mtot,avp
  real(8)             ::Tin_f, Uin_f


  
     do Bm = 1,domain%nb
        Tin_f = obj_RF(1)%Tin_f
        Uin_f = obj_RF(1)%Uin_f

        hin_f = f_tabT(Tin_f,2,h_tab)
        hin_f = hin_f + half*Uin_f*Uin_f
        
        
        n1=domain%Blk(Bm)%dim(1)
        n2=domain%Blk(Bm)%dim(2)
        n3=domain%Blk(Bm)%dim(3)

        !$omp do  &
        !$omp private(Mfvx,x_i,rotot,Rtot,temp,u) &
        !$omp private(vel_mag,yfr,yi,om_base,omf1_ijk,omfr_ijk,omf_ijk,omo_ijk,omp_ijk,flux_diff,i,QQ) 
        do k=1,n3
            do i=1,n1
                do j=1,n2

                
                  omo_ijk  = 0.d0
                  omf1_ijk = 0.d0
                  omfr_ijk = 0.d0
                  omf_ijk  = 0.d0
                  omp_ijk  = 0.d0
                  flux_diff = 0.d0
                  QQ = 0.d0



                  u = domain%Blk(Bm)%P(nu,i,j,k)
                  v = domain%Blk(Bm)%P(nv,i,j,k)
                  w = domain%Blk(Bm)%P(nw,i,j,k)
                  if (nprim>np) then
                  Mfvx = domain%Blk(Bm)%P(np+1,i,j,k)
                  else
                  Mfvx = 0.d0
                  endif
                  
                  x_i = (domain%Blk(Bm)%dl(i,j,k)%c(1)*half + domain%Blk(Bm)%node(i-1,j,k)%c(1))
                  call co_rotot_Rtot(domain%Blk(Bm)%P(1:nsc,i,j,k),rotot,Rtot)
                  temp = EOS(domain%Blk(Bm)%P(np,i,j,k),rotot,Rtot)
              
                    if (i.gt.1) then
                        vel_mag = sqrt(u*u + v*v + w*w)
                        yfr = Mfvx / (rotot+Mfvx)
                        if (yfr.lt.0) yfr = 0.d0
            
                        if(x_i.ge.obj_chem%xcomb) then
                          
                          do is=1,nsc
                            yi(is) = domain%Blk(Bm)%P(is,i,j,k)/(rotot+Mfvx) ! mass fractions
                            if (yi(is).lt.0) yi(is) = 0.d0
                          enddo
                          
                        
        
                          om_base = rotot*vel_mag/obj_chem%lcomb*yi(1)*exp(-obj_chem%Tref/temp)
                          
                
                          if (yi(1).lt.1.d-6) om_base = 0.d0
                
                            omf1_ijk = om_base*yi(2)
                            omfr_ijk = om_base*yfr

                            if (yi(2).lt.1.d-6) omf1_ijk = 0.d0
                            if (yfr.lt.1.d-6) omfr_ijk = 0.d0
                
                            omf_ijk  = omf1_ijk+omfr_ijk 
                            omo_ijk  = omf_ijk*obj_chem%o2f
                            omp_ijk  = omf_ijk + omo_ijk
                            
                            domain%Blk(Bm) % r(1,i,j,k) = domain%Blk(Bm) % r(1,i,j,k) +omo_ijk*domain%Blk(Bm)%vol(i,j,k)                       ! oxidizer
                            domain%Blk(Bm) % r(2,i,j,k) = domain%Blk(Bm) % r(2,i,j,k) +omf1_ijk*domain%Blk(Bm)%vol(i,j,k)                      ! fuel
                            domain%Blk(Bm) % r(3,i,j,k) = domain%Blk(Bm) % r(3,i,j,k) -omp_ijk*domain%Blk(Bm)%vol(i,j,k)                       ! combustion products
                            domain%Blk(Bm) % r(np+1,i,j,k) = domain%Blk(Bm) % r(np+1,i,j,k) +omfr_ijk*domain%Blk(Bm)%vol(i,j,k)
                          
                    endif
                
                    flux_diff = rotot*vel_mag*yfr*exp(obj_chem%Dm*(x_i/obj_RF(1)%chmb_len))
                              
              
                    if (yfr.lt.1.d-6) flux_diff = 0.d0
                  
                    
                        domain%Blk(Bm) % r(2,i,j,k) = domain%Blk(Bm) % r(2,i,j,k) - flux_diff*domain%Blk(Bm)%vol(i,j,k) 
                        domain%Blk(Bm) % r(np+1,i,j,k) = domain%Blk(Bm) % r(np+1,i,j,k) + flux_diff*domain%Blk(Bm)%vol(i,j,k)
                        
                        
                        
                        QQ =  (-omo_ijk*obj_chem%dhf_o -omf_ijk*obj_chem%dhf_f +omp_ijk*obj_chem%dhf_p)* domain%Blk(Bm)%vol(i,j,k)
                        
                        domain%Blk(Bm) % r(np,i,j,k) = domain%Blk(Bm) % r(np,i,j,k) + QQ + (-omfr_ijk-flux_diff)*hin_f*domain%Blk(Bm)%vol(i,j,k)
                            
                        domain%Blk(Bm) % r(nu,i,j,k) = domain%Blk(Bm) % r(nu,i,j,k) + (-omfr_ijk-flux_diff)*Uin_f*domain%Blk(Bm)%vol(i,j,k)
                    

                     endif

                enddo
            enddo
        enddo

          !$omp end do
    enddo

  !Fake combustion--------------------------------------------------------------------------------

end subroutine combustion


end module hydraMM_Source_Terms3D
