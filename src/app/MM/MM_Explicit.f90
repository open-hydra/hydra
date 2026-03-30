module hydraMM_Mod_Explicit
    use iso_fortran_env, only: I4 => int32, R8 => real64
  
    implicit none
  
  contains
  
  
    subroutine Explicit_Hybrid (  MOSE_domain,MOSKA_domain, MOSEconfig, MOSKAconfig )
    use MOSKA_Types_m,  only:  MOSKA_config_type => config_type, MOSKA_1Ddomain_type
    use MOSKA_Lib_Extrapolation
    use MOSKA_Mod_Fluxes, only: MOSKA_Fluxes => Fluxes1D
    use MOSKA_Mod_BC
    use MOSE_Global_m, only: MOSKAnres => nres
    use MOSKA_Lib_Newstate
    use MOSKA_Lib_Chemistry, only: MOSKA_Chemistry_Newstate => Chemistry_Newstate
    use MOSKA_Mod_Diagnostic, only: MOSKA_Compute_Residual => Compute_Residual
    use MOSKA_Mod_dt, only: MOSKA_Compute_dt => Compute_dt1D
    use MOSE_Mod_dt, only: MOSE_Compute_dt => Compute_dt
    use Source_Terms
    use MOSE_Types_m,  only: MOSE_config_type => config_type , MOSE_domain_type
    use MOSE_Global_m, only: MOSEnres => nres
    use MOSE_Lib_Strang, only: N_strang
    use MOSE_Lib_RK, only: obj_rk
    use MOSE_Lib_Extrapolation
    use MOSE_Mod_Fluxes,  only: MOSE_Fluxes => Fluxes
    use MOSE_Mod_BC
    use MOSE_Lib_Newstate
    use MOSE_Lib_Chemistry, only: MOSE_Chemistry_Newstate => Chemistry_Newstate
    use MOSE_Mod_Diagnostic, only: MOSE_Compute_Residual => Compute_Residual
    use MOSE_IO_Solution
    use hydraMM_Mod_Coupling
    use hydra_Types_m
    use hydraMM_Source_Terms3D, only: source3D => source

   
      implicit none
      type(MOSE_domain_type), intent(inout) :: MOSE_domain
      type(MOSKA_1Ddomain_type), intent(inout) :: MOSKA_domain
      type(MOSE_config_type), intent(inout)      :: MOSEconfig
      type(MOSKA_config_type), intent(inout)      :: MOSKAconfig
      ! Local
      logical  :: endsim, iosim
      integer  :: i_rk, i_strang, b
      real(R8) :: averageMOSE(MOSEnres),averageMOSKA(MOSKAnres)
      real(8) :: t_start,t_end
     
     
      MOSE_domain % iter = MOSE_domain % iter + 1
      MOSKA_domain % iter = MOSKA_domain % iter + 1

      MOSEconfig % iter_from_call = MOSEconfig % iter_from_call + 1
      MOSKAconfig % iter_from_call = MOSKAconfig % iter_from_call + 1

      MOSEconfig % iter_general   = MOSEconfig % iter_general + 1
      MOSKAconfig % iter_general  = MOSKAconfig % iter_general + 1
  
      if (MOSKAconfig%time_accurate .OR. MOSEconfig%time_accurate) then

        call Fix_Global_dt ( MOSKA_domain, MOSE_domain, MOSKAconfig, MOSEconfig )
        MOSE_domain % time = MOSE_domain % time + MOSE_domain % dtglobal
        MOSKA_domain % time = MOSKA_domain % time + MOSKA_domain % dtglobal
   
      else
        MOSE_domain % dtglobal = 1d5
        MOSKA_domain % dtglobal = 1d5
        call MOSE_Compute_dt ( MOSE_domain, MOSEconfig )  ! Compute local and minimum time step
        call MOSKA_Compute_dt ( MOSKA_domain, MOSKAconfig )  ! Compute local and minimum time step
       
      endif
     
      !$omp parallel
      strang: do i_strang = 1, N_strang ! Strang splitting: L = L_rk(dt/2) * L_chem(dt) * L_rk(dt/2)
       
        call Copy_State ( MOSE_domain )
       
       
       
        call Copy_State1D ( MOSKA_domain )
         
       
        
        rk: do i_rk = 1, obj_rk%n_RK
          
          call Fill_Ghost_Cell1D ( MOSKA_domain )
          call Fill_Ghost_Cell ( MOSE_domain )
          
         
          call Fill_MOSExMOSKA  ( MOSKA_domain, MOSE_domain )
          call MOSKA_Fluxes(MOSKA_domain)
         
        
         
          call MOSE_Fluxes(MOSE_domain)
          call BC_Fluxes1D(MOSKA_domain)
            
         
          call BC_Fluxes(MOSE_domain)                
          call Flux_MOSExMOSKA ( MOSKA_domain, MOSE_domain ) !Filling Ghost and Computing fluxes cells with 1D-3D connection
               
          

          call source1D ( MOSKA_domain, i_rk, i_strang, MOSKA_domain % iter)       ! ( source terms) 1D

          if (associated(source3d)) then
            
            call source3D ( MOSE_domain, i_rk)      ! (source terms) 3D

          else 
          !skip
          endif

          
         
          call RK_Newstate1D ( MOSKA_domain, i_rk )          ! State update 1D
          call RK_Newstate ( MOSE_domain, i_rk )             ! State update 3D
       
  
          
          
        enddo rk
        
        ! Chemistry, only if i_strang=1 and N_strang=2
        if ( i_strang + N_strang == 3 ) then
          !No combustion in 1D (Injector Hypothesis)
          !  do b = 1, MOSKA_domain % nb
              
          !    call MOSKA_Chemistry_Newstate ( MOSKA_domain % blk(b) % P,MOSKA_domain % blk(b) % dtlocal, MOSKA_domain % blk(b) % dim )
             
          !  enddo
         
          do b = 1, MOSE_domain % nb
            
              call MOSE_Chemistry_Newstate  ( MOSE_domain % blk(b) % P, MOSE_domain  % blk(b) % dtlocal, MOSE_domain  % blk(b) % dim )
            
            enddo

        endif
  
      enddo strang
      !$omp end parallel
  
      ! Save residuals
      if ( (MOSEconfig % iter_from_call == 1) .or. (mod (MOSE_domain % iter, obj_io%res_diter) == 0) ) then
        MOSEconfig % residuotot = 0.d0
        MOSKAconfig % residuotot = 0.d0
        do b = 1, MOSE_domain % nb
          
          call MOSE_Compute_Residual ( new=MOSE_domain%blk(b)%p, &
                                  old=MOSE_domain%blk(b)%po, &
                                  dt=MOSE_domain%blk(b)%dtlocal, &
                                  n=MOSE_domain%blk(b)%dim, &
                                  average=averageMOSE, &
                                  total=MOSEconfig%residuotot )
                                  
        enddo


         do b = 1, MOSKA_domain % nb
          
          call MOSKA_Compute_Residual ( new=MOSKA_domain%blk(b)%p, &
                                  old=MOSKA_domain%blk(b)%po, &
                                  dt=MOSKA_domain%blk(b)%dtlocal, &
                                  n=MOSKA_domain%blk(b)%dim, &
                                  average=averageMOSKA, &
                                  total=MOSKAconfig%residuotot )
                                  
        enddo

        MOSKAconfig % residuotot = sqrt ( MOSKAconfig % residuotot ) ! L2 norm time derivative
        MOSEconfig % residuotot = sqrt ( MOSEconfig % residuotot ) ! L2 norm time derivative
      endif
      
      endsim = ( MOSEconfig % iter_from_call >= MOSEconfig % itermax ) &
          .or. ( MOSEconfig % residuotot(1) <= MOSEconfig % restop ) &
          .or. ( MOSE_domain % time >= MOSEconfig % timemax )
      iosim  = ( mod (MOSE_domain % iter, obj_io%sol_diter) == 0) &
          .or. ( mod (MOSE_domain % iter, obj_io%bck_diter) == 0) &
          .or. ( MOSE_domain % time >= MOSEconfig % time_from_call + obj_io%sol_dtime )
  
      if ( endsim ) then
        MOSEconfig % TODO = 3
        MOSKAconfig % TODO = 3
      elseif ( iosim ) then
        MOSEconfig % TODO = 2
        MOSKAconfig % TODO = 2
      else
        MOSEconfig % TODO = 1
        MOSKAconfig % TODO = 1
      endif
        
    end subroutine Explicit_Hybrid
  
  
    subroutine Copy_State1D ( domain )
      use MOSKA_Types_m
  
      implicit none
      type(MOSKA_1Ddomain_type), intent(inout) :: domain
      ! Local
      integer :: i, b
  
      do b = 1, domain % nb
        !$omp do collapse(1)
       
        do i = 1, domain % blk(b) % dim
          
          domain % blk(b) % PO (:,i) = domain % blk(b) % P (:,i)
        
        enddo; 
        !$omp end do
      enddo
  
    end subroutine Copy_State1D

    subroutine Copy_State ( domain )
      use MOSE_Types_m

      implicit none
      type(MOSE_domain_type), intent(inout) :: domain
      ! Local
      integer :: i, j, k, b

      do b = 1, domain % nb
        !$omp do collapse(3)
        do k = 1, domain % blk(b) % dim(3)
        do j = 1, domain % blk(b) % dim(2)
        do i = 1, domain % blk(b) % dim(1)
          
          domain % blk(b) % PO (:,i,j,k) = domain % blk(b) % P (:,i,j,k)
        
        enddo; enddo; enddo
        !$omp end do
      enddo

   end subroutine Copy_State
  
  
  end module hydraMM_Mod_Explicit
  
