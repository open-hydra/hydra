module hydraMQ2D_Mod_Explicit
    use iso_fortran_env, only: I4 => int32, R8 => real64  
    implicit none

    private
    public :: Explicit_Step

contains   
    subroutine Explicit_Step ( Q2D_domain, MOSKA_domain, Q2D_config, MOSKA_config )
    ! MOSKA modules
    use MOSKA_Types_m, only:  MOSKA_config_type => config_type, MOSKA_1Ddomain_type
    use MOSKA_Global_m, only: MOSKAnres => nres
    use MOSKA_Mod_dt, only: MOSKA_Compute_dt => Compute_dt1D
    use MOSKA_Lib_Extrapolation
    use MOSKA_Mod_Fluxes, only: MOSKA_Fluxes => Fluxes1D
    use MOSKA_Mod_BC
    use MOSKA_Lib_Newstate
    use MOSKA_Lib_Chemistry, only: MOSKA_Chemistry_Newstate => Chemistry_Newstate
    use MOSKA_Mod_Diagnostic, only: MOSKA_Compute_Residual => Compute_Residual
    use Source_Terms
    ! Q2D modules
    use Q2D_Types_m,  only: Q2D_config_type => config_type , Q2D_domain_type
    use Q2D_Global_m, only: Q2Dnres => nres
    use Q2D_Mod_dt, only: Q2D_Compute_dt => Compute_dt
    use Q2D_Lib_Extrapolation
    use Q2D_Mod_Fluxes,  only: Q2D_Fluxes => Fluxes
    use Q2D_Mod_BC
    use Q2D_Lib_RANS
    use Q2D_Lib_Newstate
    use Q2D_Lib_Chemistry, only: Q2D_Chemistry_Newstate => Chemistry_Newstate
    use Q2D_Mod_Diagnostic, only: Q2D_Compute_Residual => Compute_Residual
    use Q2D_IO_Solution
    use Q2D_Lib_Strang, only: N_strang
    use Q2D_Lib_RK, only: obj_rk
    use Q2D_Mod_Sources
    ! HYDRA modules
    use hydra_Types_m
    use hydraMQ2D_Lib_Extrapolation, only: Compute_Extrapolation_Coupled 
    use hydraMQ2D_Mod_Coupling, only: Set_Global_dt, BC_Coupled_Fluxes

    implicit none
    type(Q2D_domain_type), intent(inout)     :: Q2D_domain
    type(MOSKA_1Ddomain_type), intent(inout) :: MOSKA_domain
    type(Q2D_config_type), intent(inout)     :: Q2D_config
    type(MOSKA_config_type), intent(inout)   :: MOSKA_config
    ! Local
    logical  :: endsim, iosim
    integer  :: i_rk, i_strang, b
    real(R8) :: averageQ2D(Q2Dnres), averageMOSKA(MOSKAnres)
    real(R8)  :: t_start, t_end

    Q2D_config % change_MG = .false.
    Q2D_domain % iter = Q2D_domain % iter + 1
    MOSKA_domain % iter = MOSKA_domain % iter + 1

    Q2D_config % iter_from_call = Q2D_config % iter_from_call + 1
    MOSKA_config % iter_from_call = MOSKA_config % iter_from_call + 1

    Q2D_config % iter_general = Q2D_config % iter_general + 1
    MOSKA_config % iter_general = MOSKA_config % iter_general + 1

    ! Controllo time-accurate a livello hydra
    if (Q2D_config%time_accurate.and.MOSKA_config%time_accurate) then
        call Set_Global_dt ( Q2D_domain, MOSKA_domain, Q2D_config, MOSKA_config )
        Q2D_domain % time = Q2D_domain % time + Q2D_domain % dtglobal
        MOSKA_domain % time = MOSKA_domain % time + MOSKA_domain % dtglobal
    else
        Q2D_domain % dtglobal = 1d5
        MOSKA_domain % dtglobal = 1d5
        call Q2D_Compute_dt ( Q2D_domain, Q2D_config )        ! Compute local and minimum time step
        call MOSKA_Compute_dt ( MOSKA_domain, MOSKA_config )  ! Compute local and minimum time step			
    endif

		!$omp parallel
        strang: do i_strang = 1, N_strang ! Strang splitting: L = L_rk(dt/2) * L_chem(dt) * L_rk(dt/2)

            call Copy_State ( Q2D_domain ) 
            call Copy_State1D ( MOSKA_domain )

            rk: do i_rk = 1, obj_rk%n_RK

                call Fill_Ghost_Cell ( Q2D_domain )
                call Fill_Ghost_Cell1D ( MOSKA_domain )
                call Compute_Extrapolation_Coupled ( Q2D_domain, MOSKA_domain )

                call Q2D_Fluxes(Q2D_domain)
                call MOSKA_Fluxes(MOSKA_domain)

                call BC_Coupled_Fluxes ( Q2D_domain, MOSKA_domain )
                call BC_Fluxes(Q2D_domain)                
                call BC_Fluxes1D(MOSKA_domain)

                call source1D ( MOSKA_domain, i_rk, i_strang, MOSKA_domain % iter)  ! ( source terms) 1D

                !if (associated(source3d)) then
                !	call source3D ( Q2D_domain, i_rk)      ! (source terms) 3D
                !else 
                !!skip
                !endif

                if (RANS .or. SOURCE_TERM) call Extrapolate_Wall(Q2D_domain)
                if (RANS)        call RANS_Source_Terms(Q2D_domain)
                if (SOURCE_TERM) call Compute_Sources(Q2D_domain)

                call RK_Newstate ( Q2D_domain, i_rk )             ! State update 3D
                call RK_Newstate1D ( MOSKA_domain, i_rk )          ! State update 1D

            enddo rk

            ! Chemistry, only if i_strang=1 and N_strang=2
            if ( i_strang + N_strang == 3 ) then
                !No combustion in 1D (Injector Hypothesis)
                !  do b = 1, MOSKA_domain % nb						
                !    call MOSKA_Chemistry_Newstate ( MOSKA_domain % blk(b) % P,MOSKA_domain % blk(b) % dtlocal, MOSKA_domain % blk(b) % dim )						
                !  enddo				
                do b = 1, Q2D_domain % nb
                    if (Q2D_domain%blk(b)%no_chem) cycle
                    call Q2D_Chemistry_Newstate  ( Q2D_domain % blk(b) % P, Q2D_domain  % blk(b) % dtlocal, Q2D_domain  % blk(b) % dim )
                enddo
            endif

        enddo strang
		!$omp end parallel

        ! Save residuals
        if ( (Q2D_config % iter_from_call == 1) .or. (mod (Q2D_domain % iter, obj_io%res_diter) == 0) ) then
            Q2D_config % residuotot = 0.d0
            MOSKA_config % residuotot = 0.d0
            do b = 1, Q2D_domain % nb
                call Q2D_Compute_Residual ( new=Q2D_domain%blk(b)%p, &
                                            old=Q2D_domain%blk(b)%po, &
                                            dt=Q2D_domain%blk(b)%dtlocal, &
                                            n=Q2D_domain%blk(b)%dim, &
                                            average=averageQ2D, &
                                            total=Q2D_config%residuotot )
            enddo
            do b = 1, MOSKA_domain % nb
                call MOSKA_Compute_Residual ( new=MOSKA_domain%blk(b)%p, &
                                              old=MOSKA_domain%blk(b)%po, &
                                              dt=MOSKA_domain%blk(b)%dtlocal, &
                                              n=MOSKA_domain%blk(b)%dim, &
                                              average=averageMOSKA, &
                                              total=MOSKA_config%residuotot )
            enddo
            MOSKA_config % residuotot = sqrt ( MOSKA_config % residuotot ) ! L2 norm time derivative
            Q2D_config % residuotot = sqrt ( Q2D_config % residuotot ) ! L2 norm time derivative
        endif

        endsim = ( Q2D_config % iter_from_call >= Q2D_config % itermax ) &
                .or. ( Q2D_config % residuotot(1) <= Q2D_config % restop ) &
                .or. ( Q2D_domain % time >= Q2D_config % timemax )
        iosim  = ( mod (Q2D_domain % iter, obj_io%sol_diter) == 0) &
                .or. ( mod (Q2D_domain % iter, obj_io%bck_diter) == 0) &
                .or. ( Q2D_domain % time >= Q2D_config % time_from_call + obj_io%sol_dtime )

        if ( endsim ) then
            Q2D_config % TODO = 3
            MOSKA_config % TODO = 3
        elseif ( iosim ) then
            Q2D_config % TODO = 2
            MOSKA_config % TODO = 2
        else
            Q2D_config % TODO = 1
            MOSKA_config % TODO = 1
        endif

    end subroutine Explicit_Step

    
    subroutine Copy_State1D ( domain )
        use MOSKA_Types_m

        implicit none
        type(MOSKA_1Ddomain_type), intent(inout) :: domain
        ! Local
        integer :: i, b

        do b = 1, domain % nb
			!$omp do 
            do i = 1, domain % blk(b) % dim

                domain % blk(b) % PO (:,i) = domain % blk(b) % P (:,i)

            enddo; 
			!$omp end do
        enddo

    end subroutine Copy_State1D

    subroutine Copy_State ( domain )
        use Q2D_Types_m

        implicit none
        type(Q2D_domain_type), intent(inout) :: domain
        ! Local
        integer :: i, j, k, b

        do b = 1, domain % nb
			!$omp do collapse(2)
            do j = 1, domain % blk(b) % dim(2)
            do i = 1, domain % blk(b) % dim(1)
    
                domain % blk(b) % PO (:,i,j) = domain % blk(b) % P (:,i,j)

            enddo; enddo
			!$omp end do
        enddo

    end subroutine Copy_State

end module hydraMQ2D_Mod_Explicit

