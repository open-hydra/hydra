module hydraMF_Lib_BC_Coupled_Fluxes
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Compute_FUSS_Flux, Compute_MOSE_Flux
  
contains

  subroutine Compute_FUSS_Flux ( FUSS_domain, MOSE_domain )
    use hydraMF_Lib_BC_Coupled_Wall_Temperature, only: FUSS_Wall_Temperature
    use FUSS_Advanced_Types_m
    use MOSE_Advanced_Types_m
    use MOSE_Global_m, only: np
    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    ! Local
    integer  :: b, f, lower, upper, i
    integer  :: Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs
    integer  :: b_, f_, lower_, upper_, i_
    integer  :: Bm_, Im_, Jm_, Km_, Fm_, Bs_, Is_, Js_, Ks_, Fs_
    real(R8) :: FUSS_new_flux

    !$omp parallel
    upper = 0
    do b = 1, FUSS_domain%nb
      do f = 1, 6
        lower = upper + 1
        upper = upper + FUSS_domain%n_bf(b,f)
        !$omp do schedule ( dynamic ) private (i,Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs) &
        !$omp private (b_, f_, lower_, upper_, i_) &
        !$omp private (Bm_, Im_, Jm_, Km_, Fm_, Bs_, Is_, Js_, Ks_, Fs_) &
        !$omp private (FUSS_new_flux)
        do i = lower, upper
          select case ( FUSS_domain%bc(i)%type )
            case (1000)
              Bm  = FUSS_domain % bc(i) % b
              Im  = FUSS_domain % bc(i) % i
              Jm  = FUSS_domain % bc(i) % j 
              Km  = FUSS_domain % bc(i) % k 
              Fm  = FUSS_domain % bc(i) % f
              Bs  = FUSS_domain % bc(i) % bs
              Is  = FUSS_domain % bc(i) % is
              Js  = FUSS_domain % bc(i) % js
              Ks  = FUSS_domain % bc(i) % ks
              Fs  = FUSS_domain % bc(i) % fs
              call FUSS_Wall_Temperature( Im, Jm, Km, Fm, FUSS_domain % blk(Bm), &
                                          Is, Js, Ks, Fs, MOSE_domain % blk(Bs), &
                                          FUSS_domain%bc(i)%Mg(1), FUSS_domain%bc(i)%Tg, FUSS_new_flux )
              FUSS_domain % bc(i) % ext_flux = FUSS_domain % bc(i) % ext_flux + 0.5 * FUSS_new_flux
              
              upper_ = 0
              do b_ = 1, MOSE_domain%nb
                do f_ = 1, 6
                  lower_ = upper_ + 1
                  upper_ = upper_ + MOSE_domain%n_bf(b_,f_)
                  do i_ = lower_, upper_
                    select case ( MOSE_domain%bc(i_)%type )
                      case (1000)
                        Bm_  = MOSE_domain % bc(i_) % b
                        Im_  = MOSE_domain % bc(i_) % i
                        Jm_  = MOSE_domain % bc(i_) % j 
                        Km_  = MOSE_domain % bc(i_) % k 
                        Fm_  = MOSE_domain % bc(i_) % f
                        Bs_  = MOSE_domain % bc(i_) % bs
                        Is_  = MOSE_domain % bc(i_) % is
                        Js_  = MOSE_domain % bc(i_) % js
                        Ks_  = MOSE_domain % bc(i_) % ks
                        Fs_  = MOSE_domain % bc(i_) % fs
                        if (Bs == Bm_ .and. Is == Im_ .and. Js == Jm_ .and. Ks == Km_ .and. Fs == Fm_) then
                          MOSE_domain % bc(i_) % ext_flux(np) = MOSE_domain % bc(i_) % ext_flux(np) - 0.5 * FUSS_new_flux
                        endif
                    end select
                  enddo
                enddo
              enddo
          end select
        enddo
      enddo
    enddo
    !$omp end parallel

  end subroutine Compute_FUSS_Flux


  subroutine Compute_MOSE_Flux ( FUSS_domain, MOSE_domain )
    use hydraMF_Lib_BC_Coupled_Wall_Temperature, only: MOSE_Wall_Temperature
    use FUSS_Advanced_Types_m
    use MOSE_Advanced_Types_m
    use MOSE_Global_m, only: nprim, np
    implicit none
    type(FUSS_domain_type), intent(inout) :: FUSS_domain
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    ! Local
    integer  :: b, f, lower, upper, i
    integer  :: Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs
    integer  :: b_, f_, lower_, upper_, i_
    integer  :: Bm_, Im_, Jm_, Km_, Fm_, Bs_, Is_, Js_, Ks_, Fs_
    real(R8) :: MOSE_new_flux(nprim)

    !$omp parallel
    upper = 0
    do b = 1, MOSE_domain%nb
      do f = 1, 6
        lower = upper + 1
        upper = upper + MOSE_domain%n_bf(b,f)
        !$omp do schedule ( dynamic ) private (i,Bm, Im, Jm, Km, Fm, Bs, Is, Js, Ks, Fs) &
        !$omp private (b_, f_, lower_, upper_, i_) &
        !$omp private (Bm_, Im_, Jm_, Km_, Fm_, Bs_, Is_, Js_, Ks_, Fs_) &
        !$omp private (MOSE_new_flux)
        do i = lower, upper
          select case ( MOSE_domain%bc(i)%type )
            case (1000)
              Bm  = MOSE_domain % bc(i) % b
              Im  = MOSE_domain % bc(i) % i
              Jm  = MOSE_domain % bc(i) % j 
              Km  = MOSE_domain % bc(i) % k 
              Fm  = MOSE_domain % bc(i) % f
              Bs  = MOSE_domain % bc(i) % bs
              Is  = MOSE_domain % bc(i) % is
              Js  = MOSE_domain % bc(i) % js
              Ks  = MOSE_domain % bc(i) % ks
              Fs  = MOSE_domain % bc(i) % fs
              call MOSE_Wall_Temperature( Im, Jm, Km, Fm, MOSE_domain % blk(Bm), &
                                          Is, Js, Ks, Fs, FUSS_domain % blk(Bs), &
                                          MOSE_domain%bc(i)%Mg(1), MOSE_domain%bc(i)%Pg, MOSE_new_flux )
              MOSE_domain % bc(i) % ext_flux = MOSE_domain % bc(i) % ext_flux + 0.5 * MOSE_new_flux
              
              upper_ = 0
              do b_ = 1, FUSS_domain%nb
                do f_ = 1, 6
                  lower_ = upper_ + 1
                  upper_ = upper_ + FUSS_domain%n_bf(b_,f_)
                  do i_ = lower_, upper_
                    select case ( FUSS_domain%bc(i_)%type )
                      case (1000)  ! Multi-Solver Connection
                        Bm_  = FUSS_domain % bc(i_) % b
                        Im_  = FUSS_domain % bc(i_) % i
                        Jm_  = FUSS_domain % bc(i_) % j 
                        Km_  = FUSS_domain % bc(i_) % k 
                        Fm_  = FUSS_domain % bc(i_) % f
                        Bs_  = FUSS_domain % bc(i_) % bs
                        Is_  = FUSS_domain % bc(i_) % is
                        Js_  = FUSS_domain % bc(i_) % js
                        Ks_  = FUSS_domain % bc(i_) % ks
                        Fs_  = FUSS_domain % bc(i_) % fs
                        if (Bs == Bm_ .and. Is == Im_ .and. Js == Jm_ .and. Ks == Km_ .and. Fs == Fm_) then
                          FUSS_domain % bc(i_) % ext_flux = FUSS_domain % bc(i_) % ext_flux - 0.5 * MOSE_new_flux(np)
                        endif
                    end select
                  enddo
                enddo
              enddo
          end select
        enddo
      enddo
    enddo
    !$omp end parallel

  end subroutine Compute_MOSE_Flux

end module hydraMF_Lib_BC_Coupled_Fluxes