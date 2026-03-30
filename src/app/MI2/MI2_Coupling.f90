module hydraMI2_Mod_Coupling
  implicit none
  private

  public:: ICExMOSE
  public:: MOSExICE
  public:: IGLOOxMOSE
  public:: MOSExIGLOO
  public:: initializeResBlock

  type, public :: resField
    real(8), dimension(:,:,:,:), allocatable :: resBlock
  end type resField

contains

  !! ICE coupling

  
  subroutine MOSExICE ( MOSE_domain, ICE_grid )
    use MOSE_Global_m, only: nsc, nu, nw, np
    use MOSE_Types_m
    use FLINT_Lib_Thermodynamic, only: EOS, co_rotot_Rtot, f_gamma, co_k_mi_lam_Wilke
    use ICE_Mod_Type, only: obj_grid
    implicit none
    type(MOSE_domain_type), intent(in)   :: MOSE_domain
    type(obj_grid), intent(inout)        :: ICE_grid
    integer(kind=I4) :: b, i, j, k 
    real(8) :: rhoi(nsc), rho, R, gam, T, mu, kappa, p

    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP SHARED(MOSE_domain,ICE_grid,nsc,nu,nw,np), &
    !$OMP PRIVATE(rhoi,rho,p,kappa,mu,R,gam,T)
    do b = 1, size(MOSE_domain%blk)
      !$OMP DO COLLAPSE(3)
      do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
          do i = 1, MOSE_domain%blk(b)%dim(1)

          rhoi = MOSE_domain%blk(b)%P(1:nsc,i,j,k)
          p = MOSE_domain%blk(b)%P(np,i,j,k)
          call co_rotot_Rtot(rhoi,rho,R)
          gam = f_gamma ( rhoi, p, rho, R )
          T = EOS(p,rho,R)
          call co_k_mi_lam_Wilke(rhoi,rho,T,mu,kappa)

          ! Density
          ICE_grid%block(b)%gas_phase%prim(1,i,j,k) = rho
          ! Velocity (3 components)
          ICE_grid%block(b)%gas_phase%prim(2:4,i,j,k) = MOSE_domain%blk(b)%P(nu:nw,i,j,k)
          ! Temperature
          ICE_grid%block(b)%gas_phase%prim(5,i,j,k) = T
          ! gamma
          ICE_grid%block(b)%gas_phase%gam(i,j,k) = gam
          ! R
          ICE_grid%block(b)%gas_phase%R(i,j,k) = R
          ! mu
          ICE_grid%block(b)%gas_phase%mu(i,j,k) = mu
          ! k
          ICE_grid%block(b)%gas_phase%k(i,j,k) = kappa

          enddo
        enddo
      enddo
      !$OMP END DO
    enddo
    !$OMP END PARALLEL

  end subroutine MOSExICE


  subroutine ICExMOSE ( MOSE_domain, ICE_grid )
    use MOSE_Types_m
    use MOSE_Global_m, only: nsc, nu, nw, np
    use FLINT_Lib_Thermodynamic, only: EOS, co_rotot_Rtot, f_gamma, co_k_mi_lam_Wilke
    use ICE_Lib_Model, only: assign_all
    use ICE_Mod_Parameters, only: ngroups, ncond
    use ICE_Mod_Sources
    use ICE_Mod_Type, only: obj_grid
    implicit none
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    type(obj_grid), intent(in)   :: ICE_grid
    integer :: b, i, j, k, g
    real(8) :: partial_residuals(nw+1)
    real(8) :: rhoi(nsc), rho, u(3), T, p
    real(8) :: R, gam, kappa, mu
    real(8) :: sources(6), dummy


    do b = 1, size(MOSE_domain%blk)

      do g = 1, ngroups

        call assign_all(g)

        !$OMP DO COLLAPSE(3) PRIVATE(partial_residuals,rhoi,rho,u,p,T,R,gam,kappa,mu,sources,i,j,k)
        do k = 1, MOSE_domain % blk(b) % dim(3)
        do j = 1, MOSE_domain % blk(b) % dim(2)
        do i = 1, MOSE_domain % blk(b) % dim(1)

          rhoi = MOSE_domain % blk(b) % P(1:nsc,i,j,k)
          u = MOSE_domain % blk(b) % P(nu:nw,i,j,k)
          p = MOSE_domain % blk(b) % P(np,i,j,k)
          call co_rotot_Rtot(rhoi, rho, R)
          gam = f_gamma ( rhoi, p, rho, R )
          T = EOS(p,rho,R)
          call co_k_mi_lam_Wilke(rhoi,rho,T,mu,kappa)

          call compute_source_ ( ICE_grid%block(b)%cond_phase(g)%prim(1:ncond(g),i,j,k), &
                                 dummy, &
                                 [rho,u,T], &
                                 R, gam, kappa, mu, &
                                 sources )

          partial_residuals(1:nsc) = 0._R8
          partial_residuals(nu:nw) = - sources(2:4)
          partial_residuals(np)    = - sources(ncond(g)-1)

          MOSE_domain % blk(b) % R(:,i,j,k) = MOSE_domain % blk(b) % R(:,i,j,k) - partial_residuals * MOSE_domain % blk(b) % vol(i,j,k)

        enddo ; enddo ; enddo
        !$OMP END DO
      enddo

    enddo

  end subroutine ICExMOSE


  !! IGLOO coupling


  subroutine MOSExIGLOO ( MOSE_domain, IGLOO )
    use MOSE_Global_m, only: nsc, nu, nw, np
    use MOSE_Types_m
    use FLINT_Lib_Thermodynamic, only: EOS, co_rotot_Rtot, f_gamma, co_k_mi_lam_Wilke
    use IGLOO_module, only: obj_IGLOO
    implicit none
    type(MOSE_domain_type), intent(in)   :: MOSE_domain
    type(obj_IGLOO), intent(inout)       :: IGLOO
    integer(kind=I4) :: b, i, j, k 
    real(8) :: rhoi(nsc), rho, R, gam, T, mu, kappa, p

    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP SHARED(MOSE_domain,IGLOO,nsc,nu,nw,np), &
    !$OMP PRIVATE(rhoi,rho,p,kappa,mu,R,gam,T)
    do b = 1, size(MOSE_domain%blk)
      !$OMP DO COLLAPSE(3)
      do k = 1, MOSE_domain%blk(b)%dim(3)
        do j = 1, MOSE_domain%blk(b)%dim(2)
          do i = 1, MOSE_domain%blk(b)%dim(1)

          rhoi = MOSE_domain%blk(b)%P(1:nsc,i,j,k)
          p = MOSE_domain%blk(b)%P(np,i,j,k)
          call co_rotot_Rtot(rhoi,rho,R)
          gam = f_gamma ( rhoi, p, rho, R )
          T = EOS(p,rho,R)
          call co_k_mi_lam_Wilke(rhoi,rho,T,mu,kappa)

          IGLOO%block(b)%density(1:nsc,i,j,k) = rhoi                               ! Density
          IGLOO%block(b)%velocity(1:3,i,j,k) = MOSE_domain%blk(b)%P(nu:nw,i,j,k)   ! Velocity (3 components)
          IGLOO%block(b)%temperature(i,j,k) = T                                    ! Temperature
          IGLOO%block(b)%gam(i,j,k) = gam                                          ! gamma
          IGLOO%block(b)%R(i,j,k) = R                                              ! R
          IGLOO%block(b)%mil(i,j,k) = mu                                           ! mu
          IGLOO%block(b)%kl(i,j,k) = kappa                                         ! k

          enddo
        enddo
      enddo
      !$OMP END DO
    enddo
    !$OMP END PARALLEL

  end subroutine MOSExIGLOO


  subroutine IGLOOxMOSE ( MOSE_domain, IGLOO, eulCoupling, relaxation, iter, iMax, relCoeff, oldResidual)
    use MOSE_Types_m
    use MOSE_Global_m, only: nsc, nu, nw, np
    use FLINT_Lib_Thermodynamic, only: EOS, co_rotot_Rtot, f_gamma, co_k_mi_lam_Wilke
    use IGLOO_module, only: obj_IGLOO
    use IGLOO_data_block, only: obj_sourceblock
    use IGLOO_variables, only: nfam
    use Lib_Equations, only: FdQd_euler
    implicit none
    type(MOSE_domain_type), intent(inout) :: MOSE_domain
    type(resField),         intent(inout) :: oldResidual(:)
    type(obj_IGLOO),        intent(in)    :: IGLOO
    real(8),                intent(in)    :: relCoeff
    integer,                intent(in)    :: iter, iMax
    logical,                intent(in)    :: relaxation, eulCoupling
    real(8), parameter :: toll=1.d-15
    logical :: useRelax
    integer :: b, i, j, k, m, fam
    real(8) :: partial_residuals(nw+1), sources(4)
    real(8) :: rhoi(nsc), rho, u(3), T, p
    real(8) :: R, gam, kappa, mu, coeff

    useRelax = (relCoeff>toll .and. relCoeff<(1.d0-toll)) .and. relaxation
    if (iter < iMax*(1.d0-relCoeff) .and. useRelax) then
      coeff = real(iter) / iMax
      coeff = coeff*(1.d0/(1.d0-relCoeff)-exp(-(coeff/relCoeff))**2)
      coeff = max(0.d0, min(1.d0, coeff))
    elseif ((iMax==1) .and. useRelax) then
    ! if (useRelax) then
      coeff = relCoeff
    else
      coeff = 1.d0
      useRelax = .false.
    endif

    do b = 1, size(MOSE_domain%blk)

      !$OMP DO COLLAPSE(3) PRIVATE(partial_residuals,sources,rhoi,rho,u,p,T,R,gam,kappa,mu,i,j,k,m)
      do k = 1, MOSE_domain % blk(b) % dim(3)
      do j = 1, MOSE_domain % blk(b) % dim(2)
      do i = 1, MOSE_domain % blk(b) % dim(1)

        if (eulCoupling) then

          rhoi = MOSE_domain % blk(b) % P(1:nsc,i,j,k)
          u    = MOSE_domain % blk(b) % P(nu:nw,i,j,k)
          p    = MOSE_domain % blk(b) % P(np,i,j,k)
          call co_rotot_Rtot(rhoi, rho, R)
          gam = f_gamma ( rhoi, p, rho, R )
          T = EOS(p,rho,R)
          call co_k_mi_lam_Wilke(rhoi,rho,T,mu,kappa)
          sources = 0._R8
          do m = 1,  nfam
            sources = sources + FdQd_euler([rho,u,T,gam,R,mu,kappa], IGLOO%euler(b,m)%velocity(:,i,j,k), &
                                IGLOO%euler(b,m)%temperature(i,j,k), IGLOO%euler(b,m)%density(i,j,k), IGLOO%euler(b,m)%np(i,j,k), 1)
                                !!!!>>>> TO DOOOO :: DEFINISCI MATERIALIIII COME CRISTO COMANDA
          enddo
          partial_residuals(:) = 0._R8
          ! partial_residuals(1:nsc) = 0._R8
          partial_residuals(nu:nw) = - sources(1:3)
          partial_residuals(np)    = - sources(4)
          partial_residuals = partial_residuals * MOSE_domain % blk(b) % vol(i,j,k)
          
        else

          partial_residuals(:) = 0._R8
          ! partial_residuals(1:nsc) = IGLOO%source(b)%sourceMass(:,i,j,k)
          partial_residuals(nu:nw) = IGLOO%source(b)%sourceMom(:,i,j,k)
          partial_residuals(np)    = IGLOO%source(b)%sourceEn(i,j,k)

        endif

        if (useRelax) then
          partial_residuals(1:nsc) = coeff * partial_residuals(1:nsc) + (1._R8 - coeff) * oldResidual(b)%resBlock(1:nsc,i,j,k)
          partial_residuals(nu:nw) = coeff * partial_residuals(nu:nw) + (1._R8 - coeff) * oldResidual(b)%resBlock(nu:nw,i,j,k) 
          partial_residuals(np)    = coeff * partial_residuals(np)    + (1._R8 - coeff) * oldResidual(b)%resBlock( np  ,i,j,k)
          oldResidual(b)%resBlock(:,i,j,k) = partial_residuals(:)
        endif

        MOSE_domain % blk(b) % R(:,i,j,k) = MOSE_domain % blk(b) % R(:,i,j,k) - partial_residuals

      enddo ; enddo ; enddo
      !$OMP END DO

    enddo

  end subroutine IGLOOxMOSE


  subroutine initializeResBlock(MOSE_simulation, oldResidual)
    use MOSE_Types_m,  only: MOSE_simulation_type
    use MOSE_Global_m, only: nsc, nu, nw, np
    implicit none
    type(MOSE_simulation_type),  intent(in)  :: MOSE_simulation
    type(resField), allocatable, intent(out) :: oldResidual(:)
    integer :: Ni, Nj, Nk, b

    allocate(oldResidual(1:size(MOSE_simulation%IOfield(1)%block)))
    do b = 1, size(oldResidual)
      Ni = MOSE_simulation%IOfield(1)%block(b)%Ni
      Nj = MOSE_simulation%IOfield(1)%block(b)%Nj
      Nk = MOSE_simulation%IOfield(1)%block(b)%Nk
      if (.not.allocated(oldResidual(b)%resBlock)) then
        allocate(oldResidual(b)%resBlock(1:nw+1, 1:Ni, 1:Nj, 1:Nk))
        oldResidual(b)%resBlock(:,:,:,:) = 0.d0
      endif
    enddo

  end subroutine initializeResBlock

end module hydraMI2_Mod_Coupling
