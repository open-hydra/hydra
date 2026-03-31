module hydra_Parameters_m

  implicit none

  integer, parameter  :: hlen_=1024
  integer, parameter  :: llen_=256
  integer, parameter  :: clen_=16
  integer, dimension(6,3), parameter ::  guide_  = reshape( [ 1, 0, 0, &
                                                             -1, 0, 0, &
                                                              0, 1, 0, &
                                                              0,-1, 0, &
                                                              0, 0, 1, &
                                                              0, 0,-1  ] , shape(guide_), order=[2,1] )

end module hydra_Parameters_m