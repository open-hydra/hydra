module hydra_Global_m

    integer, parameter :: gc_ = 2
    integer, dimension(6,3), parameter ::  guide_  = reshape((/ 1, 0, 0, &
                                                            -1, 0, 0, &
                                                             0, 1, 0, &
                                                             0,-1, 0, &
                                                             0, 0, 1, &
                                                             0, 0,-1 /), shape(guide_), order=(/2,1/) )

end module hydra_Global_m
