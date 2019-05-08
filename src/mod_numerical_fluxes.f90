MODULE mod_numerical_fluxes

    USE mod_euler_flux_jacobian
    USE mod_constants
    USE mod_exact_riemann

    IMPLICIT NONE

CONTAINS

    SUBROUTINE godunov_flux_x(w, F)
        IMPLICIT NONE

        REAL(dp_kind), DIMENSION(:,:,:), INTENT(IN) :: w
        REAL(dp_kind), DIMENSION(:,:,:), INTENT(INOUT) :: F

        REAL(dp_kind), DIMENSION(4) :: wl, wr, wCl, wCr, w_vert, lambda_l, lambda_r, lambda

        REAL(dp_kind) :: vl, ul,  vr, ur,  vCl, vCr, uC,  ss ! shock speed

        INTEGER :: i, j

        LOGICAL, PARAMETER :: dir = .TRUE.

        REAL(dp_kind), PARAMETER :: rel_tol = 1.0d-8

        DO j = 1, SIZE(w, 2) ! ciclo esterno su y
            DO i = 1, SIZE(w,1)-1 ! ciclo interno su x

                ! recuperiamo lo stato sinistro wl e lo stato destro wl
                wl = w(i,j,i_cons)
                wr = w(i+1,j,i_cons)

                ! controlliamo che ci sia salto
                IF (SQRT(SUM((wr - wl)**2))  <  SQRT(SUM((wl + wr)**2)) * rel_tol) THEN

                    w_vert = wl     ! idem wr

                ELSE

                    vl = 1/wl(1);   ul = wl(2)/wl(1)
                    vr = 1/wr(1);   ur = wr(2)/wr(1)

                    CALL exact_riemann (wl, wr,  wCl, wCr, dir )

                    WRITE(*,*) 'F: i,j = ', i,j

                    uC = wCl(2)/wCl(1) ! ;  PC = Pi(wCl)

                    IF (uC > 0) THEN  ! The contact discontinuity propagates to
                                       ! the right: the third wave is not relevant

                        vCl = 1/wCl(1)

                        IF (vCl < vl) THEN ! LEFT SHOCK

                            ! Compute the shock speed
                            ss = (vl*uC - vCl*ul)/(vl - vCl)

                            IF (ss > 0) THEN ! right propagating shock
                                w_vert = wl
                            ELSE             ! left propagating shock
                                w_vert = wCl
                            ENDIF

                        ELSE ! LEFT RAREFACTION

                            lambda_l = eigenvalues(wl, dir )
                            lambda   = eigenvalues(wCl, dir )

                            IF (lambda_l(1) >= 0) THEN ! right propagating fan
                                w_vert = wl
                            ELSEIF (lambda(1) <= 0) THEN ! left propagating fan
                                w_vert = wCl
                            ELSE ! transonic rarefaction for the first wave
                                w_vert = ws_sonic_state(1, wl, dir )
                            ENDIF


                        ENDIF ! alternative shock  vs  fan for the first wave


                    ELSE ! (uC < 0) The contact discontinuity propagates to
                          !          the left: the first wave is not relevant

                        vCr = 1/wCr(1)

                        IF (vCr < vr) THEN ! RIGHT SHOCK

                            ! Compute the shock speed
                            ss = (vr*uC - vCr*ur)/(vr - vCr)

                            IF (ss > 0) THEN ! right propagating shock
                                w_vert = wCr
                            ELSE             ! left propagating shock
                                w_vert = wr
                            ENDIF


                        ELSE ! RIGHT RAREFACTION

                            lambda   = eigenvalues(wCr, dir )
                            lambda_r = eigenvalues(wr, dir )

                            IF (lambda(4) >= 0) THEN ! right propagating fan
                                w_vert = wCr
                            ELSEIF (lambda_r(4) <= 0) THEN ! left propagating fan
                                w_vert = wr
                            ELSE ! transonic rarefaction for the third wave
                                w_vert = ws_sonic_state(4, wr, dir )
                            ENDIF


                        ENDIF ! alternative shock  vs  fan for the fourth wave


                    ENDIF ! alternative uC >  or  uC <= 0


                ENDIF ! il salto è rilevante

                F(i,j,:) = flux_x(w_vert)

            END DO

        END DO

    END SUBROUTINE godunov_flux_x


    SUBROUTINE godunov_flux_y(w, G)
        IMPLICIT NONE

        REAL(dp_kind), DIMENSION(:,:,:), INTENT(IN) :: w
        REAL(dp_kind), DIMENSION(:,:,:), INTENT(INOUT) :: G

        REAL(dp_kind), DIMENSION(4) :: wl, wr, wCl, wCr, w_vert, lambda_l, lambda_r, lambda

        REAL(dp_kind) :: vl, ul,  vr, ur,  vCl, vCr, uC,  ss ! shock speed

        INTEGER :: i, j

        LOGICAL, PARAMETER :: dir = .FALSE.

        REAL(dp_kind), PARAMETER :: rel_tol = 1.0d-8

        DO i = 1, SIZE(w, 2) ! ciclo esterno su y
            DO j = 1, SIZE(w,1)-1 ! ciclo interno su x

                ! recuperiamo lo stato sinistro wl e lo stato destro wl
                wl = w(i,j,i_cons)
                wr = w(i,j+1,i_cons)

                ! controlliamo che ci sia salto
                IF (SQRT(SUM((wr - wl)**2))  <  SQRT(SUM((wl + wr)**2)) * rel_tol) THEN

                    w_vert = wl     ! idem wr

                ELSE

                    vl = 1/wl(1);   ul = wl(3)/wl(1)
                    vr = 1/wr(1);   ur = wr(3)/wr(1)

                    CALL exact_riemann (wl, wr,  wCl, wCr, dir )

                    WRITE(*,*) 'G: i,j = ', i,j

                    uC = wCl(3)/wCl(1) ! ;  PC = Pi(wCl)

                    IF (uC > 0) THEN  ! The contact discontinuity propagates to
                                       ! the right: the third wave is not relevant

                        vCl = 1/wCl(1)

                        IF (vCl < vl) THEN ! LEFT SHOCK

                            ! Compute the shock speed
                            ss = (vl*uC - vCl*ul)/(vl - vCl)

                            IF (ss > 0) THEN ! right propagating shock
                                w_vert = wl
                            ELSE             ! left propagating shock
                                w_vert = wCl
                            ENDIF

                        ELSE ! LEFT RAREFACTION

                            lambda_l = eigenvalues(wl, dir )
                            lambda   = eigenvalues(wCl, dir )

                            IF (lambda_l(1) >= 0) THEN ! right propagating fan
                                w_vert = wl
                            ELSEIF (lambda(1) <= 0) THEN ! left propagating fan
                                w_vert = wCl
                            ELSE ! transonic rarefaction for the first wave
                                w_vert = ws_sonic_state(1, wl, dir )
                            ENDIF


                        ENDIF ! alternative shock  vs  fan for the first wave


                    ELSE ! (uC < 0) The contact discontinuity propagates to
                          !          the left: the first wave is not relevant

                        vCr = 1/wCr(1)

                        IF (vCr < vr) THEN ! RIGHT SHOCK

                            ! Compute the shock speed
                            ss = (vr*uC - vCr*ur)/(vr - vCr)

                            IF (ss > 0) THEN ! right propagating shock
                                w_vert = wCr
                            ELSE             ! left propagating shock
                                w_vert = wr
                            ENDIF


                        ELSE ! RIGHT RAREFACTION

                            lambda   = eigenvalues(wCr, dir )
                            lambda_r = eigenvalues(wr, dir )

                            IF (lambda(4) >= 0) THEN ! right propagating fan
                                w_vert = wCr
                            ELSEIF (lambda_r(4) <= 0) THEN ! left propagating fan
                                w_vert = wr
                            ELSE ! transonic rarefaction for the third wave
                                w_vert = ws_sonic_state(4, wr, dir )
                            ENDIF


                        ENDIF ! alternative shock  vs  fan for the fourth wave


                    ENDIF ! alternative uC >  or  uC <= 0


                ENDIF ! il salto è rilevante

                G(i,j,:) = flux_y(w_vert)

            END DO

        END DO


    END SUBROUTINE godunov_flux_y


END MODULE mod_numerical_fluxes
