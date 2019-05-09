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

                ! recuperiamo lo stato sinistro wl e lo stato destro wr
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


    SUBROUTINE HLLC_flux_x(w,F)
        ! HLLC numerical flux according to Toro 10.6 pagg. 331-332

        IMPLICIT NONE

        REAL(dp_kind), DIMENSION(:,:,:), INTENT(IN) :: w
        REAL(dp_kind), DIMENSION(:,:,:), INTENT(INOUT) :: F

        REAL(dp_kind), DIMENSION(SIZE(w,1)-1,SIZE(w,2)) :: p_star, p_pvrs, rho_bar, c_bar, q_l, q_r, S_l, S_r, S_star

        REAL(dp_kind), DIMENSION(4) :: w_star_k

        INTEGER :: r_min, r_max, l_min, l_max

        INTEGER :: i, j

        REAL(dp_kind), PARAMETER :: rel_tol = 1.0d-8

        r_min = 2
        r_max = num_cells

        l_min = 1
        l_max = num_cells - 1

        rho_bar = HALF*( w(l_min:l_max,:,i_rho) + w(r_min:r_max,:,i_rho) )

        c_bar =  HALF*( w(l_min:l_max,:,i_c) + w(r_min:r_max,:,i_c) )

        p_pvrs = HALF*( w(l_min:l_max,:,i_P) + w(r_min:r_max,:,i_P) ) - &
                 HALF*( w(r_min:r_max,:,i_u) - w(l_min:l_max,:,i_u) ) * rho_bar * c_bar

        WHERE (p_pvrs > 0)
            p_star = p_pvrs
        ELSEWHERE
            p_star = 0.0d0
        END WHERE

        q_l = ONE;    q_r = ONE

        WHERE (p_star > w(l_min:l_max,:,i_P))
            q_l = SQRT(ONE + (gamma + ONE)/(2.0*gamma)*( p_star/w(l_min:l_max,:,i_P) - ONE)  )
        END WHERE

        WHERE (p_star > w(r_min:r_max,:,i_P))
            q_r = SQRT(ONE + (gamma + ONE)/(2.0*gamma)*( p_star/w(r_min:r_max,:,i_P) - ONE)  )
        END WHERE


        S_l = w(l_min:l_max,:,i_u) - w(l_min:l_max,:,i_c)*q_l
        S_r = w(r_min:r_max,:,i_u) + w(r_min:r_max,:,i_c)*q_r

        ! S* = (pR - pL + rhoL uL (SL - uL )  - rhoR uR (SR - uR) ) / ( rhoL(SL -uL ) - rhoR (SR - uR ) )

        S_star = ( w(r_min:r_max,:,i_P) - w(l_min:l_max,:,i_P) + &
                   w(l_min:l_max,:,i_rho) * w(l_min:l_max,:,i_u) * (S_l - w(l_min:l_max,:,i_u) ) - &
                   w(r_min:r_max,:,i_rho) * w(r_min:r_max,:,i_u) * (S_r - w(r_min:r_max,:,i_u) ) &
                 ) / ( w(l_min:l_max,:,i_rho) * (S_l - w(l_min:l_max,:,i_u) ) - &
                       w(r_min:r_max,:,i_rho) * (S_r - w(r_min:r_max,:,i_u) ) )

        ! flux for left and right boundary
        DO j = 1, SIZE(w,2)
            F(1,j,:) = flux_x(w(1,j,i_cons))
            F(num_cells+1,j,:) = flux_x(w(num_cells,j,i_cons))
        END DO

        DO j = 1, SIZE(w,2)
            DO i = 1, SIZE(w,1) - 1
                IF (S_l(i,j) >= 0.0 ) THEN
                    F(i+1,j,:) = flux_x(w(i,j,i_cons)) ! FL
                ELSEIF ( S_star(i,j) >= 0.0 ) THEN
                    ! inizializzo le var note di w_star_k (U_star_k nel testo)
                    w_star_k(1) = ONE
                    w_star_k(2) = S_star(i,j)
                    w_star_k(3) = w(i,j,i_v)
                    w_star_k(4) = w(i,j,i_eT)/w(i,j,i_rho) + &                ! Ek/rho_k +
                                  (S_star(i,j) - w(i,j,i_u) ) *  &            ! (S_star - u_k) *
                                  (S_star(i,j) + w(i,j,i_P)/( w(i,j,i_rho)* ( S_l(i,j) - w(i,j,i_u)  ) )  )! ( S_star + p_K/(rho_k*(S_k - u_k)))

                    w_star_k = w_star_k*w(i,j,i_rho)*( (S_l(i,j) - w(i,j,i_u) ) / ( S_l(i,j) - S_star(i,j)  )   )

                    F(i+1,j,:) = flux_x(w(i,j,i_cons)) + S_l(i,j)*(w_star_k - w(i,j,i_cons))  ! F_star_L
                ELSEIF ( S_r(i,j ) <= 0.0 ) THEN
                    F(i+1,j,:) = flux_x(w(i+1,j,i_cons)) ! FR
                ELSE
                    ! inizializzo le var note di w_star_k (U_star_k nel testo)
                    w_star_k(1) = ONE
                    w_star_k(2) = S_star(i,j)
                    w_star_k(3) = w(i+1,j,i_v)
                    w_star_k(4) = w(i+1,j,i_eT)/w(i+1,j,i_rho) +   &            ! Ek/rho_k +
                                  (S_star(i,j) - w(i+1,j,i_u) ) *    &           ! (S_star - u_k ) *
                                  (S_star(i,j) + w(i+1,j,i_P)/( w(i+1,j,i_rho)*( S_r(i,j) - w(i+1,j,i_u) ) ) ) ! ( S_star + p_K/(rho_k*(S_k - u_k)))

                    w_star_k = w_star_k*w(i+1,j,i_rho)*( (S_r(i,j) - w(i+1,j,i_u) ) / ( S_r(i,j) - S_star(i,j)  )   )

                    F(i+1,j,:) = flux_x(w(i+1,j,i_cons)) + S_r(i,j)*(w_star_k - w(i+1,j,i_cons))  ! F_star_L
                END IF
            END DO
        END DO



    END SUBROUTINE HLLC_flux_x



    SUBROUTINE HLLC_flux_y(w,G)
        ! HLLC numerical flux according to Toro 10.6 pagg. 331-332

        IMPLICIT NONE

        REAL(dp_kind), DIMENSION(:,:,:), INTENT(IN) :: w
        REAL(dp_kind), DIMENSION(:,:,:), INTENT(INOUT) :: G

        REAL(dp_kind), DIMENSION(SIZE(w,1),SIZE(w,2)-1) :: p_star, p_pvrs, rho_bar, c_bar, q_l, q_r, S_l, S_r, S_star

        REAL(dp_kind), DIMENSION(4) :: w_star_k

        INTEGER :: r_min, r_max, l_min, l_max

        INTEGER :: i, j

        REAL(dp_kind), PARAMETER :: rel_tol = 1.0d-8

        r_min = 2
        r_max = num_cells

        l_min = 1
        l_max = num_cells - 1

        rho_bar = HALF*( w(:,l_min:l_max,i_rho) + w(:,r_min:r_max,i_rho) )

        c_bar =  HALF*( w(:,l_min:l_max,i_c) + w(:,r_min:r_max,i_c) )

        p_pvrs = HALF*( w(:,l_min:l_max,i_P) + w(:,r_min:r_max,i_P) ) - &
                 HALF*( w(:,r_min:r_max,i_v) - w(:,l_min:l_max,i_v) ) * rho_bar * c_bar

        WHERE (p_pvrs > 0)
            p_star = p_pvrs
        ELSEWHERE
            p_star = 0.0d0
        END WHERE

        q_l = ONE;    q_r = ONE

        WHERE (p_star > w(:,l_min:l_max,i_P))
            q_l = SQRT(ONE + (gamma + ONE)/(2.0*gamma)*( p_star/w(:,l_min:l_max,i_P) - ONE)  )
        END WHERE

        WHERE (p_star > w(:,r_min:r_max,i_P))
            q_r = SQRT(ONE + (gamma + ONE)/(2.0*gamma)*( p_star/w(:,r_min:r_max,i_P) - ONE)  )
        END WHERE


        S_l = w(:,l_min:l_max,i_v) - w(:,l_min:l_max,i_c)*q_l
        S_r = w(:,r_min:r_max,i_v) + w(:,r_min:r_max,i_c)*q_r

        ! S* = (pR - pL + rhoL uL (SL - uL )  - rhoR uR (SR - uR) ) / ( rhoL(SL -uL ) - rhoR (SR - uR ) )

        S_star = ( w(:,r_min:r_max,i_P) - w(:,l_min:l_max,i_P) + &
                   w(:,l_min:l_max,i_rho) * w(:,l_min:l_max,i_v) * (S_l - w(:,l_min:l_max,i_v) ) - &
                   w(:,r_min:r_max,i_rho) * w(:,r_min:r_max,i_v) * (S_r - w(:,r_min:r_max,i_v) ) &
                 ) / ( w(:,l_min:l_max,i_rho) * (S_l - w(:,l_min:l_max,i_v) ) - &
                       w(:,r_min:r_max,i_rho) * (S_r - w(:,r_min:r_max,i_v) ) )

        ! flux for lower and upper boundary
        DO i = 1, SIZE(w,1)
            G(i,1,:) = flux_y(w(i,1,i_cons))
            G(i,num_cells+1,:) = flux_y(w(i,num_cells,i_cons))
        END DO


        DO j = 1, SIZE(w,2) - 1
            DO i = 1, SIZE(w,1)
                IF (S_l(i,j) >= 0.0 ) THEN
                    G(i,j+1,:) = flux_y(w(i,j,i_cons)) ! FL
                ELSEIF ( S_star(i,j) >= 0.0 ) THEN
                    ! inizializzo le var note di w_star_k (U_star_k nel testo)
                    w_star_k(1) = ONE
                    w_star_k(2) = w(i,j,i_u)
                    w_star_k(3) = S_star(i,j)
                    w_star_k(4) = w(i,j,i_eT)/w(i,j,i_rho) + &                ! Ek/rho_k +
                                  (S_star(i,j) - w(i,j,i_v) ) *  &            ! (S_star - u_k) *
                                  (S_star(i,j) + w(i,j,i_P)/( w(i,j,i_rho)* ( S_l(i,j) - w(i,j,i_v)  ) )  )! ( S_star + p_K/(rho_k*(S_k - u_k)))

                    w_star_k = w_star_k*w(i,j,i_rho)*( (S_l(i,j) - w(i,j,i_v) ) / ( S_l(i,j) - S_star(i,j)  )   )

                    G(i,j+1,:) = flux_y(w(i,j,i_cons)) + S_l(i,j)*(w_star_k - w(i,j,i_cons))  ! F_star_L
                ELSEIF ( S_r(i,j ) <= 0.0 ) THEN
                    G(i,j+1,:) = flux_y(w(i,j+1,i_cons)) ! FR
                ELSE
                    ! inizializzo le var note di w_star_k (U_star_k nel testo)
                    w_star_k(1) = ONE
                    w_star_k(2) = w(i,j+1,i_u)
                    w_star_k(3) = S_star(i,j)
                    w_star_k(4) = w(i,j+1,i_eT)/w(i,j+1,i_rho) +   &            ! Ek/rho_k +
                                  (S_star(i,j) - w(i,j+1,i_v) ) *    &           ! (S_star - u_k ) *
                                  (S_star(i,j) + w(i,j+1,i_P)/( w(i,j+1,i_rho)*( S_r(i,j) - w(i,j+1,i_v) ) ) ) ! ( S_star + p_K/(rho_k*(S_k - u_k)))

                    w_star_k = w_star_k*w(i,j+1,i_rho)*( (S_r(i,j) - w(i,j+1,i_v) ) / ( S_r(i,j) - S_star(i,j)  )   )

                    G(i,j+1,:) = flux_y(w(i,j+1,i_cons)) + S_r(i,j)*(w_star_k - w(i,j+1,i_cons))  ! F_star_L
                END IF
            END DO
        END DO



    END SUBROUTINE HLLC_flux_y




END MODULE mod_numerical_fluxes
