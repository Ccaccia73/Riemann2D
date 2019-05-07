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

        REAL(dp_kind), DIMENSION(4) :: wl, wr, w_vert, lambda_l, lambda_r, wCl, wCr

        INTEGER :: i, j

        REAL(dp_kind), PARAMETER :: rel_tol = 1.0d-8

        DO j = 1, SIZE(w, 2) ! ciclo esterno su y
            DO i = 1, SIZE(w,1)-1 ! ciclo interno su x

                ! recuperiamo lo stato sinistro wl e lo stato destro wl
                wl = w(i,j,i_cons)
                wr = w(i+1,j,i_cons)

                ! controlliamo che ci sia salto
                IF (SQRT(SUM((wr - wl)**2))  <  SQRT(SUM((wl + wr)**2)) * rel_tol) THEN

                    w_vert = (wl + wr)/2

                ELSE
                    ! calcoliamo gli autovalori, funzione eigenvalues data dall'esterno
                    lambda_l = eigenvalues(wl, .TRUE.)
                    lambda_r = eigenvalues(wr, .TRUE.)

                    ! escludiamo i casi semplici e, se serve, risolviamo riemann
                    IF (lambda_l(1) > 0  .AND.  lambda_r(1) > 0) THEN

                        w_vert = wl ! tira tutto a destra

                    ELSEIF (lambda_l(4) < 0  .AND.  lambda_r(4) < 0) THEN

                        w_vert = wr ! tira tutto a sinistra

                    ELSE ! necessario risolvere il problema di Riemann

                        CALL exact_riemann (wl, wr,  wCl, wCr, .TRUE. )

                        w_vert = vertical_state(wl, wCl, wCr, wr)

                    ENDIF

                ENDIF

                F(i,j,:) = flux_x(w_vert)

            END DO

        END DO

    END SUBROUTINE godunov_flux_x































    SUBROUTINE godunov_flux_y(w, G)
        IMPLICIT NONE

        REAL(dp_kind), DIMENSION(:,:,:), INTENT(IN) :: w
        REAL(dp_kind), DIMENSION(:,:,:), INTENT(INOUT) :: G




    END SUBROUTINE godunov_flux_y


END MODULE mod_numerical_fluxes
