MODULE mod_thermodynamics

    USE mod_constants

    IMPLICIT NONE

CONTAINS

    SUBROUTINE conservatives(w)
        ! writes conservative variables from primitive ones
        IMPLICIT NONE

        REAL(dp_kind), DIMENSION(:,:,:), INTENT(INOUT) :: w

        ! write x, y momentum
        ! rho*u
        w(:,:,i_mx) = w(:,:,i_rho)*w(:,:,i_u)
        ! rho*v
        w(:,:,i_my) = w(:,:,i_rho)*w(:,:,i_v)

        ! specific internal energy
        w(:,:,i_e) = specific_energy(w(:,:,i_P),w(:,:,i_rho))

        ! eT = rho*(e + 0.5*(u**2 + v**2) )
        w(:,:,i_eT) = w(:,:,i_rho)*( w(:,:,i_e) + 0.5d0 * ( w(:,:,i_u)**2 + w(:,:,i_v)**2 ) )

    END SUBROUTINE conservatives


    SUBROUTINE primitive(w)
        ! writes primitive variables from conservative ones
        IMPLICIT NONE

        REAL(dp_kind), DIMENSION(:,:,:), INTENT(INOUT) :: w

        ! write velocities
        ! u
        w(:,:,i_u) = w(:,:,i_mx)/w(:,:,i_rho)
        ! v
        w(:,:,i_v) = w(:,:,i_my)*w(:,:,i_rho)

        ! specific internal energy: e = eT/rho - 0.5(u**2 +v**2)
        w(:,:,i_e) = w(:,:,i_eT)/w(:,:,i_rho) - 0.5d0 * ( w(:,:,i_u)**2 + w(:,:,i_v)**2 )
        ! Pressure
        w(:,:,i_P) = pressure(w(:,:,i_e),w(:,:,i_rho))

    END SUBROUTINE primitive


    SUBROUTINE pseudoschlieren(w)
        ! computes pseudo Schlieren of density rho
        IMPLICIT NONE

        REAL(dp_kind), DIMENSION(:,:,:), INTENT(INOUT) :: w

        REAL(dp_kind), ALLOCATABLE, DIMENSION(:,:) :: drho_dx, drho_dy, tmp_schl

        INTEGER :: istat

        ALLOCATE (drho_dx(num_cells,num_cells),STAT=istat)
        ALLOCATE (drho_dy(num_cells,num_cells),STAT=istat)
        ALLOCATE (tmp_schl(num_cells,num_cells),STAT=istat)

        drho_dx(1,:) = 0.0d0
        drho_dy(:,1) = 0.0d0

        drho_dx(2:num_cells,:) = ((w(2:num_cells,:,i_rho)-w(1:num_cells-1,:,i_rho))/Dx)**2
        drho_dy(:,2:num_cells) = ((w(:,2:num_cells,i_rho)-w(:,1:num_cells-1,i_rho))/Dy)**2

        tmp_schl = drho_dx + drho_dy

        w(:,:,i_schl) = tmp_schl / MAXVAL(tmp_schl)

        DEALLOCATE(drho_dx)
        DEALLOCATE(drho_dy)
        DEALLOCATE(tmp_schl)

    END SUBROUTINE pseudoschlieren




    ELEMENTAL FUNCTION specific_energy(P, rho) RESULT(e)

        IMPLICIT NONE

        REAL(dp_kind), INTENT(IN) :: P, rho
        REAL(dp_kind) :: e

        e = P / ((gamma - 1) * rho)

    END FUNCTION specific_energy

    ELEMENTAL FUNCTION pressure(e, rho) RESULT(P)

        IMPLICIT NONE

        REAL(dp_kind), INTENT(IN)  :: e, rho
        REAL(dp_kind) :: P

        P = (gamma - 1) * rho * e

    END FUNCTION pressure


END MODULE mod_thermodynamics
