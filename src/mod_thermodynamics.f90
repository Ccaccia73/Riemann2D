MODULE mod_thermodynamics

    USE mod_constants

    IMPLICIT NONE

CONTAINS


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
