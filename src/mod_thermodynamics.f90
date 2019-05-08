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

    ELEMENTAL FUNCTION sound_speed(e, rho) RESULT(c)

        IMPLICIT NONE

        REAL(dp_kind), INTENT(IN) :: e
        REAL(dp_kind), INTENT(IN), OPTIONAL :: rho
        REAL(dp_kind) :: c

        c = SQRT(gamma * (gamma - 1) * e)

    END FUNCTION sound_speed


    FUNCTION pressure_RAR(P_, v_, v,  dP_dv) RESULT(P)

        IMPLICIT NONE

        REAL(dp_kind), INTENT(IN) :: P_, v_, v
        REAL(dp_kind), INTENT(OUT), OPTIONAL :: dP_dv
        REAL(dp_kind) :: P

        P = P_ * (v_ / v)**gamma

        IF (PRESENT(dP_dv)) dP_dv = - gamma * P/v

    END FUNCTION pressure_RAR

    FUNCTION integral_c_v_isentropic(P_, v_, v) RESULT(integ)

        ! calcolo di int_{v_}^v c(s_, v)/v dv

        IMPLICIT NONE

        REAL(dp_kind), INTENT(IN) :: P_, v_, v
        REAL(dp_kind) :: integ

        integ =  (2*SQRT(gamma*P_*v_)/(gamma-1)) * (1 - (v_/v)**((gamma-1)/2))

    END FUNCTION integral_c_v_isentropic



    FUNCTION sound_speed_isentropic(P_, v_, v) RESULT(c_isent)

        IMPLICIT NONE

        REAL(dp_kind), INTENT(IN) :: P_, v_, v
        REAL(dp_kind) :: c_isent

        c_isent = SQRT(gamma*P_*v_) * (v_/v)**((gamma-1)/2)

    END FUNCTION sound_speed_isentropic



    FUNCTION pressure_RH(P_, v_, v,  dP_dv) RESULT(P)

        IMPLICIT NONE

        REAL(dp_kind), INTENT(IN) :: P_, v_, v
        REAL(dp_kind), INTENT(OUT), OPTIONAL :: dP_dv
        REAL(dp_kind) :: P

        REAL(dp_kind) :: NUM, DEN

        ! Controllo che v non sia inferiore all'asintoto (Vedi testo)

        NUM = (gamma + 1) * v_  -  (gamma - 1) * v
        DEN = (gamma + 1) * v   -  (gamma - 1) * v_

        P = P_ * NUM/DEN

        IF (PRESENT(dP_dv)) dP_dv = -4 * gamma * P_ * v_/DEN**2

    END FUNCTION pressure_RH

END MODULE mod_thermodynamics
