MODULE mod_thermodynamics
    IMPLICIT NONE
    ! ratio of specific heats
    REAL(kind=8), PARAMETER :: gamma = 1.4d00



    ! index of variables in the solution matrix
    INTEGER(kind=4), PARAMETER :: i_rho   =  1
    INTEGER(kind=4), PARAMETER :: i_u     =  2
    INTEGER(kind=4), PARAMETER :: i_v     =  3
    INTEGER(kind=4), PARAMETER :: i_P     =  4
    INTEGER(kind=4), PARAMETER :: i_rho_u =  5  ! x momentum
    INTEGER(kind=4), PARAMETER :: i_rho_v =  6  ! y momentum
    INTEGER(kind=4), PARAMETER :: i_eT    =  7  ! total energy XXX ???
    INTEGER(kind=4), PARAMETER :: i_T     =  8  ! ?
    INTEGER(kind=4), PARAMETER :: i_schl  =  9
    !    INTEGER(kind=4), PARAMETER :: i_e    = 1

    ! number of variables
    INTEGER(kind=4), PARAMETER :: nvar      = 9

CONTAINS


    ELEMENTAL FUNCTION specific_energy(P, rho) RESULT(e)

        IMPLICIT NONE

        REAL(KIND=8), INTENT(IN) :: P, rho
        REAL(KIND=8) :: e

        e = P / ((gamma - 1) * rho)

    END FUNCTION specific_energy

    ELEMENTAL FUNCTION pressure(e, rho) RESULT(P)

        IMPLICIT NONE

        REAL(KIND=8), INTENT(IN)  :: e, rho
        REAL(KIND=8) :: P

        P = (gamma - 1) * rho * e

    END FUNCTION pressure


END MODULE mod_thermodynamics
