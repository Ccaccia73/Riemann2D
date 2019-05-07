MODULE mod_constants
    IMPLICIT NONE

    INTEGER, PARAMETER :: sp_kind   = kind(1.0)  !< single precision
    INTEGER, PARAMETER :: dp_kind   = kind(1.d0) !< double precision
    INTEGER, PARAMETER :: sh_kind   = selected_int_kind(2)
    INTEGER, PARAMETER :: int_kind  = selected_int_kind(4)
    INTEGER, PARAMETER :: long_kind  = selected_int_kind(8)


    ! ratio of specific heats
    REAL(dp_kind), PARAMETER :: gamma = 1.4d00



    ! index of variables in the solution matrix
    INTEGER(sh_kind), PARAMETER :: i_rho   =  1
    INTEGER(sh_kind), PARAMETER :: i_u     =  2
    INTEGER(sh_kind), PARAMETER :: i_v     =  3
    INTEGER(sh_kind), PARAMETER :: i_P     =  4
    INTEGER(sh_kind), PARAMETER :: i_mx    =  5  ! x momentum
    INTEGER(sh_kind), PARAMETER :: i_my    =  6  ! y momentum
    INTEGER(sh_kind), PARAMETER :: i_eT    =  7  ! total energy
    INTEGER(sh_kind), PARAMETER :: i_e     =  8  ! internal energy
    INTEGER(sh_kind), PARAMETER :: i_c     =  9  ! sound speed
    INTEGER(sh_kind), PARAMETER :: i_schl  =  10
    !    INTEGER(kind=4), PARAMETER :: i_e    = 1

    ! indexes of conservative variables
    INTEGER(sh_kind), PARAMETER, DIMENSION(4) :: i_cons = (/i_rho, i_mx, i_my, i_eT/)

    ! number of variables
    INTEGER(sh_kind), PARAMETER :: nvar      = 10

    ! geometric values
    INTEGER :: num_cells,half_cells  ! number of cells per side: total num (n x n)

    REAL(dp_kind) :: Dx, Dy                   ! cell size

    REAL(dp_kind), ALLOCATABLE, DIMENSION(:,:,:) :: x,y,z



END MODULE mod_constants
