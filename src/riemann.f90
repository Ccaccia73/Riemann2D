PROGRAM riemann

    USE mod_euler_flux_jacobian
    USE mod_constants
    USE mod_thermodynamics

    USE mod_write_vtk

    IMPLICIT NONE

    INTEGER :: i , j ,k                           ! general iteration variables
    INTEGER :: istat                                !
    LOGICAL :: test                                 ! general test variable
    CHARACTER(len=32) :: arg                        ! command line argument
    INTEGER :: count                                ! number of CL arguments
    INTEGER, PARAMETER :: n_params = 10             ! number of params in config file
    LOGICAL, DIMENSION(n_params) :: is_parsed = .FALSE.    ! check if every param has been parsed

    ! simulation parameters to be parsed in the config file
    REAL(dp_kind) :: t_final, dt                     ! total time and delta time between writes
    REAL(dp_kind), DIMENSION(4) :: p1, p2, p3, p4    ! initial condition in quadrants

    INTEGER(sh_kind) :: wr                           ! write data to vtk file
    CHARACTER(len=24) :: dirname, outname           ! folder and file where to write results

    ! simulation variables
    REAL(dp_kind), ALLOCATABLE, DIMENSION(:,:,:) :: w0   ! matrix containing all variables
    REAL(dp_kind) :: lambda_x, lambda_Y             ! maximum eigenvalue for F and G

    ! time related variables
    REAL(dp_kind) :: t_curr                         ! current time
    REAL(dp_kind) :: dt_cfl                         ! delta time given by CFL condition
    REAL(dp_kind) :: act_dt                         ! actual delta time
    REAL(dp_kind) :: t2next_w                       ! time to next write
    REAL cpu_time_begin, cpu_time_end               ! time elapsed in computations
    LOGICAL :: writestep                            ! write solution file at current step
    INTEGER :: i_step = 0                           ! count step

    ! start parsing input file
    count = command_argument_count()

    IF (count < 1 ) THEN
        PRINT*,'count = ',count
        WRITE(*,*) 'Wrong number of input parameters... Exiting'
        STOP
    END IF

    CALL get_command_argument(1, arg)

    WRITE (*,*) CHAR(10),'Using configuration file: ',TRIM(arg), CHAR(10)
    CALL parse_file(arg)

    ! test if all simulation required parameters have been correctly parsed
    IF (ALL(is_parsed)) THEN
        PRINT *, CHAR(10),'All arguments have been parsed',CHAR(10)
    ELSE
        PRINT *, CHAR(10),'NOT All arguments have been parsed... Exiting',CHAR(10)
        STOP
    END IF

    ! test if folder for results exists and create if necessary
    IF (wr == 1) THEN
        INQUIRE(FILE='./'//TRIM(dirname)//'/.', EXIST=test)

        IF (test) THEN
            PRINT *, CHAR(10)//'Directory '//TRIM(dirname)//' exists...'//CHAR(10)
        ELSE
            PRINT *, CHAR(10)//'Directory '//TRIM(dirname)//' does not exist. Creating...'//CHAR(10)
            CALL execute_command_line('mkdir -p '//TRIM(dirname))
        END IF
    END IF


    ! we want to round up to next even integer (if n is odd)
    IF (num_cells > (num_cells/2)*2 ) THEN
        num_cells = num_cells+1
    END IF
    half_cells = num_cells/2

    ! compute cell dimension
    Dx = 1.0 / REAL(num_cells)
    Dy = Dx

    ! allocate memory for solution
    ALLOCATE (w0(num_cells,num_cells,nvar),STAT=istat)
    IF (istat /= 0) THEN
        PRINT*, "Failed to allocate variables"
        PRINT*, "Error code: ", istat
        PRINT*, "Exiting..."
        STOP
    END IF

    !allocate memory for point coordinatee
    ALLOCATE (x(0:num_cells,0:num_cells,0:0),STAT=istat)
    ALLOCATE (y(0:num_cells,0:num_cells,0:0),STAT=istat)
    ALLOCATE (z(0:num_cells,0:num_cells,0:0),STAT=istat)

    DO k=0, 0
        DO j=0, num_cells
            DO i=0, num_cells
                x(i, j, k) = i*Dx
                y(i, j, k) = j*Dx
            END DO
        END DO
    END DO

    z = 0.0

    ! Intialize known variables
    ! 1
    w0(half_cells+1:num_cells,half_cells+1:num_cells,i_rho) = p1(i_rho)
    w0(half_cells+1:num_cells,half_cells+1:num_cells,i_u)   = p1(i_u)
    w0(half_cells+1:num_cells,half_cells+1:num_cells,i_v)   = p1(i_v)
    w0(half_cells+1:num_cells,half_cells+1:num_cells,i_P)   = p1(i_P)
    ! 2
    w0(1:half_cells,half_cells+1:num_cells,i_rho) = p2(i_rho)
    w0(1:half_cells,half_cells+1:num_cells,i_u)   = p2(i_u)
    w0(1:half_cells,half_cells+1:num_cells,i_v)   = p2(i_v)
    w0(1:half_cells,half_cells+1:num_cells,i_P)   = p2(i_P)
    ! 3
    w0(1:half_cells,1:half_cells,i_rho) = p3(i_rho)
    w0(1:half_cells,1:half_cells,i_u)   = p3(i_u)
    w0(1:half_cells,1:half_cells,i_v)   = p3(i_v)
    w0(1:half_cells,1:half_cells,i_P)   = p3(i_P)
    ! 4
    w0(half_cells+1:num_cells,1:half_cells,i_rho) = p4(i_rho)
    w0(half_cells+1:num_cells,1:half_cells,i_u)   = p4(i_u)
    w0(half_cells+1:num_cells,1:half_cells,i_v)   = p4(i_v)
    w0(half_cells+1:num_cells,1:half_cells,i_P)   = p4(i_P)

    ! write conservative variables
    CALL conservatives(w0)

    ! compute auxiliary variables (internal energy, speed of sound, schieren
    CALL auxiliaries(w0)

    ! write first output file
    IF (wr == 1) THEN
        CALL write_vtk(w0,dirname,outname,i_step)
    END IF

    ! start computation
    t_curr = 0.0d0
    t2next_w = dt
    writestep = .FALSE.

    DO WHILE (t_curr <= t_final)

        ! CFL computation
        lambda_x = max_eig(w0(:,:,i_u),w0(:,:,i_c))
        lambda_y = max_eig(w0(:,:,i_v),w0(:,:,i_c))

        dt_cfl = -1.0

        IF (dt_cfl < t2next_w) THEN
            ! delta time from CFL smaller than time to next write
            writestep = .FALSE.
            act_dt = dt_cfl
            t2next_w = t2next_w - act_dt
        ELSE
            ! delta time from cfl bigger than time to next write
            writestep = .TRUE.
            act_dt = t2next_w
            t2next_w = dt
        END IF
        ! update time
        t_curr = t_curr + act_dt
        CALL CPU_TIME(cpu_time_begin)

        ! perform some veeeeeeeery difficult computation
        ! CALL SLEEP(2)
        ! TODO





        IF (wr == 1 .and. writestep ) THEN
            i_step = i_step + 1
            CALL write_vtk(w0,dirname,outname,i_step)
        END IF

        CALL CPU_TIME(cpu_time_end)

        WRITE(*,111) "t: ", t_curr, "Execution time: ", cpu_time_end - cpu_time_begin, "sec"

    END DO


    DEALLOCATE(w0)

    PRINT*,"End of computation!"

111 FORMAT(1X, A3, 1X, ES13.3, 1X, A16, ES13.3, 1X, A3)


CONTAINS

    SUBROUTINE parse_file(filename)
        IMPLICIT NONE
        CHARACTER(LEN=32),INTENT(IN) :: filename
        CHARACTER(len=100) :: buffer, label
        INTEGER :: pos
        INTEGER, PARAMETER :: fh = 15
        INTEGER :: ios = 0
        INTEGER :: line = 0

        OPEN(fh, file=filename)

        ! ios is negative if an end of record condition is encountered or if
        ! an endfile condition was detected.  It is positive if an error was
        ! detected.  ios is zero otherwise.

        DO WHILE (ios == 0)
            READ(fh, '(A)', iostat=ios) buffer
            IF (ios == 0) THEN
                line = line + 1

                IF (buffer(1:1) == '#') THEN
                    PRINT*, "comment: ",buffer
                    CYCLE
                END IF
                ! Find the first instance of whitespace.  Split label and data.
                pos = SCAN(buffer, '    ')
                label = buffer(1:pos)
                buffer = buffer(pos+1:)

                SELECT CASE (label)
                    CASE ('T')
                        READ(buffer, *, iostat=ios) t_final
                        IF (ios == 0) THEN
                            WRITE(*,101) 'Read T:', t_final
                            is_parsed(1) = .TRUE.
                        ELSE
                            PRINT *, 'Error parsing T'
                        END IF

                    CASE ('dt')
                        READ(buffer, *, iostat=ios) dt
                        IF (ios == 0) THEN
                            WRITE(*,101) 'Read dt:', dt
                            is_parsed(2) = .TRUE.
                        ELSE
                            PRINT *, 'Error parsing dt'
                        END IF

                    CASE ('n')
                        READ(buffer, *, iostat=ios) num_cells
                        IF (ios == 0) THEN
                            WRITE(*,102) 'Read num cells:', num_cells
                            is_parsed(3) = .TRUE.
                        ELSE
                            PRINT *, 'Error parsing num cells'
                        END IF

                    CASE ('p1')
                        READ(buffer, *, iostat=ios) p1
                        IF (ios == 0) THEN
                            WRITE(*,103) 'Read p1:', p1
                            is_parsed(4) = .TRUE.
                        ELSE
                            PRINT *, 'Error parsing p1'
                        END IF

                    CASE ('p2')
                        READ(buffer, *, iostat=ios) p2
                        IF (ios == 0) THEN
                            WRITE(*,103) 'Read p2:', p2
                            is_parsed(5) = .TRUE.
                        ELSE
                            PRINT *, 'Error parsing p2'
                        END IF

                    CASE ('p3')
                        READ(buffer, *, iostat=ios) p3
                        IF (ios == 0) THEN
                            WRITE(*,103) 'Read p3:', p3
                            is_parsed(6) = .TRUE.
                        ELSE
                            PRINT *, 'Error parsing p3'
                        END IF

                    CASE ('p4')
                        READ(buffer, *, iostat=ios) p4
                        IF (ios == 0) THEN
                            WRITE(*,103) 'Read p4:', p4
                            is_parsed(7) = .TRUE.
                        ELSE
                            PRINT *, 'Error parsing p4'
                        END IF

                    CASE ('write')
                        READ(buffer, *, iostat=ios) wr
                        IF (ios==0) THEN
                            WRITE(*,102) 'Read write condition:', wr
                            is_parsed(8) = .TRUE.
                            IF (wr == 0) THEN
                                is_parsed(9:10) = .TRUE.
                            ELSEIF (wr /= 1) THEN
                                wr = 1
                            END IF
                        ELSE
                            PRINT *, 'Error parsing write'
                        END IF
                    CASE ('outdir')
                        READ(buffer, *, iostat=ios) dirname
                        IF (ios==0) THEN
                            WRITE(*, 104) 'Read output directory:', dirname
                            is_parsed(9) = .TRUE.
                        ELSE
                            PRINT *, 'Error parsing output directory'
                        END IF

                    CASE ('outfile')
                        READ(buffer, *, iostat=ios) outname
                        IF (ios==0) THEN
                            WRITE(*, 104) 'Read output filename:', outname
                            is_parsed(10) = .TRUE.
                        ELSE
                            PRINT *, 'Error parsing output filename'
                        END IF

                    CASE DEFAULT
                        PRINT *, 'Skipping invalid label at line', line
                END SELECT
            END IF
        END DO

101     FORMAT(1X, A25, 3X, ES13.3)
102     FORMAT(1X, A25, 3X, I5)
103     FORMAT(1X, A25, 3X, 4(ES16.6, 2X))
104     FORMAT(1X, A25, 3X, A24)
    END SUBROUTINE parse_file


END PROGRAM riemann
