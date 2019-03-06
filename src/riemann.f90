PROGRAM riemann
    IMPLICIT NONE

    !INTEGER :: i !, j ,k                           ! general iteration variables
    LOGICAL :: test                                 ! general test variable
    CHARACTER(len=32) :: arg                        ! command line argument
    INTEGER :: count                                ! number of CL arguments
    INTEGER, PARAMETER :: n_params = 10             ! number of params in config file
    LOGICAL, DIMENSION(n_params) :: is_parsed = .FALSE.    ! check if every param has been parsed

    ! simulation parameters to be parsed in the config file
    REAL(kind=8) :: T, dt                           ! total time and delta time between writes
    REAL(kind=8), DIMENSION(4) :: p1, p2, p3, p4    ! initial condition in quadrants
    INTEGER :: n                                    ! number of cells per side: total num (n x n)
    INTEGER(kind=4) :: wr                           ! write data to vtk file
    CHARACTER(len=24) :: dirname, outname           ! folder and file where to write results



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

    IF (ALL(is_parsed)) THEN
        PRINT *, CHAR(10),'All arguments have been parsed',CHAR(10)
    ELSE
        PRINT *, CHAR(10),'NOT All arguments have been parsed... Exiting',CHAR(10)
        STOP
    END IF

    ! test if folder exists and create
    INQUIRE(FILE='./'//TRIM(dirname)//'/.', EXIST=test)

    IF (test) THEN
        PRINT *, CHAR(10)//'Directory '//TRIM(dirname)//' exists...'//CHAR(10)
    ELSE
        PRINT *, CHAR(10)//'Directory '//TRIM(dirname)//' does not exist...'//CHAR(10)
        CALL execute_command_line('mkdir -p '//TRIM(dirname))
    END IF


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
                        READ(buffer, *, iostat=ios) T
                        IF (ios == 0) THEN
                            WRITE(*,101) 'Read T:', T
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
                        READ(buffer, *, iostat=ios) n
                        IF (ios == 0) THEN
                            WRITE(*,102) 'Read n:', n
                            is_parsed(3) = .TRUE.
                        ELSE
                            PRINT *, 'Error parsing n'
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
