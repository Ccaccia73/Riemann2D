PROGRAM riemann
    IMPLICIT NONE

    !INTEGER :: i !, j ,k          ! general iteration variables
    CHARACTER(len=32) :: arg    ! command line argument
    INTEGER :: count            ! number of CL arguments

    ! simulation parameters to be parsed in the config file
    REAL(kind=8) :: T, dt
    REAL(kind=8), DIMENSION(4) :: p1, p2, p3, p4
    INTEGER :: n

    ! start parsing input file
    count = command_argument_count()

    IF (count < 1 ) THEN
        PRINT*,"count = ",count
        WRITE(*,*) "Wrong number of input parameters... Exiting"
        STOP
    END IF

    CALL get_command_argument(1, arg)

    WRITE (*,*) "Using configuration file: ",TRIM(arg)

    CALL parse_file(arg)


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

                !IF (buffer(1) == '#') THEN
                !    CYCLE
                !END IF
                ! Find the first instance of whitespace.  Split label and data.
                pos = SCAN(buffer, '    ')
                label = buffer(1:pos)
                buffer = buffer(pos+1:)

                SELECT CASE (label)
                    CASE ('T')
                        READ(buffer, *, iostat=ios) T
                        PRINT *, 'Read T: ', T
                    CASE ('dt')
                        READ(buffer, *, iostat=ios) dt
                        PRINT *, 'Read dt: ', dt
                    CASE ('n')
                        READ(buffer, *, iostat=ios) n
                        PRINT *, 'Read n: ', n
                    CASE ('p1')
                        READ(buffer, *, iostat=ios) p1
                        PRINT *, 'Read p1: ', p1
                    CASE ('p2')
                        READ(buffer, *, iostat=ios) p2
                        PRINT *, 'Read vector: ', p2
                    CASE ('p3')
                        READ(buffer, *, iostat=ios) p3
                        PRINT *, 'Read vector: ', p3
                    CASE ('p4')
                        READ(buffer, *, iostat=ios) p4
                        PRINT *, 'Read vector: ', p4
                    CASE DEFAULT
                        PRINT *, 'Skipping invalid label at line', line
                END SELECT
            END IF
        END DO

    END SUBROUTINE parse_file


END PROGRAM riemann
