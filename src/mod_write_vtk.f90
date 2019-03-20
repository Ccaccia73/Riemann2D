MODULE mod_write_vtk

    USE mod_constants
    USE mod_thermodynamics
    USE vtk_fortran, only : vtk_file

    IMPLICIT NONE


CONTAINS
    SUBROUTINE write_vtk(act_w, dirname, prefix, step)

        IMPLICIT NONE
        REAL(dp_kind),DIMENSION(:,:,:),INTENT(IN) :: act_w
        INTEGER, INTENT(IN) :: step
        CHARACTER(len=*), INTENT(IN) :: dirname, prefix

        CHARACTER(len=48) :: filename
        CHARACTER(len=4) :: nstring

        type(vtk_file) :: a_vtk_file                        ! A VTK file.

        !        INTEGER, PARAMETER :: u = 81
        INTEGER :: rstat !, i,j,k

        !Dx = 1./REAL(n)



        WRITE(nstring,'(I4.4)') step

        filename = TRIM(dirname)//'/'//TRIM(prefix)//'_'//nstring//'.vts'


        rstat = a_vtk_file%initialize(format='binary', filename=filename, mesh_topology='StructuredGrid', &
                                      nx1=0, nx2=num_cells, ny1=0, ny2=num_cells, nz1=0, nz2=0)
        rstat = a_vtk_file%xml_writer%write_piece(nx1=0, nx2=num_cells, ny1=0, ny2=num_cells, nz1=0, nz2=0)
        rstat = a_vtk_file%xml_writer%write_geo(n=(num_cells+1)*(num_cells+1), x=x, y=y, z=z)

        rstat = a_vtk_file%xml_writer%write_dataarray(location='cell', action='open')

        ! rho
        rstat = a_vtk_file%xml_writer%write_dataarray(data_name='rho', x=act_w(:,:,i_rho), one_component=.true.)
        ! P
        rstat = a_vtk_file%xml_writer%write_dataarray(data_name='P', x=act_w(:,:,i_P), one_component=.true.)
        ! velocity vector
        rstat = a_vtk_file%xml_writer%write_dataarray(data_name='v', &
                                x=reshape(act_w(:,:,i_u),(/num_cells*num_cells/) ), &
                                y=reshape(act_w(:,:,i_v),(/num_cells*num_cells/) ), &
                                z=0.0*reshape(act_w(:,:,i_u),(/num_cells*num_cells/) ) )
        ! x momentum
        rstat = a_vtk_file%xml_writer%write_dataarray(data_name='mx', x=act_w(:,:,i_mx), one_component=.true.)
        ! y momentum
        rstat = a_vtk_file%xml_writer%write_dataarray(data_name='my', x=act_w(:,:,i_my), one_component=.true.)
        ! total energy
        rstat = a_vtk_file%xml_writer%write_dataarray(data_name='tot energy', x=act_w(:,:,i_eT), one_component=.true.)
        ! internal energy
        rstat = a_vtk_file%xml_writer%write_dataarray(data_name='int energy', x=act_w(:,:,i_e), one_component=.true.)
        ! pseudo schlieren
        rstat = a_vtk_file%xml_writer%write_dataarray(data_name='sound speed', x=act_w(:,:,i_c), one_component=.true.)
        ! pseudo schlieren
        rstat = a_vtk_file%xml_writer%write_dataarray(data_name='schlieren', x=act_w(:,:,i_schl), one_component=.true.)

        rstat = a_vtk_file%xml_writer%write_dataarray(location='cell', action='close')
        rstat = a_vtk_file%xml_writer%write_piece()
        rstat = a_vtk_file%finalize()

        !PRINT*,x
        !PRINT*,y
        !PRINT*,z

    !        IMPLICIT NONE
    !        REAL(kind=8),DIMENSION(:,:,:),INTENT(IN) :: act_w
    !        CHARACTER(len=*), INTENT(IN) :: dirname, prefix
    !        INTEGER, INTENT(IN) :: step
    !        INTEGER, INTENT(IN) :: n
    !        CHARACTER(len=48) :: filename
    !        CHARACTER(len=4) :: nstring
    !
    !        INTEGER, PARAMETER :: u = 81
    !        INTEGER :: rstat, i,j
    !
    !        WRITE(nstring,'(I4.4)') step
    !
    !        filename = TRIM(dirname)//'/'//TRIM(prefix)//'_'//nstring//'.vtk'
    !
    !        OPEN(UNIT=u, FILE=filename, STATUS='REPLACE', ACTION='WRITE', IOSTAT=rstat)
    !        IF (rstat /= 0) THEN
    !            PRINT*, "Errore durante l’apertura del file: ", filename
    !            STOP
    !        END IF
    !
    !        WRITE(u,'(A)') '# vtk DataFile Version 2.0'
    !        WRITE(u,'(A,I4)') 'Riemann 2D step ',step
    !        WRITE(u,'(A)') 'ASCII'
    !        WRITE(u,*) ' '
    !        WRITE(u,'(A)') 'DATASET STRUCTURED_POINTS'
    !        WRITE(u,'(A,3(1X,I10))') 'DIMENSIONS',n+1,n+1,1
    !        WRITE(u,202) 'SPACING',1.0/REAL(n),1.0/REAL(n),1.0/REAL(n)
    !        WRITE(u,202) 'ORIGIN',0.0,0.0,0.0
    !        WRITE(u,*) ' '
    !        WRITE(u,203) 'CELL_DATA', n*n
    !        WRITE(u,204) 'SCALARS','rho','double',1
    !        WRITE(u,'(A)') 'LOOKUP_TABLE default'
    !
    !        i = 1
    !        DO WHILE (i <=n)
    !            WRITE(u,'(*(ES16.6,1X))') act_w(i,:,i_rho)
    !            i = i+1
    !        END DO
    !
    !        WRITE(u,*) ' '
    !        WRITE(u,203) 'CELL_DATA', n*n
    !        WRITE(u,205) 'VECTORS','velocity','double'
    !        i = 1; j = 1
    !        DO WHILE (j<=n)
    !            DO WHILE (i<=n)
    !                WRITE(u,206) act_w(i,j,i_u:i_v),0.0
    !            END DO
    !        END DO
    !
    !        CLOSE(u)
    !
    !!201     FORMAT(A48)
    !202     FORMAT(A, 3(1X, F8.5))
    !203     FORMAT(A,I10)
    !204     FORMAT(3(A,1X),I1)
    !205     FORMAT(2(A,1X),A)
    !206     FORMAT(2(ES16.6,1X),F5.3)
    END SUBROUTINE write_vtk
END MODULE mod_write_vtk
