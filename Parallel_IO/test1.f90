! Parallel I/O Test 1:
! Each MPI ranks reads a portion from a common file.

PROGRAM test1

IMPLICIT NONE

INCLUDE 'mpif.h'


INTEGER, PARAMETER :: filesize = 20 ! number of items in the file

INTEGER :: myrank, np, nints, buffsize, file_handle, numprocs(1), status(MPI_STATUS_SIZE), ierr, i
INTEGER(KIND=MPI_OFFSET_KIND) :: offset
REAL*8, ALLOCATABLE :: buff(:)
CHARACTER(LEN=100) :: filename 


CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs(1), ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
np = numprocs(1)


! each file will read an equal chink from the file
buffsize = filesize/np
ALLOCATE(buff(buffsize))

! open parallel file
filename = TRIM('output.dat')
CALL MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, ierr)

! set the file offset (i.e. number of bytes to skip from the start of the file) 
offset = myrank*buffsize*8  ! 8 bytes per data item
CALL MPI_FILE_READ_AT(file_handle, offset, buff, buffsize, MPI_DOUBLE_PRECISION, status, ierr)

CALL MPI_FILE_CLOSE(file_handle, ierr)

DO i = 0, np-1
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
    IF(myrank .EQ. i) THEN
        PRINT*,''
        PRINT*,'Myrank =',myrank
        PRINT*,'Data read from file: ',buff
        PRINT*,''
        CALL SLEEP(1)
    END IF    
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
END DO

IF(myrank .EQ. 0) PRINT*,'Done!'


CALL MPI_FINALIZE(ierr)


CONTAINS

END PROGRAM test1
