! Parallel I/O Test 2:
! Each MPI ranks writes a portion to a common file.

PROGRAM test2

IMPLICIT NONE

INCLUDE 'mpif.h'


INTEGER, PARAMETER :: filesize = 20 ! number of items in the file

INTEGER :: myrank, np, nints, buffsize, file_handle, numprocs(1), status(MPI_STATUS_SIZE), ierr, i
INTEGER(KIND=MPI_OFFSET_KIND) :: offset
REAL*8, ALLOCATABLE :: buff(:)
CHARACTER(LEN=100) :: filename 
REAL*8 :: t1, t2, t_io, t_all

CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs(1), ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
np = numprocs(1)


! each file will read an equal chink from the file
buffsize = filesize/np
ALLOCATE(buff(buffsize))

DO i = 1, buffsize
    buff(i) = myrank*buffsize + i 
END DO

! open parallel file
filename = TRIM('output2.dat')

t1 = MPI_WTIME()
CALL MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, file_handle, ierr)

! set the file offset (i.e. number of bytes to skip from the start of the file) 
offset = myrank*buffsize*8  ! 8 bytes per data item
CALL MPI_FILE_WRITE_AT(file_handle, offset, buff, buffsize, MPI_DOUBLE_PRECISION, status, ierr)

CALL MPI_FILE_CLOSE(file_handle, ierr)
t2 = MPI_WTIME()

CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

t_io = t2-t1

CALL MPI_ALLREDUCE(t_io, t_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)



IF(myrank .EQ. 0) THEN

    PRINT*,'Average Parallel file write time per rank (ms) = ',1e3*t_all/np
    PRINT*,'Done!'

END IF

CALL MPI_FINALIZE(ierr)


CONTAINS

END PROGRAM test2
