! Program for testing MPI collective parallel file write. Each MPI rank needs to write 
! a variable amount of data into a common file, which means that for offsets to be computed,
! each rank needs to communicate with all of its preceding ranks to find out their counts, then
! iuse that information to compute the file offset.   


PROGRAM collective_IO

USE MPI

IMPLICIT NONE

!INCLUDE 'mpif.h'



REAL(8), ALLOCATABLE :: my_buff(:)
INTEGER :: myrank, np, nints, buffsize, file_handle, numprocs(1), status(MPI_STATUS_SIZE), ierr
INTEGER(KIND=MPI_OFFSET_KIND) :: offset
CHARACTER(LEN=300) :: filename 
REAL(8) :: t1, t2, t_io, t_all, reduced_size
INTEGER :: i, mysize


CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs(1), ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
np = numprocs(1)


! allocate memory for my data
mysize = 1+myrank
ALLOCATE(my_buff(mysize))

! sample data
DO i = 1, mysize
    my_buff(i) = 10*myrank+i 
END DO

GO TO 110
DO i = 0, np-1 


    CALL SLEEP(1)

    IF(i .EQ. myrank) THEN 

    PRINT*,''
    PRINT*,'#################'
    PRINT*,'Myrank = ',myrank
    PRINT*,'Mydata = ',my_buff
    PRINT*,'#################'
    PRINT*,''

    END IF

    CALL SLEEP(1)

END DO
110 CONTINUE



! perform partial sum of data sizes from rank 0..my rank
CALL MPI_SCAN(DBLE(mysize), reduced_size, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)


! compute the offset
offset = (reduced_size-mysize)*8  ! 8 bytes per data item

GO TO 111
DO i = 0, np-1 


    CALL SLEEP(1)

    IF(i .EQ. myrank) THEN 

    PRINT*,''
    PRINT*,'#################'
    PRINT*,'Myrank = ',myrank
    PRINT*,'Mysize = ',mysize
    PRINT*,'Reduced size = ',INT(reduced_size)
    PRINT*,'File offset = ',offset
    PRINT*,'#################'
    PRINT*,''

    END IF

    CALL SLEEP(1)

END DO
111 CONTINUE




! open parallel file
filename = TRIM('output_parallel.dat')
CALL MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, file_handle, ierr)

! write to file
CALL MPI_FILE_WRITE_AT_ALL(file_handle, offset, my_buff, mysize, MPI_DOUBLE_PRECISION, status, ierr)

CALL MPI_FILE_CLOSE(file_handle, ierr)



CALL MPI_FINALIZE(ierr)


PRINT*,''
PRINT*,'Parallel file write completed for rank: ',myrank
PRINT*,''



CONTAINS





END PROGRAM collective_IO