PROGRAM workshare

USE MPI
IMPLICIT NONE


INTEGER :: this_process, n_process, err_num, n 
INTEGER :: status(MPI_STATUS_SIZE)
INTEGER, ALLOCATABLE :: x(:)
INTEGER, PARAMETER :: factor = 5
INTEGER :: i,j,k
INTEGER :: start, end, recv_start


! initialize mpi
CALL MPI_INIT(err_num)  
! get total # of proceses
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, n_process, err_num)
! get # of this MPI rank
CALL MPI_COMM_RANK(MPI_COMM_WORLD, this_process, err_num)


! each process will work on part of the array x (# of elements per rank = factors)
n = factor*n_process
ALLOCATE(x(1:n))
x = 0 

! specify the portion of the array that this process will work on
start = 1+this_process*factor
end = (1+this_process)*factor

PRINT*,'Rank# =',this_process
PRINT*,'start index, end index=', start,end
PRINT*,''

! work on the designated portion of the array
DO  i = start, end
    x(i) = i*factor
END DO


! print array contents
DO i=1,n
    PRINT*,'i,x(i)=',i,x(i) 
END DO
PRINT*,''

IF(this_process .EQ. 0) THEN
    ! have root rank receive array contents from the other worker ranks
    DO i = 1, n_process-1
        recv_start = 1+factor*i ! starting index of worker rank's array portion
        CALL MPI_RECV(x(recv_start),factor,MPI_INTEGER,i,1,MPI_COMM_WORLD,status, err_num)
    END DO
ELSE
    ! have worker ranks send their array contents to the root rank
    CALL MPI_SEND(x(start),factor, MPI_INTEGER,0,1,MPI_COMM_WORLD,err_num)

END IF


! print final array contents
IF(this_process .EQ. 0) THEN
    ! have root rank receive array contents from the other worker ranks
    DO i = 1, n
        PRINT*,'i,x(i)=',i,x(i) 
    END DO
    PRINT*,''
END IF

CALL MPI_FINALIZE(err_num)


END PROGRAM workshare