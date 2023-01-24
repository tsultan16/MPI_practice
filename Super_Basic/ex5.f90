! Collective Communication #1: Rank 0 will broadcast a number
!                             to all other ranks.

PROGRAM ex5

USE MPI
IMPLICIT NONE

! declare all varialbes and arrays
INTEGER :: ierr, myid, numprocs, itag, irc
INTEGER :: count, root, tag, status(MPI_STATUS_SIZE)
INTEGER :: comm, datatype
REAL :: x,y

! Initialize MPI
CALL MPI_INIT(ierr)
! Get my rank = myid
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
! Get total number of ranks = numprocs
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)


! communication call arguments
count = 1 ! no. of elements in message buffer
tag = 2 ! message tag
root = 0
comm = MPI_COMM_WORLD
datatype = MPI_REAL

x = 0.0 

! set value of x on process 0
IF(myid .EQ. root) THEN
    x = 2.0 
    PRINT*,'Send buffer:', x
    PRINT*,''
END IF

! broadcast value of x to all other processes
CALL MPI_BCAST(x, count, datatype, root, comm, ierr)
         
IF(myid .NE. root)THEN
    y = x**(myid+1)
    PRINT*,'Process #',myid, ', y=',y
END IF
PRINT*,''


! terminte MPI
CALL MPI_FINALIZE(irc)





END PROGRAM ex5
