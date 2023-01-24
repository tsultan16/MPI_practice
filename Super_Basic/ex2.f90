! Basic P2P Communication #1: Rank 0 will send an array
!                             to Rank 1.

PROGRAM ex2

USE MPI
IMPLICIT NONE

! declare all varialbes and arrays
INTEGER :: ierr, myid, numprocs, itag, irc
INTEGER :: count, src, dest, tag, status(MPI_STATUS_SIZE)
INTEGER :: comm, datatype
REAL :: x(5),y, buff

! Initialize MPI
CALL MPI_INIT(ierr)
! Get my rank = myid
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
! Get total number of ranks = numprocs
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)


! communication call arguments
count = 5 ! no. of elements in message buffer
tag = 2 ! message tag
comm = MPI_COMM_WORLD
datatype = MPI_REAL

x = (/ 1., 2., 3., 4., 5. /)
y = 0.0
buff = 0.0

IF(myid .EQ. 0) THEN
    dest = 1
    CALL MPI_SEND(x(1), count, datatype, dest, tag, comm, ierr)
    PRINT*,'Send buffer:', x
ELSE IF(myid .EQ. 1) THEN
    src = 0 
    CALL MPI_RECV(y(1),count, datatype, src, tag, comm, status, ierr) ! note: recv call an an extra 'status' argument (which contains message tag, src rank# and other info) 
    PRINT*,'Received buffer:', y
END IF

PRINT*,''

! terminte MPI
CALL MPI_FINALIZE(irc)





END PROGRAM ex2
