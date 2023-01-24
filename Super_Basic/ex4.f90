! Non-blocking communication Communication ex#1: 
!
! Rank 0 will send a number to Rank 1. After calling
! non-blocking send, it will carry on other computations, 
! then it will wait for communication to complete and then 
! overwrite the send buffer and carry out some more computations.                        

PROGRAM ex4

USE MPI
IMPLICIT NONE

! declare all varialbes and arrays
INTEGER :: ierr, myid, numprocs, itag, irc, i
INTEGER :: count, src, des, tag, status(MPI_STATUS_SIZE), request
INTEGER :: comm, datatype
REAL :: x,y,z,a,buff

! Initialize MPI
CALL MPI_INIT(ierr)
! Get my rank = myid
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
! Get total number of ranks = numprocs
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)


! communication call arguments
tag = 1 ! message tag
count = 1
des = 0 
comm = MPI_COMM_WORLD
datatype = MPI_REAL

x = 5.0
y = 0.0
z = 0.0
a = 0.0

IF(myid .EQ. 0) THEN
    des = 1
    ! non-blocking send
    CALL MPI_ISEND(x, count, datatype, des, tag, comm, request, ierr)

    ! compute other stuff while communication takes place
    a = 2.0
    y = a**2

    PRINT*,'Communication in progress. Computing other stuff...' 


    ! wait for communication to complete
    CALL MPI_WAIT(request, status, ierr)

    PRINT*,'Communication complete...'

    ! finish remaining computations after communication is complete
    a = x
    x = x**2

    PRINT*,'a,x,y=',a,x,y
ELSE     
    src = 0
    ! blocking recieve
    CALL MPI_RECV(z,count, datatype, src, tag, comm, status, ierr) !note: non-blocking recv has a request handle which will be used later to block computation until receive completes

    PRINT*,'Rank #',myid,', received z=',z

END IF

PRINT*,''

! terminte MPI
CALL MPI_FINALIZE(irc)





END PROGRAM ex4
