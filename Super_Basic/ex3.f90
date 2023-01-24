! Basic P2P Communication #2: Rank 0 will recieve numbers from
!                             all other ranks, then compute average of those numbers.

PROGRAM ex3

USE MPI
IMPLICIT NONE

! declare all varialbes and arrays
INTEGER :: ierr, myid, numprocs, itag, irc, i
INTEGER :: count, src, des, tag, status(MPI_STATUS_SIZE)
INTEGER :: comm, datatype
REAL :: x(6),y, buff

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

x = (/ 1., 2., 3., 4., 5., 6. /)
y = 0.0
buff = 0.0

IF(numprocs .LT. 7) THEN
    PRINT*,'Error: Need at least 7 ranks...'
    STOP
END IF

IF(myid .EQ. 0) THEN
    ! receive numbers
    DO i = 1, numprocs-1 
        src = i
        CALL MPI_RECV(buff, count, datatype, src, tag, comm, status, ierr)
        y = y+buff
    END DO
    ! compute average
    y = y/REAL(numprocs-1)
    PRINT*,'Rank #:', myid,', Average: ',y
ELSE     
    CALL MPI_SEND(x(myid),count, datatype, des, tag, comm, ierr) ! note: recv call an an extra 'status' argument (which contains message tag, src rank# and other info) 
    PRINT*,'Rank #:', myid,', Message: ',x(myid)
END IF

PRINT*,''

! terminte MPI
CALL MPI_FINALIZE(irc)





END PROGRAM ex3
