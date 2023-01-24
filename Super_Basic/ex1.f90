PROGRAM ex1

USE MPI
IMPLICIT NONE

! declare all varialbes and arrays
INTEGER :: ierr, myid, numprocs, itag, irc
REAL :: x,y

! Initialize MPI
CALL MPI_INIT(ierr)
! Get my rank = myid
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
! Get total number of ranks = numprocs
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)


x = 5.0
PRINT*,'x=',x
PRINT*,''

IF(myid .EQ. 0) THEN
    y = x**2
ELSE IF(myid .EQ. 1) THEN
    y = x**3
ELSE IF(myid .EQ. 2) THEN
    y = x**4
END IF

PRINT*,'This is Process#',myid, ', y=',y
PRINT*,''

! terminte MPI
CALL MPI_FINALIZE(irc)

END PROGRAM ex1
