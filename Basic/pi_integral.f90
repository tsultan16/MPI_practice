PROGRAM mpi_pi

! This eaxmple program parallelizes the numerical computation of an integral using MPI.
! The program is run on multiple processes. Process with rank 0 is designated as the "manager",
! the others are designated as "workers". Manager gets value of n (numerical stencil size) from the user,
! then broadcasts it to the workers. The manager and workers compute part of the sum.
! Partial sums from the workers are then collected by the manager and added to get the final result. 

USE MPI
IMPLICIT NONE

REAL(8), PARAMETER :: PI25DT = 3.141592653589793238462643d0 ! pi exact to 25 digits
REAL(8) :: mypi, pi, h, sum, x, f, a, ti, tf
INTEGER :: n, myid, numprocs, i, ierr

! funciton to be integrated
f(a) = 4.d0 / (1.d0 + a*a)

! MPI initialization
CALL MPI_INIT(ierr)

CALL MPI_COMM_RANK(MPI_COMM_WORLD, &  ! (default) communicator handle
                   myid          , &  ! rank of current process within our group
                   ierr          )    ! error code

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, &  ! (default) communicator handle
                   numprocs      , &  ! total number of processes in our group
                   ierr          )    ! error code


IF(ierr .EQ. MPI_SUCCESS) THEN

DO
    IF(myid .EQ. 0)THEN
        !prompt user to enter n
        PRINT*,'Enter n (0 to QUIT):'
        READ(*,*) n
    END IF

    ti = MPI_WTIME()

    ! broadcast n to all worker ranks
    CALL MPI_BCAST(n,               &  ! broadcast message buffer
                   1,               &  ! message buffer size
                   MPI_INTEGER,     &  ! message data type
                   0,               &  ! message source process
                   MPI_COMM_WORLD,  &  ! communicator handle
                   ierr     )          ! error code

    IF(n .LE. 0) EXIT

    ! compute interval size
    h = 1.d0/n
    sum = 0.d0

    ! compute partial sum
    DO i = myid+1, n, numprocs
        x = h * (DBLE(i)-0.5d0)
        sum = sum + f(x)
    END DO
    mypi = sum*h

    ! collect all partial sums in the manager process and add them
    CALL MPI_REDUCE(mypi,                 &  ! source buffer
                    pi,                   &  ! result buffer
                    1,                    &  ! source buffer size
                    MPI_DOUBLE_PRECISION, &  ! data type 
                    MPI_SUM,              &  ! reduction operation
                    0 ,                   &  ! rank of process holding the reduction result
                    MPI_COMM_WORLD,       &  ! communicator handle
                    ierr )
 
    tf = MPI_WTIME()
                   
    ! manager prints the result
    IF(myid .EQ. 0) THEN
        PRINT*,'pi = ',pi,'. Error = ',abs(pi-PI25DT)
        PRINT*,'Number of processes=', numprocs,'. Time elapsed(sec)=',tf-ti
    END IF
END DO

ELSE
    PRINT*,'MPI initialization was unsuccessful.'
END IF



CALL MPI_FINALIZE(ierr)


END PROGRAM mpi_pi
