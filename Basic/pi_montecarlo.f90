PROGRAM mpi_pi_com

! This eaxmple program parallelizes the monte carlo computation of pi using MPI.
! The program is run on multiple processes. One of the processes is designated as the "server",
! the others are designated as "workers". A custom communicator will be defined for the group of 
! workers only (i.e. excluding the server). The server generates random numbers and sends these to the workers. 
! Each worker performs Monte Carlo sampling of it's random numbers. Sampling results from all workers are then 
! collected, added and used to compute pi. If the error exceeds the tolerance, then workers request more 
! random numbers from the server and carry out further iterations. 

USE MPI
IMPLICIT NONE

INTEGER :: maxnums = 1000000
REAL(8), PARAMETER :: PI25DT = 3.141592653589793238462643d0 ! pi exact to 25 digits
INTEGER, PARAMETER :: REQUEST = 1, REPLY = 2 ! message tags
INTEGER, PARAMETER :: CHUNKSIZE = 1000
REAL(8) :: x, y, pi, error, epsilon, rands(CHUNKSIZE) 
REAL(8) :: ti, tf
REAL :: p
INTEGER :: server, myid, numprocs, ierr, status(MPI_STATUS_SIZE)
INTEGER :: workerid, totalin, totalout
INTEGER :: iter, in, out, i, iters, max, ix, iy, ranks(1), done, temp
INTEGER :: request_msg
! custom communicator variables
INTEGER :: world, workers
INTEGER :: world_group, worker_group


! MPI initialization
CALL MPI_INIT(ierr)

world = MPI_COMM_WORLD  ! designate 'world' as the default communicator 

CALL MPI_COMM_RANK(world,          &  !  communicator handle
                   myid,           &  ! rank of current process within our group
                   ierr          )    ! error code

CALL MPI_COMM_SIZE(world,          &  ! (default) communicator handle
                   numprocs,       &  ! total number of processes in our group
                   ierr          )    ! error code

! designate the last process in the group as the server
server = numprocs - 1

PRINT*,''

IF(ierr .EQ. MPI_SUCCESS) THEN

    IF(myid .EQ. 0) THEN
        ! set error tolerance
        epsilon = 1.0e-12    
    END IF

    ! broadcast epsilon to all other processes in world group
    CALL MPI_BCAST(epsilon,               &  ! address of first message buffer element 
                   1,                     &  ! message count
                   MPI_DOUBLE_PRECISION,  &  ! message data type
                   0,                     &  ! message source process
                   world, ierr) 

    ! form new worker group
    CALL MPI_COMM_GROUP(world, world_group, ierr) ! extract group of processes from 'world' communicator
    ranks(1) = server
    CALL MPI_GROUP_EXCL(world_group, 1, ranks, worker_group, ierr) ! make a new group 'worker_group' that excludes the server process
    CALL MPI_COMM_CREATE(world, worker_group, workers, ierr) ! create communicator handle for this new group
    CALL MPI_GROUP_FREE(worker_group, ierr) ! now that the ne communicator has been created, free up the group 

    ti = MPI_WTIME()

    IF(myid .EQ. server) THEN  ! server process

        DO 
            ! receive request message from workers    
            CALL MPI_RECV(request_msg,    &  ! message buffer
                          1,              &  ! message count
                          MPI_INTEGER,    &  ! data type
                          MPI_ANY_SOURCE, &  ! wildcard source
                          REQUEST,        &  ! message tag
                          world,          & 
                          status,         &  
                          ierr)  
            IF(request_msg .EQ. 1)THEN
                ! generate random numbers and send them
                DO i = 1, CHUNKSIZE
                    CALL RANDOM_NUMBER(p)
                    rands(i) = p
                END DO

                CALL MPI_SEND(rands,                 &  ! message buffer
                              CHUNKSIZE,             &  ! message count
                              MPI_DOUBLE_PRECISION,  &  ! data type
                              status(MPI_SOURCE),    &  ! destination
                              REPLY,                 &  ! message tag
                              world, ierr)   
            ELSE 
                EXIT
            END IF
        END DO

    ELSE  ! worker process
        request_msg = 1
        done = 0
        in = 0
        out = 0
        
        ! request random numbers form server
        CALL MPI_SEND(request_msg,       &  ! message buffer
                 1,             &  ! message count
                 MPI_INTEGER,   &  ! data type
                 server,        &  ! destination
                 REQUEST,       &  ! message tag
                 world, ierr)   
        
        ! get current worker rank
        CALL MPI_COMM_RANK(workers,   &  !  communicator handle
                           workerid,  &  ! rank of current process within our group
                           ierr)         ! error code

        iter = 0
        DO WHILE (done .NE. 1)
            iter = iter + 1   
            request_msg = 1         
            ! receive random numbers from server
            CALL MPI_RECV(rands,                 &  ! message buffer
                          CHUNKSIZE,             &  ! message count
                          MPI_DOUBLE_PRECISION,  &  ! data type
                          server,                &  ! source
                          REPLY,                 &  ! tag 
                          world,                 & 
                          MPI_STATUS_IGNORE,     &  ! ignore status information
                          ierr) 
  
            DO i  = 1, CHUNKSIZE, 2
                x = rands(i)
                y = rands(i+1)
                IF( x*x+y*y .LT. 1.d0)THEN
                    in = in + 1
                ELSE
                    out = out + 1
                END IF
            END DO

            ! sum up 'in' and 'out' values from all servers and make the result available to them
            CALL MPI_ALLREDUCE(in,             &  ! source buffer
                               totalin,        &  ! result buffer
                               1,              &  ! source buffer size
                               MPI_INTEGER,    &  ! data type 
                               MPI_SUM,        &  ! reduction operation
                               workers,        &  ! communicator handle
                               ierr )

            CALL MPI_ALLREDUCE(out,            &  ! source buffer
                               totalout,       &  ! result buffer
                               1,              &  ! source buffer size
                               MPI_INTEGER,    &  ! data type 
                               MPI_SUM,        &  ! reduction operation
                               workers,        &  ! communicator handle
                               ierr )
 
            CALL MPI_ALLREDUCE(iter,           &  ! source buffer
                               iters,          &  ! result buffer
                               1,              &  ! source buffer size
                               MPI_INTEGER,    &  ! data type 
                               MPI_SUM,        &  ! reduction operation
                               workers,        &  ! communicator handle
                               ierr )
 
            ! compute pi
            pi = 4.d0*totalin/DBLE(totalin+totalout)
            error = abs(pi-PI25DT)

            ! print out current value of pi
            !IF(myid .EQ. 0) THEN
            !    PRINT*,'Current pi,error = ',pi,error
            !END IF

            ! send more requests if error exceeds the tolerance
            IF(error .LT. epsilon .OR. (totalin+totalout) .GT. maxnums) THEN
                request_msg = 0 
                done = 1
            ELSE
                request_msg = 1
                done  = 0
            END IF

            CALL MPI_SEND(request_msg,  &  ! message buffer
                     1,                 &  ! message count
                     MPI_INTEGER,       &  ! data type
                     server,            &  ! destination
                     REQUEST,           &  ! message tag
                     world, ierr)   
           

        END DO

        CALL MPI_COMM_FREE(workers, ierr)

    END IF

    tf = MPI_WTIME()

ELSE 
    PRINT*,'MPI initialization was unsuccessful..'
END IF

IF(myid .EQ. 0) THEN
    PRINT*,'Total in, Total out=',totalin, totalout
    PRINT*,'Total Iterations = ',iters
    PRINT*,'pi,error =',pi,error
    PRINT*,'Time elapsed(s) = ',tf-ti
END IF


CALL MPI_FINALIZE(ierr)



END PROGRAM mpi_pi_com
