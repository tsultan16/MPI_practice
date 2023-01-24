PROGRAM mpi_matvec

! This eaxmple program parallelizes the multiplication of a vector 'b' by a matrix 'A' using MPI.
! The program is run on multiple processes. Process with rank 0 is designated as the "manager",
! the others are designated as "workers". Manager first broadcasts the vector 'b' to all the workers. 
! The manager then sends each worker one row from the matrix 'A'. Then start a loop that runs over 
! the rows. Each worker computes the dot product of it's row of A with b, then sends the result back to the manager.
! Whenever the manager receives a dot product from a worker, it sends the next row to that worker. The loop ends when
! dot products for all rows have been received. The manager then sends out termination messages to the workers.
! In side each worker process, a loop is initiated. Inside this loop, rows are received from the manager, and dot product
! is computed and the result send back to manager. This loop terminates upon receiving the termination message.

USE MPI
IMPLICIT NONE

INTEGER, PARAMETER :: MAX_ROWS = 10000, MAX_COLS = 10000
INTEGER :: rows, cols 
REAL(8) :: A(MAX_ROWS,MAX_COLS), b(MAX_COLS) 
REAL(8) :: c(MAX_ROWS), buffer(MAX_COLS)
REAL(8) :: ti, tf
INTEGER :: manager, myid, numprocs, tag, dest, ierr, status(MPI_STATUS_SIZE)
INTEGER :: i, j, numsent, sender
INTEGER :: anstype, numrow, ans


! MPI initialization
CALL MPI_INIT(ierr)

CALL MPI_COMM_RANK(MPI_COMM_WORLD, &  ! (default) communicator handle
                   myid          , &  ! rank of current process within our group
                   ierr          )    ! error code

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, &  ! (default) communicator handle
                   numprocs      , &  ! total number of processes in our group
                   ierr          )    ! error code

manager = 0
rows = 10000
cols = 10000

PRINT*,''

IF(ierr .EQ. MPI_SUCCESS .AND. numprocs .GT. 1) THEN

    ! manager broadcasts vector b, then sends matrix A rows to workers
    IF(myid .EQ. manager) THEN
        ! first, initialize vector b and matrix A
        DO j = 1, cols
            b(j) = 1
            DO i = 1, rows
                A(i,j) = i
            END DO
        END DO

        ti = MPI_WTIME()

        ! broadcast vector b to all workers
        CALL MPI_BCAST(b,                     &  ! message buffer 
                       cols,                  &  ! message count
                       MPI_DOUBLE_PRECISION,  &  ! message data type
                       manager,               &  ! message source process
                       MPI_COMM_WORLD, ierr)

        numsent = 0

        ! send a row to each worker rank, tag with row number
        DO i = 1, MIN(rows, numprocs-1)
            buffer(:) = A(i,:)
            tag = i
            dest = i
            CALL MPI_SEND(buffer,                &  ! message buffer
                          cols,                  &  ! message count
                          MPI_DOUBLE_PRECISION,  &  ! data type
                          dest,                  &  ! destination
                          tag,                   &  ! message tag (row number)
                          MPI_COMM_WORLD, ierr)        
            numsent = numsent + 1    
        END DO

        ! receive dot products from workers
        DO  i= 1, rows 
            CALL MPI_RECV(ans,                   &  ! message buffer
                          1,                     &  ! message count
                          MPI_DOUBLE_PRECISION,  &  ! data type
                          MPI_ANY_SOURCE,        &  ! wildcard source
                          MPI_ANY_TAG,           &  ! wildcard tag (row number)
                          MPI_COMM_WORLD,        & 
                          status,                &  ! contains information about message (source, tag, etc..)
                          ierr)  

            ! find out message source and tag 
            sender = status(MPI_SOURCE)
            numrow = status(MPI_TAG)  ! tag value is the row
            
            ! store received dot product in result array
            c(numrow) = ans
            
            ! loop over remaining rows (if any) and send to next avilable worker (i.e. the most recent worker to send a dot product)
            IF(numsent .LT. rows)Then
                buffer(:) = A(numsent + 1,:)
                dest = sender
                tag = numsent + 1
                CALL MPI_SEND(buffer,                &  ! message buffer
                              cols,                  &  ! message count
                              MPI_DOUBLE_PRECISION,  &  ! data type
                              dest,                  &  ! destination
                              tag,                   &  ! message tag (row number)
                              MPI_COMM_WORLD, ierr)        
                numsent = numsent + 1 
           
            ! send termination message to worker (message of count=0 and tag=0)
            ELSE
                dest = sender
                tag = 0
                CALL MPI_SEND(MPI_BOTTOM,            &  ! message buffer
                              0,                     &  ! message count
                              MPI_DOUBLE_PRECISION,  &  ! data type
                              dest,                  &  ! destination
                              tag,                   &  ! message tag (row number)
                              MPI_COMM_WORLD, ierr)        
                
            END IF
 
        END DO

        tf = MPI_WTIME()

        ! print the final result
        !PRINT*,'A*b = '
        !DO i = 1, cols
        !    PRINT*,'  ',c(i)
        !END DO
        PRINT*,''
        PRINT*,'Number of processes = ',numprocs,'. Time elapsed(s) = ',tf-ti
        PRINT*,''

    ! workers receive b and rows, compute dot product and send back to manager, until they get termination message
    ELSE
        CALL MPI_BCAST(b,                     &  ! message buffer 
                       cols,                  &  ! message count
                       MPI_DOUBLE_PRECISION,  &  ! message data type
                       manager,               &  ! message source process
                       MPI_COMM_WORLD, ierr)

        ! skip if more processes than work
        IF(myid .LE. rows) THEN
            DO
                CALL MPI_RECV(buffer,                &  ! message buffer
                              cols,                  &  ! message count
                              MPI_DOUBLE_PRECISION,  &  ! data type
                              manager,               &  ! source is manager
                              MPI_ANY_TAG,           &  ! wildcard message tag (row number)
                              MPI_COMM_WORLD,        & 
                              status,                &  ! contains information about message (source, tag, etc..)
                              ierr)

                numrow = status(MPI_TAG)

                IF(numrow .EQ. 0) EXIT
                
                ! compute dot product
                ans = 0.0
                DO i = 1, cols
                    ans = ans + b(i) * buffer(i)
                END DO

                ! send result back to manager
                dest = manager
                tag = numrow
                CALL MPI_SEND(ans,                   &  ! message buffer
                              1,                     &  ! message count
                              MPI_DOUBLE_PRECISION,  &  ! data type
                              dest,                  &  ! destination
                              tag,                   &  ! message tag (row number)
                              MPI_COMM_WORLD, ierr)        
              

            END DO

        END IF

    END IF

ELSE 
    IF(ierr .NE. MPI_SUCCESS) THEN
        PRINT*,'MPI initialization was unsuccessful.'
    END IF
   
    IF(numprocs .LT. 2) THEN
        PRINT*,'Need at least two processes.'
    END IF
    PRINT*,''
   
END IF



CALL MPI_FINALIZE(ierr)


END PROGRAM mpi_matvec
