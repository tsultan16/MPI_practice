! Overview:  This program numerically solves the 2D Poisson's equation (nabla(u(x,y)) = f(x,y)) using a 5-pt finite 
! difference approximation and Jacobi iterations on an (N+2) by (N+2) square grid: 
!
!      u_i,j^(k+1) = (1/4)*(u_i+1,j^(k) + u_i-1,j^(k) + u_i,j+1^(k) + u_i,j-1^(k) - h^2 * f_i,j ) 
!
! 1D Domain decomposition:  The 2d computational grid is decomposed into equally sized chunks.
! The decomposition is 1d (the chunks are cut out along the y direction).
! A process is assigned to each chunk. We use an MPI_CART communicator to
! to define the decomposition geometry and handle communication between processes.  
! Each sub-domain has a layer of boundary/ghost points which is 1-point deep.
! Ghost point data is exchanged between neighboring processes. 
!
! To ensure that processes can perform in parallel, we can use: 
! 1) ordered blocking sends and recvs(even numbered ranks send first, odd ones recv first). 
! 2) combined blocking sendrecvs. 
! 3) non-blocking sends and recvs.
!
! Note: Using unproperly ordered blocking sends/recvs can result in deadlock or even cause the process
! to run sequentially.

PROGRAM poisson_1

USE MPI
IMPLICIT NONE


INTEGER, PARAMETER :: exchange_option = 3  ! 1: separate send-recv 2: combined sendrecv 3: non-blocking
INTEGER :: comm1d, ierr, ndim, myrank, numprocs(1), coord(1), req(4), status(MPI_STATUS_SIZE), source, tag
LOGICAL :: isperiodic(1), reorder 
INTEGER, PARAMETER :: N = 126, nb = 1, maxiter = 1000
REAL(8), PARAMETER :: tol = 1.d-4
REAL(8) :: u(1:N+2,1:N+2), h, diff, difftot, ti, tf
REAL(8), ALLOCATABLE :: u_old(:,:),u_new(:,:), f(:,:), u_flat(:), &
                        buffer(:), sbuffer(:), rbuffer(:), &
                        sbuffertop(:), rbuffertop(:), sbufferbot(:), rbufferbot(:)
INTEGER :: ilow, ihi, jlow, jhi, i, j, iter, np
INTEGER :: bottom, top


! MPI initialization
CALL MPI_INIT(ierr)

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs(1), ierr)  

np = numprocs(1)

! create a new cartesian communicator for our 2d domain decomposition

isperiodic(1) = .FALSE. ! periodic boundaries
reorder = .TRUE. ! allow MPI rank reordering 
ndim = 1 ! 1d decomposition


CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndim, numprocs, isperiodic, reorder, comm1d, ierr)

CALL MPI_COMM_RANK(comm1d, myrank, ierr)

! get rank coordinate
CALL MPI_CART_COORDS(comm1d, myrank, ndim, coord, ierr) 

!PRINT*,''
!PRINT*,'My rank, coordinate, num processes = ',myrank, coord, numprocs
!PRINT*,''

! compute boundary indices for this rank
CALL compute_bound(ilow, ihi, jlow, jhi)

h = 1.d0/(N+1)

! initialize u and f
IF(exchange_option .EQ. 1 .OR. exchange_option .EQ. 2)THEN
    ALLOCATE(buffer(ilow:ihi))
    ALLOCATE(sbuffer(ilow:ihi),rbuffer(ilow:ihi))
ELSE IF (exchange_option .EQ. 3)THEN
    ALLOCATE(sbuffertop(ilow:ihi),rbuffertop(ilow:ihi))
    ALLOCATE(sbufferbot(ilow:ihi),rbufferbot(ilow:ihi))
END IF

ALLOCATE(u_old(ilow:ihi, jlow-nb:jhi+nb),u_new(ilow:ihi, jlow-nb:jhi+nb))
ALLOCATE(f(ilow:ihi, jlow-nb:jhi+nb))
ALLOCATE(u_flat(1:(N+2)*(N+2)))

CALL init_1d(u_old, f, ilow, ihi, jlow, jhi)
u_new = u_old
u = 0.0

OPEN(UNIT = 10, FILE = 'output.txt')
OPEN(UNIT = 11, FILE = 'horcut.txt')
OPEN(UNIT = 12, FILE = 'vercut.txt')

CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
ti = MPI_WTIME()

! Jacobi iteration loop
iter = 0
DO WHILE(iter .LT. maxiter)

    !PRINT*,'Rank # = ',myrank,' Iteration #',iter
    ! exchange ghost point data with neighbor ranks 
    !PRINT*,'Rank # = ',myrank,' Starting Exchange.'
    IF(exchange_option .EQ. 1) THEN
        CALL exchange_data_1(u_old, sbuffer, rbuffer, ilow, ihi, jlow, jhi)
    ELSE IF(exchange_option .EQ. 2) THEN
        CALL exchange_data_2(u_old, sbuffer, rbuffer, ilow, ihi, jlow, jhi)
    ELSE IF(exchange_option .EQ. 3) THEN
        CALL exchange_data_3(u_old, sbuffertop, rbuffertop, sbufferbot, rbufferbot, &
                             ilow, ihi, jlow, jhi)
    END IF
    !PRINT*,'Rank # = ',myrank,' Exchange Complete.'

    ! update u
    CALL sweep(u_old, u_new, ilow, ihi, jlow, jhi)

    !PRINT*,'Rank # = ',myrank,' Update Complete.'

    ! check convergence
    CALL compute_difference(u_old, u_new, diff, ilow, ihi, jlow, jhi)   

    !PRINT*,'Rank # = ',myrank,' Difference = ',diff

    ! add up difference from all ranks
    CALL MPI_ALLREDUCE(diff,                  &  ! source buffer
                       difftot,               &  ! result buffer
                       1,                     &  ! source buffer size
                       MPI_DOUBLE_PRECISION,  &  ! data type 
                       MPI_SUM,               &  ! reduction operation
                       comm1d,                &  ! communicator handle
                       ierr )

    !PRINT*,'Total difference across all ranks = ',difftot

    ! update u_old
    u_old = u_new

    iter = iter + 1

    IF(difftot .LE. tol) EXIT
   
END DO

CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
tf = MPI_WTIME()

!PRINT*,''
!PRINT*,'Rank #',myrank, ' Done.'
!PRINT*,'Num processes = ', numprocs
!PRINT*,''

IF(myrank .EQ. 0) THEN
    PRINT*,''
    PRINT*,'Iterations = ',iter
    PRINT*,'Average difference = ',difftot
    PRINT*,'Time elapsed(s) = ',tf-ti
    PRINT*,''
END IF


!PRINT*,'Rank# ',myrank,'u_new = '
!DO j = jhi+1,jlow-1,-1
!    DO i = ilow, ihi
!        WRITE(*,FMT='(1f8.3)',ADVANCE='no') u_new(i,j)
!    END DO
!    PRINT*,''
!END DO

20  FORMAT(1i4)


! I/O: Send flattened array back to rank 0
CALL prepare_u(u_new)

IF(myrank .EQ. 0) THEN
    DO i = 1, N+2
        DO j = 1, N+2
            WRITE(10,*) i*h, j*h, u(i,j)
        END DO
        WRITE(11,*) i*h, u(i,(N+2)/2)
        WRITE(12,*) i*h, u((N+2)/2,i)
    END DO
END IF

IF(exchange_option .EQ. 1 .OR. exchange_option .EQ. 2)THEN
    DEALLOCATE(buffer)
    DEALLOCATE(sbuffer, rbuffer)
ELSE IF(exchange_option .EQ. 3)THEN
    DEALLOCATE(sbuffertop, rbuffertop, sbufferbot, rbufferbot)
END IF 
DEALLOCATE(u_old, u_new, u_flat)

CLOSE(UNIT = 10)
CLOSE(UNIT = 11)
CLOSE(UNIT = 12)

! MPI termination
CALL MPI_FINALIZE(ierr)


CONTAINS


SUBROUTINE compute_bound(ilow, ihi, jlow, jhi)

    INTEGER, INTENT(INOUT) :: ilow, ihi, jlow, jhi

    ilow = 0
    ihi = N+1
    jlow = 0+coord(1)*(N+2)/np
    jhi = jlow-1+(N+2)/np

END SUBROUTINE compute_bound


SUBROUTINE prepare_u(u_new)

    REAL(8), INTENT(IN) :: u_new(:,:)
    INTEGER :: i, j

    IF(np .GT. 1) THEN 
        IF(myrank .NE. 0) THEN
            tag = coord(1)
            CALL flatten_array(u_new, u_flat)
            CALL MPI_SEND(u_flat, (N+2)*(N+2), MPI_DOUBLE_PRECISION, 0, tag, comm1d, ierr)
           
            !PRINT*,'Rank #',myrank,'u_my ='           
            !DO j = 2+(N+2)/np,1,-1
            !    DO i = 1, N+2
            !        WRITE(*,FMT='(1f8.3)',ADVANCE='no') u_new(i,j)
            !    END DO
            !    PRINT*,''
            !END DO
            !PRINT*,''

        ELSE

            !PRINT*,'Rank #',myrank,'u_my ='
            !DO j = 2+(N+2)/np,1,-1
            !    DO i = 1, N+2
            !        WRITE(*,FMT='(1f8.3)',ADVANCE='no') u_new(i,j)
            !    END DO
            !    PRINT*,''
            !END DO
            !PRINT*,''

            CALL flatten_array(u_new, u_flat)
            CALL unflatten_save(u_flat, 0) 
                
            DO i= 1, np-1
                tag = i
                source = i
                CALL MPI_RECV(u_flat, (N+2)*(N+2), MPI_DOUBLE_PRECISION, source, tag, comm1d, status, ierr)                
                CALL unflatten_save(u_flat, tag)
            END DO

            !PRINT*,'u ='
            !DO j = N+2,1,-1
            !    DO i = 1, N+2
            !        WRITE(*,FMT='(1f8.3)',ADVANCE='no') u(i,j)
            !    END DO
            !    PRINT*,''
            !END DO

        END IF
    ELSE 
        tag = 0
        
        PRINT*,'u='
        DO j = N+2,1,-1
            DO i = 1, N+2
                u(i,j) = u_new(i,j+1)
                !WRITE(*,FMT='(1f8.3)',ADVANCE='no') u(i,j)
            END DO
            PRINT*,''
        END DO

    END IF


END SUBROUTINE prepare_u



SUBROUTINE init_1d(u0, f, ilow, ihi, jlow, jhi)
    INTEGER, INTENT(INOUT) :: ilow, ihi, jlow, jhi
    REAL(8), INTENT(INOUT) :: u0(ilow:ihi, jlow-nb:jhi+nb),f(ilow:ihi, jlow-nb:jhi+nb)
    INTEGER :: i, j
    REAL(8) :: x, y

    u0(:,:) = 0.1
    DO i = ilow, ihi
        DO j = jlow, jhi
            !u0(i,j) = j + 1
            f(i,j) = fn(i,j)
        END DO
    END DO


END SUBROUTINE init_1d



FUNCTION fn(i,j) RESULT(fx)
    INTEGER, INTENT(IN) :: i,j 
    REAL(8) :: fx, x, y

    x = i*h
    y = j*h

    IF((x-0.5)**2 + (y-0.5)**2 .LE. (0.3)**2) THEN
        fx = -500.0      
    ELSE
        fx = 0.0
    END IF
   
END FUNCTION fn



SUBROUTINE exchange_data_1(u_this, sbuffer, rbuffer, ilow, ihi, jlow, jhi)
    INTEGER, INTENT(INOUT) :: ilow, ihi, jlow, jhi
    REAL(8), INTENT(INOUT) :: u_this(ilow:ihi, jlow-nb:jhi+nb), sbuffer(ilow:ihi), rbuffer(ilow:ihi)
 
    INTEGER :: i,tag
    ! get neighbor ranks
    CALL MPI_CART_SHIFT(comm1d, 0, 1, bottom, top, ierr)
  
    !PRINT*,'Bottom, Top Neighbor Ranks = ',bottom,top
   
    ! even numbered ranks send first, odd numbered ranks recv first 
    IF(MOD(coord(1),2) .EQ. 0) THEN
        tag = 0
        ! send top row to neighbor above
        DO i = ilow, ihi
            sbuffer(i) = u_this(i,jhi) ! load buffer
        END DO

        CALL MPI_SEND(sbuffer, N+2, MPI_DOUBLE_PRECISION, top, tag, comm1d, ierr)

        ! receive bottom ghost row from neighbor below
        rbuffer = 0.0
        CALL MPI_RECV(rbuffer, N+2, MPI_DOUBLE_PRECISION, bottom, tag, comm1d, MPI_STATUS_IGNORE, ierr)
        DO i = ilow,ihi
            u_this(i,jlow-1) = rbuffer(i) ! copy recv buffer into array
        END DO 

        tag = 1
        ! send bottom row to neighbor below
        DO i = ilow, ihi
            sbuffer(i) = u_this(i,jlow) ! load buffer
        END DO
        CALL MPI_SEND(sbuffer, N+2, MPI_DOUBLE_PRECISION, bottom, tag, comm1d, ierr)

        ! receive top ghost row from neighbor above
        rbuffer = 0.0
        CALL MPI_RECV(rbuffer, N+2, MPI_DOUBLE_PRECISION, top, tag, comm1d, MPI_STATUS_IGNORE, ierr)
        DO i = ilow, ihi
            u_this(i,jhi+1) = rbuffer(i) ! copy recv buffer into array
        END DO 

    ELSE
        tag = 0
        ! receive bottom ghost row from neighbor below
        rbuffer = 0.0
        CALL MPI_RECV(rbuffer, N+2, MPI_DOUBLE_PRECISION, bottom, tag, comm1d, MPI_STATUS_IGNORE, ierr)
        DO i = ilow, ihi
            u_this(i,jlow-1) = rbuffer(i) ! copy recv buffer into array
        END DO 

        !PRINT*,'Rank#',myrank,' u_this='
        !DO j = jhi+1,jlow-1,-1
        !    DO i = ilow, ihi
        !        WRITE(*,FMT='(1f8.3)',ADVANCE='no') u_this(i,j)
        !    END DO
        !    PRINT*,''
        !END DO

        ! send top row to neighbor above
        DO i = ilow, ihi
            sbuffer(i) = u_this(i,jhi) ! load buffer
        END DO
        CALL MPI_SEND(sbuffer, N+2, MPI_DOUBLE_PRECISION, top, tag, comm1d, ierr)


        tag = 1
        ! receive top ghost row from neighbor above
        rbuffer = 0.0
        CALL MPI_RECV(rbuffer, N+2, MPI_DOUBLE_PRECISION, top, tag, comm1d, MPI_STATUS_IGNORE, ierr)
        DO i = ilow, ihi
            u_this(i,jhi+1) = rbuffer(i) ! copy recv buffer into array
        END DO

        ! send bottom row to neighbor below
        DO i = ilow, ihi
            sbuffer(i) = u_this(i,jlow) ! load buffer
        END DO
        CALL MPI_SEND(sbuffer, N+2, MPI_DOUBLE_PRECISION, bottom, tag, comm1d, ierr)

        

    END IF

END SUBROUTINE exchange_data_1



SUBROUTINE exchange_data_2(u_this, sbuffer, rbuffer, ilow, ihi, jlow, jhi)
    INTEGER, INTENT(INOUT) :: ilow, ihi, jlow, jhi
    REAL(8), INTENT(INOUT) :: u_this(ilow:ihi, jlow-nb:jhi+nb), sbuffer(ilow:ihi), rbuffer(ilow:ihi)
 
    INTEGER :: i,tag
    ! get neighbor ranks
    CALL MPI_CART_SHIFT(comm1d, 0, 1, bottom, top, ierr)
  
    !PRINT*,'Bottom, Top Neighbor Ranks = ',bottom,top
   
    tag = 0
    ! send top row to neighbor above and receive bottom ghost row from neighbor below
    DO i = ilow, ihi
        sbuffer(i) = u_this(i,jhi) ! load buffer
    END DO
    rbuffer = 0.0
    CALL MPI_SENDRECV(sbuffer,               &  ! send buffer
                      N+2,                     &  ! send count
                      MPI_DOUBLE_PRECISION,  &  ! data type
                      top,                   &  ! dest
                      tag,                   &  ! send tag
                      rbuffer,               &  ! recv buffer
                      N+2,                     &  ! recv count
                      MPI_DOUBLE_PRECISION,  &  ! data type
                      bottom,                &  ! source
                      tag,                   &  ! recv tag
                      comm1d, MPI_STATUS_IGNORE, ierr)
    DO i = ilow, ihi
        u_this(i,jlow-1) = rbuffer(i) ! copy recv buffer into array
    END DO 

    tag = 1
    ! send bottom row to neighbor below and receive top ghost row from neighbor above
    DO i = ilow, ihi
        sbuffer(i) = u_this(i,jlow) ! load buffer
    END DO
    rbuffer = 0.0
    CALL MPI_SENDRECV(sbuffer,               &  ! send buffer
                      N+2,                     &  ! send count
                      MPI_DOUBLE_PRECISION,  &  ! data type
                      bottom,                &  ! dest
                      tag,                   &  ! send tag
                      rbuffer,               &  ! recv buffer
                      N+2,                     &  ! recv count
                      MPI_DOUBLE_PRECISION,  &  ! data type
                      top,                   &  ! source
                      tag,                   &  ! recv tag
                      comm1d, MPI_STATUS_IGNORE, ierr)
    DO i = ilow, ihi
        u_this(i,jhi+1) = rbuffer(i) ! copy recv buffer into array
    END DO 
   

END SUBROUTINE exchange_data_2



SUBROUTINE exchange_data_3(u_this, sbuffertop, rbuffertop, sbufferbot, rbufferbot, &
                           ilow, ihi, jlow, jhi)
    INTEGER, INTENT(INOUT) :: ilow, ihi, jlow, jhi
    REAL(8), INTENT(INOUT) :: u_this(ilow:ihi, jlow-nb:jhi+nb)
    REAL(8), INTENT(INOUT) :: sbuffertop(ilow:ihi), rbuffertop(ilow:ihi), &
                              sbufferbot(ilow:ihi), rbufferbot(ilow:ihi)

 
    INTEGER :: i,tag
    ! get neighbor ranks
    CALL MPI_CART_SHIFT(comm1d, 0, 1, bottom, top, ierr)
  
    !PRINT*,'Bottom, Top Neighbor Ranks = ',bottom,top
   
    tag = 0
    ! send top row to neighbor above
    DO i = ilow, ihi
        sbuffertop(i) = u_this(i,jhi) ! load buffer
    END DO
    CALL MPI_ISEND(sbuffertop, N+2, MPI_DOUBLE_PRECISION, top, tag, comm1d, req(1), ierr)

    ! receive bottom ghost row from neighbor below
    CALL MPI_IRECV(rbufferbot, N+2, MPI_DOUBLE_PRECISION, bottom, tag, comm1d, req(2), ierr)
    
    tag = 1
    ! send bottom row to neighbor below
    DO i = ilow, ihi
        sbufferbot(i) = u_this(i,jlow) ! load buffer
    END DO
    CALL MPI_ISEND(sbufferbot, N+2, MPI_DOUBLE_PRECISION, bottom, tag, comm1d, req(3), ierr)

    ! receive top ghost row from neighbor above
    CALL MPI_IRECV(rbuffertop, N+2, MPI_DOUBLE_PRECISION, top, tag, comm1d, req(4), ierr)
   
   
    ! wait for all communications to complete
    CALL MPI_WAITALL(4, req, MPI_STATUSES_IGNORE, ierr)

    ! copy received data into array
    DO i = ilow, ihi
        u_this(i,jlow-1) = rbufferbot(i) 
    END DO 

    DO i = ilow, ihi
        u_this(i,jhi+1) = rbuffertop(i) 
    END DO 

END SUBROUTINE exchange_data_3



SUBROUTINE sweep(u_old, u_new, ilow, ihi, jlow, jhi)
    INTEGER, INTENT(INOUT) :: ilow, ihi, jlow, jhi
    REAL(8), INTENT(INOUT) :: u_old(ilow:ihi, jlow-nb:jhi+nb), &
                              u_new(ilow:ihi, jlow-nb:jhi+nb)
    INTEGER :: i, j

    DO i = ilow, ihi
        DO j = jlow, jhi
             u_new(i,j) = 0.25 * (u_old(i+1,j) + u_old(i-1,j) + u_old(i,j+1) + &
                          u_old(i,j-1) - h*h*f(i,j) )
        END DO
    END DO


END SUBROUTINE sweep



SUBROUTINE compute_difference(u_old, u_new, diff, ilow, ihi, jlow, jhi)
    INTEGER, INTENT(INOUT) :: ilow, ihi, jlow, jhi
    REAL(8), INTENT(INOUT) :: u_old(ilow:ihi, jlow-nb:jhi+nb), &
                              u_new(ilow:ihi, jlow-nb:jhi+nb), diff
    INTEGER :: i, j

    diff = 0.d0
    DO i = ilow, ihi
        DO j = jlow, jhi
             diff = diff + abs(u_old(i,j) - u_new(i,j)) 
        END DO
    END DO
    diff = diff * h * h

END SUBROUTINE compute_difference



SUBROUTINE flatten_array(u_this, u_flat)
    REAL(8), INTENT(IN) :: u_this(:,:)
    REAL(8), INTENT(INOUT) :: u_flat(:)
    INTEGER :: i, j, ix
  
    DO i = 1, N+2
        DO j = 2, 2-1+(N+2)/np
            ix = i+(j-2)*(N+2)
            u_flat(ix)  = u_this(i,j)   
        END DO
    END DO

END SUBROUTINE flatten_array



SUBROUTINE unflatten_save(u_flat, tag)
    REAL(8), INTENT(IN) :: u_flat(:)
    INTEGER, INTENT(IN) :: tag

    INTEGER :: jlow, jhi
    
    INTEGER :: i, j, ix

    jlow = 1+tag*(N+2)/np
    jhi = jlow-1+(N+2)/np

    DO i = 1, N+2
        DO j = jlow, jhi
            ix = i+(j-jlow)*(N+2)
            u(i,j) = u_flat(ix)  
        END DO
    END DO
    
END SUBROUTINE unflatten_save


END PROGRAM poisson_1
