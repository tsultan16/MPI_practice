! Overview:  This program numerically solves the 2D Poisson's equation (nabla(u(x,y)) = f(x,y)) using a 5-pt finite 
! difference approximation and Jacobi iterations on an (N+2) by (N+2) square grid: 
!
!      u_i,j^(k+1) = (1/4)*(u_i+1,j^(k) + u_i-1,j^(k) + u_i,j+1^(k) + u_i,j-1^(k) - h^2 * f_i,j ) 
!
! 2D Domain decomposition:  The 2d computational grid is decomposed into equally sized chunks.
! The decomposition is 2d. A process is assigned to each chunk. We use an MPI_CART communicator to
! to define the decomposition geometry and handle communication between processes.  
! Each sub-domain has a layer of boundary/ghost points which is 1-point deep.
! Ghost point data is exchanged between neighboring processes. 
!
! To ensure that processes can perform in parallel, we used combined blocking sendrecvs. 
!

PROGRAM poisson_2

USE MPI
IMPLICIT NONE


INTEGER :: comm2d, ierr, ndim, myrank, numprocs(1), dims(2), coord(2), req(4), &
           coltype, status(MPI_STATUS_SIZE), source, tag
LOGICAL :: isperiodic(2), reorder 
INTEGER, PARAMETER :: N = 126, nb = 1, maxiter = 1000
REAL(8), PARAMETER :: tol = 1.d-4
REAL(8) :: u(1:N+2,1:N+2), h, diff, difftot, ti, tf
REAL(8), ALLOCATABLE :: u_old(:,:),u_new(:,:), f(:,:), u_flat(:)
INTEGER :: ilow, ihi, jlow, jhi, i, j, iter, np, p(2)


! MPI initialization
CALL MPI_INIT(ierr)

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs(1), ierr)  

np = numprocs(1)

! create a new cartesian communicator for our 2d domain decomposition

isperiodic(1) = .FALSE. ! periodic boundaries
isperiodic(2) = .FALSE. ! periodic boundaries
reorder = .TRUE. ! allow MPI rank reordering 
ndim = 2 ! 2d decomposition
dims(1) = 2 !SQRT(REAL(np))    ! number of ranks in x direction
dims(2) = 4 !SQRT(REAL(np)) ! number of ranks in y direction


p(1) = dims(1)
p(2) = dims(2)

CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndim, dims, isperiodic, reorder, comm2d, ierr)

CALL MPI_COMM_RANK(comm2d, myrank, ierr)

! get rank coordinate
CALL MPI_CART_COORDS(comm2d, myrank, ndim, coord, ierr) 

!PRINT*,''
!PRINT*,'My rank, coordinate, num processes = ',myrank, coord, numprocs
!PRINT*,''

! compute boundary indices for this rank
CALL compute_bound(ilow, ihi, jlow, jhi)

h = 1.d0/(N+1)

! initialize u and f
ALLOCATE(u_old(ilow-nb:ihi+nb, jlow-nb:jhi+nb),u_new(ilow-nb:ihi+nb, jlow-nb:jhi+nb))
ALLOCATE(f(ilow-nb:ihi+nb, jlow-nb:jhi+nb))
ALLOCATE(u_flat(1:((N+2)/dims(1))*((N+2)/dims(2))))

CALL init_1d(u_old, f, ilow, ihi, jlow, jhi)
u_new = u_old
u = 0.0

OPEN(UNIT = 10, FILE = 'output.txt')
OPEN(UNIT = 11, FILE = 'horcut.txt')
OPEN(UNIT = 12, FILE = 'vercut.txt')


! create MPI derived data type for column buffer
CALL MPI_TYPE_VECTOR(1+jhi-jlow,            &  ! count 
                     1,                     &  ! block length
                     3+ihi-ilow,            &  ! stride
                     MPI_DOUBLE_PRECISION,  &  ! oldtype
                     coltype,               &  ! newtype
                     ierr)
CALL MPI_TYPE_COMMIT(coltype, ierr)
 

CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
ti = MPI_WTIME()

! Jacobi iteration loop
iter = 0
DO WHILE(iter .LT. maxiter)

    !PRINT*,'Rank # = ',myrank,' Iteration #',iter
    ! exchange ghost point data with neighbor ranks 
    !PRINT*,'Rank # = ',myrank,' Starting Exchange.'

    CALL exchange_data(u_old, ilow, ihi, jlow, jhi)

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
                       comm2d,                &  ! communicator handle
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


!PRINT*,'Rank# = ',myrank,' u_new = '
!DO j = jhi+1,jlow-1,-1
!    DO i = ilow-1, ihi+1
!        WRITE(*,FMT='(1f8.5)',ADVANCE='no') u_new(i,j)
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


DEALLOCATE(u_old, u_new, u_flat)

CLOSE(UNIT = 10)
CLOSE(UNIT = 11)
CLOSE(UNIT = 12)


! free MPI derived type
CALL MPI_TYPE_FREE(coltype, ierr)
! MPI termination
CALL MPI_FINALIZE(ierr)


CONTAINS


SUBROUTINE compute_bound(ilow, ihi, jlow, jhi)

    INTEGER, INTENT(INOUT) :: ilow, ihi, jlow, jhi

    ilow = 0+coord(1)*(N+2)/dims(1)
    ihi = ilow-1+(N+2)/dims(1)
    jlow = 0+coord(2)*(N+2)/dims(2)
    jhi = jlow-1+(N+2)/dims(2)

END SUBROUTINE compute_bound


SUBROUTINE init_1d(u0, f, ilow, ihi, jlow, jhi)
    INTEGER, INTENT(INOUT) :: ilow, ihi, jlow, jhi
    REAL(8), INTENT(INOUT) :: u0(ilow-nb:ihi+nb, jlow-nb:jhi+nb),f(ilow-nb:ihi+nb, jlow-nb:jhi+nb)
    INTEGER :: i, j
    REAL(8) :: x, y

    u0(:,:) = 0.1
    f(:,:) = 0.0
    DO i = ilow, ihi
        DO j = jlow, jhi
            f(i,j) = fn(i,j)
        END DO
    END DO


END SUBROUTINE init_1d



FUNCTION fn(i,j) RESULT(fx)
    INTEGER, INTENT(IN) :: i,j 
    REAL(8) :: fx, x, y

    x = i*h
    y = j*h

    IF((x-0.5)**2 + (y-0.5)**2 .LE. (0.2)**2) THEN
        fx = -500.0    
    ELSE
        fx = 0.0
    END IF
   
END FUNCTION fn



SUBROUTINE exchange_data(u_this, ilow, ihi, jlow, jhi)
    INTEGER, INTENT(INOUT) :: ilow, ihi, jlow, jhi
    REAL(8), INTENT(INOUT) :: u_this(ilow-nb:ihi+nb, jlow-nb:jhi+nb)
    INTEGER :: bottom, top, left, right

 
    INTEGER :: i,tag,nx
    ! get neighbor ranks
    CALL MPI_CART_SHIFT(comm2d, 0, 1, left, right, ierr) 
    CALL MPI_CART_SHIFT(comm2d, 1, 1, bottom, top, ierr)
      
 
    !PRINT*,'Bottom, Top Neighbor Ranks = ',bottom,top
   
    tag = 0
    ! send top row to neighbor above and receive bottom ghost row from neighbor below
    nx = 1+ihi-ilow
    CALL MPI_SENDRECV(u_this(ilow,jhi),      &  ! send buffer
                      nx,                    &  ! send count
                      MPI_DOUBLE_PRECISION,  &  ! data type
                      top,                   &  ! dest
                      tag,                   &  ! send tag
                      u_this(ilow,jlow-1),   &  ! recv buffer
                      nx,                    &  ! recv count
                      MPI_DOUBLE_PRECISION,  &  ! data type
                      bottom,                &  ! source
                      tag,                   &  ! recv tag
                      comm2d, MPI_STATUS_IGNORE, ierr)

    tag = 1
    ! send bottom row to neighbor below and receive top ghost row from neighbor above
    CALL MPI_SENDRECV(u_this(ilow,jlow),     &  ! send buffer
                      nx,                    &  ! send count
                      MPI_DOUBLE_PRECISION,  &  ! data type
                      bottom,                &  ! dest
                      tag,                   &  ! send tag
                      u_this(ilow,jhi+1),    &  ! recv buffer
                      nx,                    &  ! recv count
                      MPI_DOUBLE_PRECISION,  &  ! data type
                      top,                   &  ! source
                      tag,                   &  ! recv tag
                      comm2d, MPI_STATUS_IGNORE, ierr)
   

    tag = 0
    ! send rightmost column to right neighbor and receive left ghost column from left neighbor 
    CALL MPI_SENDRECV(u_this(ihi,jlow),      &  ! send buffer
                      1,                     &  ! send count
                      coltype,               &  ! data type
                      right,                 &  ! dest
                      tag,                   &  ! send tag
                      u_this(ilow-1,jlow),   &  ! recv buffer
                      1,                     &  ! recv count
                      coltype,               &  ! data type
                      left,                  &  ! source
                      tag,                   &  ! recv tag
                      comm2d, MPI_STATUS_IGNORE, ierr)
   
    tag = 1
    ! send leftmost column to left neighbor and receive right ghost column from right neighbor 
    CALL MPI_SENDRECV(u_this(ilow,jlow),     &  ! send buffer
                      1,                     &  ! send count
                      coltype,               &  ! data type
                      left,                  &  ! dest
                      tag,                   &  ! send tag
                      u_this(ihi+1,jlow),    &  ! recv buffer
                      1,                     &  ! recv count
                      coltype,               &  ! data type
                      right,                 &  ! source
                      tag,                   &  ! recv tag
                      comm2d, MPI_STATUS_IGNORE, ierr)


END SUBROUTINE exchange_data



SUBROUTINE sweep(u_old, u_new, ilow, ihi, jlow, jhi)
    INTEGER, INTENT(INOUT) :: ilow, ihi, jlow, jhi
    REAL(8), INTENT(INOUT) :: u_old(ilow-nb:ihi+nb, jlow-nb:jhi+nb), &
                              u_new(ilow-nb:ihi+nb, jlow-nb:jhi+nb)
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

SUBROUTINE prepare_u(u_new)

    REAL(8), INTENT(IN) :: u_new(:,:)
    INTEGER :: i, j, s_coord(2)

    IF(np .GT. 1) THEN 
        IF(myrank .NE. 0) THEN
            tag = myrank
            CALL flatten_array(u_new, u_flat)
            CALL MPI_SEND(u_flat, ((N+2)/dims(1))*((N+2)/dims(2)), MPI_DOUBLE_PRECISION, 0, tag, comm2d, ierr)
           
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
            CALL unflatten_save(u_flat, (/ 0, 0 /)) 
                
            DO i= 1, np-1
                tag = i
                source = i
                CALL MPI_RECV(u_flat, ((N+2)/dims(1))*((N+2)/dims(2)), MPI_DOUBLE_PRECISION, source, tag, comm2d, status, ierr)   

                ! get source coordinate
                CALL MPI_CART_COORDS(comm2d, tag, ndim, s_coord, ierr) 

             
                CALL unflatten_save(u_flat, s_coord)
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
        
        !PRINT*,'u='
        !DO j = N+2,1,-1
        !    DO i = 1, N+2
        !        u(i,j) = u_new(i,j+1)
        !        !WRITE(*,FMT='(1f8.3)',ADVANCE='no') u(i,j)
        !    END DO
        !    PRINT*,''
        !END DO

    END IF


END SUBROUTINE prepare_u



SUBROUTINE flatten_array(u_this, u_flat)
    REAL(8), INTENT(IN) :: u_this(:,:)
    REAL(8), INTENT(INOUT) :: u_flat(:)
    INTEGER :: i, j, ix
  
    DO i = 2, 2-1+(N+2)/p(1)
        DO j = 2, 2-1+(N+2)/p(2)
            ix = 1+(i-2)+(j-2)*((N+2)/p(1))
            u_flat(ix)  = u_this(i,j)   
        END DO
    END DO


END SUBROUTINE flatten_array



SUBROUTINE unflatten_save(u_flat, s_coord)
    REAL(8), INTENT(IN) :: u_flat(:)
    INTEGER, INTENT(IN) :: s_coord(2)

    INTEGER :: ilow, ihi, jlow, jhi
    
    INTEGER :: i, j, ix

    ilow = 1+s_coord(1)*(N+2)/p(1)
    ihi = ilow-1+(N+2)/p(1)
    jlow = 1+s_coord(2)*(N+2)/p(2)
    jhi = jlow-1+(N+2)/p(2)

    
    DO i = ilow, ihi
        DO j = jlow, jhi
            ix = 1+(i-ilow)+(j-jlow)*(N+2)/p(1)
            u(i,j) = u_flat(ix)  
        END DO
    END DO

    
END SUBROUTINE unflatten_save



END PROGRAM poisson_2
