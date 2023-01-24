! In this example program, a 2d domain is decomposed into equally sized chunks
! The decomposition is 1d, i.e. the chunks lie along one direction.
! A process is assigned to each chunk. We use an MPI_CART communicator to
! to define the decomposition geometry and handle communication between processes.  
! Each sub-domain has a layer of boundary/ghost points which is 1-point deep.
! Ghost point data is exchanged between neighboring sub-domains.

PROGRAM domain_decomp_1d

USE MPI
IMPLICIT NONE


INTEGER :: comm1d, ierr, ndim, myrank, numprocs(1)
LOGICAL :: isperiodic(1), reorder 
INTEGER, PARAMETER :: N = 8
INTEGER :: u(0:N+1,0:N+1)
INTEGER, ALLOCATABLE :: u_this(:,:)
INTEGER :: ilow, ihi, jlow, jhi, i, j
INTEGER :: bottom, top


! MPI initialization
CALL MPI_INIT(ierr)

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs(1), ierr)  

! create a new cartesian communicator for our 2d domain decomposition

isperiodic(1) = .TRUE. ! periodic boundaries
reorder = .TRUE. ! allow MPI rank reordering 
ndim = 1 ! 1d decomposition


CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndim, numprocs, isperiodic, reorder, comm1d, ierr)

CALL MPI_COMM_RANK(comm1d, myrank, ierr)


PRINT*,''
PRINT*,'My rank = ',myrank
PRINT*,''

! compute boundary indices for this rank
CALL compute_bound()

! initialize matrix u
ALLOCATE(u_this(ilow:ihi, jlow-1:jhi+1))
u_this(:,:) = 0
DO j = jlow, jhi
    DO i = ilow, ihi
        u_this(i,j) = i+(j-1)*(N+2)
    END DO
END DO

PRINT*,'jlo,jhi=',jlow,jhi

! set ghost point values
!DO j = jlow, jhi
!    u_this(ilow-1,j) = u_this(ilow,j)
!    u_this(ihi+1,j) = u_this(ihi,j)
!END DO
!DO i = ilow, ihi
!    u_this(i,jlow-1) = u_this(i,jlow)
!    u_this(i,jhi+1) = u_this(i,jhi)
!END DO

PRINT*,'u = '
DO j = jhi+1,jlow-1,-1
    DO i = ilow, ihi
        WRITE(*,FMT =('(i4)'),ADVANCE='no') u_this(i,j)
    END DO
    PRINT*,''
END DO




! exchange ghost point data with neighbor ranks 
CALL exchange_data()

PRINT*,'Exchange complete.'

PRINT*,'u = '
DO j = jhi+1,jlow-1,-1
    DO i = ilow, ihi
        WRITE(*,FMT =('(i4)'),ADVANCE='no') u_this(i,j)
    END DO
    PRINT*,''
END DO




! MPI termination
CALL MPI_FINALIZE(ierr)


CONTAINS


SUBROUTINE compute_bound()

    ilow = 0
    ihi = N+1
    jlow = 0+myrank*(N+2)/numprocs(1)
    jhi = jlow-1+(N+2)/numprocs(1)

END SUBROUTINE compute_bound


SUBROUTINE exchange_data()
 
    INTEGER :: tag
    ! get neighbor ranks
    CALL MPI_CART_SHIFT(comm1d, 0, 1, bottom, top, ierr)
  
    PRINT*,'Bottom, Top Neighbor Ranks = ',bottom,top
   
    tag = 0
    ! send top row to neighbor above
    CALL MPI_SEND(u_this(ilow,jhi), N, MPI_DOUBLE_PRECISION, top, tag, comm1d, ierr)
    ! receive bottom ghost row from neighbor below
    CALL MPI_RECV(u_this(ilow,jlow-1), N, MPI_DOUBLE_PRECISION, bottom, tag, comm1d, MPI_STATUS_IGNORE, ierr)

    tag = 1
    ! send bottom row to neighbor below
    CALL MPI_SEND(u_this(ilow,jlow), N, MPI_DOUBLE_PRECISION, bottom, tag, comm1d, ierr) 
    ! receive top ghost row from neighbor above
    CALL MPI_RECV(u_this(ilow,jhi+1), N, MPI_DOUBLE_PRECISION, top, tag, comm1d, MPI_STATUS_IGNORE, ierr)


END SUBROUTINE exchange_data

END PROGRAM domain_decomp_1d

