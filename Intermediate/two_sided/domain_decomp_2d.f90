! In this example program, a 2d domain is decomposed into equally sized chunks.
! A process is assigned to each chunk. We use an MPI_CART communicator to
! to define the decomposition geometry and handle communication between processes.  
! Each sub-domain has a layer of boundary/ghost points which is 1-point deep.
! Ghost point data is exchanged between neighboring sub-domains.

PROGRAM domain_decomp_2d

USE MPI
IMPLICIT NONE


INTEGER :: comm2d, ierr, dims(2), ndim, coords(2), myrank
LOGICAL :: isperiodic(2), reorder 
INTEGER, PARAMETER :: N = 2
INTEGER :: u(0:N+1,0:N+1)
INTEGER, ALLOCATABLE :: u_this(:,:)
INTEGER :: ilow, ihi, jlow, jhi, i, j
INTEGER :: top, bottom, left, right


! MPI initialization
CALL MPI_INIT(ierr)

! create a new cartesian communicator for our 2d domain decomposition

dims(1) = 2 ! number of ranks in x direction
dims(2) = 2 ! number of ranks in y direction
isperiodic(:) = .TRUE. ! periodic boundaries
reorder = .TRUE. ! allow MPI rank reordering 
ndim = 2 ! 2d grid

CALL MPI_CART_CREATE(MPI_COMM_WORLD,   &  ! default communicator handle 
                     ndim,             &  ! cartesian grid dimensions
                     dims,             &  ! # of ranks in each direction 
                     isperiodic,       &  ! periodic boundary conditions setting
                     reorder,          &  ! rank reordering setting
                     comm2d,           &  ! new cart communicator handle
                     ierr)


! get MPI cartesian grid co-ordinates of this process
CALL MPI_COMM_RANK(comm2d, myrank, ierr)
CALL MPI_CART_COORDS(comm2d, myrank, ndim, coords, ierr)

PRINT*,''
PRINT*,'My rank = ',myrank
PRINT*,'Coordinates =(',coords(1),',',coords(2),')'
PRINT*,''

! compute boundary indices for this rank
CALL compute_bound()

! initialize matrix u
ALLOCATE(u_this(ilow-1:ihi+1, jlow-1:jhi+1))
u_this(:,:) = 0
DO j = jlow, jhi
    DO i = ilow, ihi
        u_this(i,j) = i+(j-1)*(N+2)
    END DO
END DO

!PRINT*,'ilo,ihi,jlo,jhi=',ilow,ihi,jlow,jhi

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
!DO j = jhi+1,jlow-1,-1
!    DO i = ilow-1, ihi+1
!        WRITE(*,FMT =('(i4)'),ADVANCE='no') u_this(i,j)
!    END DO
!    PRINT*,''
!END DO




! exchange ghost point data with neighbor ranks 
CALL exchange_data()




! MPI termination
CALL MPI_FINALIZE(ierr)


CONTAINS


SUBROUTINE compute_bound()

    ilow = 0+coords(1)*(N+2)/dims(1)
    ihi = ilow-1+(N+2)/dims(1)
    jlow = 0+coords(2)*(N+2)/dims(2)
    jhi = jlow-1+(N+2)/dims(2)

END SUBROUTINE compute_bound


SUBROUTINE exchange_data()
 
    ! get neighbor ranks
    CALL MPI_CART_SHIFT(comm2d, 0, 1, left, right, ierr)
    CALL MPI_CART_SHIFT(comm2d, 1, 1, bottom, top, ierr)

    !PRINT*,'Left, Right Neighbor Ranks = ',left,right
    !PRINT*,'Bottom, Top Neighbor Ranks = ',bottom,top

    ! sent top row to neighbor above
    CALL MPI_SEND(u_this(1,jhi), N, MPI_DOUBLE_PRECISION, top, 0, comm2d, ierr)
    ! receive from neighbor below
    CALL MPI_SEND(u_this(1,jlow-1), N, MPI_DOUBLE_PRECISION, bottom, 0, comm2d, MPI_STATUS_IGNORE, ierr)
    
    ! sent right column to right neighbor
    CALL MPI_SEND(u_this(1,jhi), N, MPI_DOUBLE_PRECISION, top, 0, comm2d, ierr)
    ! receive from left neighbor
    CALL MPI_SEND(u_this(1,jlow-1), N, MPI_DOUBLE_PRECISION, bottom, 0, comm2d, MPI_STATUS_IGNORE, ierr)
    


END SUBROUTINE exchange_data

END PROGRAM domain_decomp_2d

