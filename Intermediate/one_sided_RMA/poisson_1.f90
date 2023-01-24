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
! We use MPI RMA, i.e. one-sided communication between ranks. Ranks will have a designated memory window 
! that will serve as the RMA buffer.
!
!
! Note: Using unproperly ordered blocking sends/recvs can result in deadlock or even cause the process
! to run sequentially.

PROGRAM poisson_1


IMPLICIT NONE
INCLUDE 'mpif.h'


INTEGER :: comm1d, ierr, ndim, myrank, numprocs(1), coord(1), req(4), status(MPI_STATUS_SIZE), source, tag
LOGICAL :: isperiodic(1), reorder 
INTEGER, PARAMETER :: N = 126, nb = 1, maxiter = 1000
REAL(8), PARAMETER :: tol = 1.d-4
REAL(8) :: h, diff, difftot, ti, tf
REAL(8), ALLOCATABLE :: uold(:,:),unew(:,:), f(:,:)
REAL(8), ALLOCATABLE :: A(:,:)                       
INTEGER :: xlow, xhi, ylow, yhi, ylow_off, yhi_off, i, j, iter, np
INTEGER :: bottom, top
! RMA vars
INTEGER :: sizedouble
INTEGER :: win 


! MPI initialization
CALL MPI_INIT(ierr)

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs(1), ierr)  

np = numprocs(1)

! create a new cartesian communicator for our 1d domain decomposition
isperiodic(1) = .FALSE. ! periodic boundaries
reorder = .TRUE. ! allow MPI rank reordering 
ndim = 1 ! 1d decomposition


CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndim, numprocs, isperiodic, reorder, comm1d, ierr)

CALL MPI_COMM_RANK(comm1d, myrank, ierr)

! get rank coordinate
CALL MPI_CART_COORDS(comm1d, myrank, ndim, coord, ierr) 


! set grid spacing
h = 1.d0/DBLE(N+1)

! set domains bounds
xlow = 0
xhi = N
ylow = 0 + coord(1) * (N/np + 1)
yhi = ylow + N/np 

ylow_off = 0
yhi_off = 0
IF(coord(1) .EQ. 0) ylow_off = 1
IF(coord(1) .EQ. np-1) yhi_off = 1

! allocate memory for work arrays and source term
ALLOCATE(uold(xlow:xhi,ylow:yhi), unew(xlow:xhi,ylow:yhi))
ALLOCATE(f(0:N,0:N))

uold = 0.d0
unew = 0.d0
f = 0.d0

! set source term
CALL setf()


! set initial values
uold = 0.d0
unew = 0.d0

! create RMA buffer and window object
ALLOCATE(A(xlow:xhi, ylow:yhi))

! determine array element size in bytes    
CALL MPI_SIZEOF(A(xlow,ylow), sizedouble, ierr)

! create MPI RMA window object
CALL MPI_WIN_CREATE(A,(1+xhi-xlow)*(1+yhi-ylow)*sizedouble, sizedouble, MPI_INFO_NULL, MPI_COMM_WORLD, win, ierr)



! start the iteration loop
DO iter = 1, maxiter

    ! exachange ghost cells
    CALL exchange1()

    ! update u
    DO j = ylow+ylow_off, yhi-yhi_off
        DO i = xlow+1, xhi-1
    
            unew(i,j) = 0.25d0 * (uold(i-1,j) + uold(i+1,j) + uold(i,j-1) + uold(i,j+1)) - &
                        0.25d0 * (h**2) * f(i,j) 
        
        END DO
    END DO

    ! save to file


END DO



CONTAINS


SUBROUTINE exchange1()

    INTEGER :: ierr


END SUBROUTINE exchange1




    





SUBROUTINE setf()

   INTEGER :: i,j
   
   
   DO j = 0, N
       DO i = 0, N
            IF(i .GE. N/2 - 5 .AND. i .LE. N/2 + 5 .AND. &
               j .GE. N/2 - 5 .AND. j .LE. N/2 + 5  ) THEN
            
                f = 1.d0
                
            ELSE
            
                f = 0.d0
                
            END IF
           
       END DO
   END DO    
   

END SUBROUTINE setf



END PROGRAM poisson_1