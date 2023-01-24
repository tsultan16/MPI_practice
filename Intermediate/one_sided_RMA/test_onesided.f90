
PROGRAM test_onesided

IMPLICIT NONE

INCLUDE 'mpif.h'



INTEGER :: comm1d, ierr, ndim, myrank, numprocs(1), coord(1), req(4), status(MPI_STATUS_SIZE), source, tag
LOGICAL :: isperiodic(1), reorder 
INTEGER, PARAMETER :: n = 3, nb = 1, maxiter = 1000
REAL(8), PARAMETER :: tol = 1.d-4
REAL(8) :: h, diff, difftot, ti, tf
REAL(8), ALLOCATABLE :: A(:,:), buffer_top(:), buffer_bottom(:)
INTEGER :: i, j, iter, np
INTEGER :: bottom, top
! RMA vars
INTEGER :: sizedouble, winsize
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

IF(myrank .EQ. 0) PRINT*,'# processes = ',np


! get rank coordinate
CALL MPI_CART_COORDS(comm1d, myrank, ndim, coord, ierr) 

! get neighbor ranks
CALL MPI_CART_SHIFT(comm1d, 0, 1, bottom, top, ierr)
    
! set grid spacing
h = 1.d0/DBLE(N+1)


! RMA window
ALLOCATE(A(1-nb:n+nb,1-nb:n+nb), buffer_top(n+2*nb), buffer_bottom(n+2*nb))
A = 0.d0
buffer_top = 0.d0
buffer_bottom = 0.d0

! determine array element size in bytes    
CALL MPI_SIZEOF(A(1,1), sizedouble, ierr)
winsize = (n+2*nb)*(n+2*nb)*sizedouble

IF(myrank .EQ. 0) PRINT*,'sizedouble, winsize = ',sizedouble, winsize 


! create rma local memory window
CALL MPI_Win_create(A, winsize, sizedouble, MPI_INFO_NULL, comm1d, win, ierr)


! put some test values inside A
IF(coord(1) .EQ. 1) THEN
    A(:,1) = (/ -1,1,2,3,-1 /) 
    A(:,n) = (/ -1,11,12,13,-1 /) 
    
    buffer_top = A(:,n)
    buffer_bottom = A(:,1)
END IF


PRINT*,'myrank, coord= ',myrank, coord(1)
PRINT*,''
PRINT*,''

CALL MPI_BARRIER(ierr)


IF(coord(1) .EQ. 1 ) THEN
PRINT*,''
PRINT*,'myrank =',myrank
DO  j = n+1, 0, -1
    DO i = 0, n+1
       WRITE(*,FMT='(f5.1)',ADVANCE='NO') A(i,j)
    END DO
    PRINT*,''
END DO    
PRINT*,''
END IF

CALL MPI_BARRIER(ierr)




! exchange these values with neighbors
CALL exchange1()



PRINT*,'Exchange Complete.'

IF(coord(1) .EQ. 2) THEN

PRINT*,''
PRINT*,'Myrank = ',myrank
PRINT*,''
PRINT*,'A='

DO  j = n+1, 0, -1
    DO i = 0, n+1
       WRITE(*,FMT='(f5.1)',ADVANCE='NO') A(i,j)
    END DO
    PRINT*,''
END DO    


PRINT*,''
PRINT*,''

END IF

CALL MPI_BARRIER(ierr)


IF(coord(1) .EQ. 0) THEN

PRINT*,''
PRINT*,'Myrank = ',myrank
PRINT*,''
PRINT*,'A='

DO  j = n+1, 0, -1
    DO i = 0, n+1
       WRITE(*,FMT='(f5.1)',ADVANCE='NO') A(i,j)
    END DO
    PRINT*,''
END DO    


PRINT*,''
PRINT*,''

END IF


!CALL MPI_WIN_FREE(win,ierr)
CALL MPI_FINALIZE(ierr)

DEALLOCATE(A, buffer_top, buffer_bottom)

CONTAINS


! ghost cell exchange with neighboring ranks
SUBROUTINE exchange1()

    INTEGER :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND) :: bottom_ghost_disp, top_ghost_disp 
    
    ! region of RMA load/stores bounded by fences
    CALL MPI_WIN_FENCE(0, win, ierr)
            
        !PRINT*,'Neighbors: bottom, top =', bottom, top
        
        ! put bottom layer of interior into bottom neighbor's top ghost cells 
        top_ghost_disp = 1-nb + (n+2*nb)*(n+2*nb - 1)
        CALL MPI_PUT(buffer_bottom, n+2*nb, MPI_DOUBLE_PRECISION, bottom, &
                     top_ghost_disp, n+2*nb, MPI_DOUBLE_PRECISION, win, ierr)
    
        ! put top layer of interior into top neighbor's bottom ghost cells 
        bottom_ghost_disp = 1-nb 
        CALL MPI_PUT(buffer_top, n+2*nb, MPI_DOUBLE_PRECISION, top, &
                     bottom_ghost_disp, n+2*nb, MPI_DOUBLE_PRECISION, win, ierr)
    
    
    CALL MPI_WIN_FENCE(0, win, ierr)

    
END SUBROUTINE exchange1


END PROGRAM  test_onesided