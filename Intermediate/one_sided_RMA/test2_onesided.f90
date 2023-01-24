! This example program demonstrates halo exchange with 1d domain decomposition
! using one-sided RMA. Each process directly puts data from it's top and 
! bottom interior layers into the neighbor's RMA window. 'Fences' are set up around 
! the RMA put operation for active synchronization. After put's are completed, each process
! retrieves the data placed in it's RMA window by the neighbors.


PROGRAM test2_onesided

USE ISO_C_BINDING

IMPLICIT NONE

INCLUDE 'mpif.h'



INTEGER :: comm1d, ierr, ndim, myrank, numprocs(1), coord(1), req(4), status(MPI_STATUS_SIZE), source, tag
LOGICAL :: isperiodic(1), reorder 
INTEGER, PARAMETER :: n = 5
REAL(8), ALLOCATABLE :: A(:), u(:,:), buffer(:)
INTEGER :: i, j, np
INTEGER :: bottom, top
! RMA vars
INTEGER :: sizedouble, winsize
INTEGER :: win 


! MPI initialization
CALL MPI_INIT(ierr)

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs(1), ierr)  

np = numprocs(1)

IF(np .NE. 3) THEN
    PRINT*,'Need 3 ranks for this test...'
    STOP
END IF

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
    

! allocate RMA window and work arrays
ALLOCATE(A(2*n), u(n,n), buffer(n))
A = 0.d0
u = 0.d0
buffer = 0.d0


! determine array element size in bytes    
CALL MPI_SIZEOF(A(1), sizedouble, ierr)
winsize = (2*n)*sizedouble

IF(myrank .EQ. 0) PRINT*,'sizedouble, winsize = ',sizedouble, winsize 


! create rma local memory window
CALL MPI_Win_create(A, winsize, sizedouble, MPI_INFO_NULL, comm1d, win, ierr)


! put some test values inside u
IF(coord(1) .EQ. 1) THEN
    u(:,1) = (/ -1,1,2,3,-1 /) 
    u(:,n) = (/ -1,11,12,13,-1 /)     
END IF


PRINT*,'myrank, coord= ',myrank, coord(1)
PRINT*,''
PRINT*,''

CALL MPI_BARRIER(ierr)


IF(coord(1) .EQ. 1 ) THEN
PRINT*,''
PRINT*,'myrank =',myrank
DO  j = n, 1, -1
    DO i = 1, n
       WRITE(*,FMT='(f5.1)',ADVANCE='NO') u(i,j)
    END DO
    PRINT*,''
END DO    
PRINT*,''
END IF

CALL MPI_BARRIER(ierr)


! exchange these values with neighbors
CALL exchange1()

! unpack data from recv buffer
u(:,1) = A(1:n) 
u(:,n) = A(n+1:2*n) 

PRINT*,'Exchange Complete.'

IF(coord(1) .EQ. 2) THEN

PRINT*,''
PRINT*,'Myrank = ',myrank
PRINT*,''
PRINT*,'u='

DO  j = n, 1, -1
    DO i = 1, n
       WRITE(*,FMT='(f5.1)',ADVANCE='NO') u(i,j)
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
PRINT*,'u='

DO  j = n, 1, -1
    DO i = 1, n
       WRITE(*,FMT='(f5.1)',ADVANCE='NO') u(i,j)
    END DO
    PRINT*,''
END DO    


PRINT*,''
PRINT*,''

END IF


!CALL MPI_WIN_FREE(win,ierr)
CALL MPI_FINALIZE(ierr)

DEALLOCATE(u, A, buffer)

CONTAINS


! ghost cell exchange with neighboring ranks
SUBROUTINE exchange1()

    INTEGER :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND) :: bottom_ghost_disp, top_ghost_disp 
    
    ! region of RMA load/stores bounded by fences
    CALL MPI_WIN_FENCE(0, win, ierr)
            
        !PRINT*,'Neighbors: bottom, top =', bottom, top
        
        ! put bottom layer of interior into bottom neighbor's top ghost cells
        
        top_ghost_disp = n 
        buffer = u(:,1) 
        CALL MPI_PUT(buffer, n, MPI_DOUBLE_PRECISION, bottom, &
                     top_ghost_disp, n, MPI_DOUBLE_PRECISION, win, ierr)
    
        ! put top layer of interior into top neighbor's bottom ghost cells 
        
        bottom_ghost_disp = 0 ! displacement unit = 0 means put at the starting location of destination window
        buffer = u(:,n) 
        CALL MPI_PUT(buffer, n, MPI_DOUBLE_PRECISION, top, &
                     bottom_ghost_disp, n, MPI_DOUBLE_PRECISION, win, ierr)
    
    
    CALL MPI_WIN_FENCE(0, win, ierr)

    
END SUBROUTINE exchange1


END PROGRAM  test2_onesided