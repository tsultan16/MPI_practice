! This example program demonstrates halo exchange with 2d domain decomposition
! using one-sided RMA. Each process directly puts data from it's top and 
! bottom interior layers into the neighbor's RMA window. 
! We use passive target RMA and acheive synchronization via MPI_WIN_FLUSH. 

PROGRAM test6_onesided

USE ISO_C_BINDING
!USE MPI

IMPLICIT NONE

INCLUDE 'mpif.h'


!*********************************************************************************************************************************
INTEGER, PARAMETER :: n = 5
INTEGER, PARAMETER :: nranks_x = 2  ! # of ranks in x-direction
INTEGER, PARAMETER :: nranks_y = 2  ! # of ranks in y-direction
INTEGER, PARAMETER :: ndim = 2     ! dimensionality of domain decomposition


INTEGER :: comm2d, ierr, myrank, dims(2), numprocs(1), coord(2), req(4), status(MPI_STATUS_SIZE), source, tag
LOGICAL :: isperiodic(2), reorder 
REAL(8), ALLOCATABLE :: u(:,:), buffer(:)
INTEGER :: i, j, k, np
INTEGER :: bottom, top, left, right, neighbor
! RMA vars
INTEGER :: sizedouble, winsize
INTEGER(KIND=MPI_ADDRESS_KIND) :: disp
INTEGER :: win 
TYPE(C_PTR) :: rma_cmem  ! C pointer to start of RMA window memory block
REAL(8), POINTER :: A(:) ! Fortran pointer that will be associated with the RMA window's C pointer


!*********************************************************************************************************************************



! MPI initialization

CALL rma_init()

! open RMA exposure epoch
CALL open_exposure_epoch()

CALL MPI_BARRIER(ierr)


! put some test values inside u
IF(myrank .EQ. 0) THEN
    u = -1
    u(2:n-1,2) = (/ 1,1,1 /)     
    u(2:n-1,3) = (/ 2,2,2 /)     
    u(2:n-1,n-1) = (/ 3,3,3 /)     
    

    PRINT*,''
    PRINT*,'Myrank = ',myrank
    PRINT*,''
    PRINT*,'u ='
    DO j = n, 1, -1
    DO i = 1, n
        WRITE(*,FMT='(f5.1)', ADVANCE='NO') u(i,j)
    END DO    
    PRINT*,''
    END DO
    
    CALL sleep(1)
    
END IF

CALL MPI_BARRIER(ierr)

IF(myrank .EQ. 1) THEN
    
    u = -2
    u(2:n-1,2) = (/ 4,4,4 /)     
    u(2:n-1,3) = (/ 5,5,5 /)     
    u(2:n-1,n-1) = (/ 6,6,6 /)  
    

    PRINT*,''
    PRINT*,'Myrank = ',myrank
    PRINT*,''
    PRINT*,'u ='
    DO j = n, 1, -1
    DO i = 1, n
        WRITE(*,FMT='(f5.1)', ADVANCE='NO') u(i,j)
    END DO
    PRINT*,''
    END DO

    CALL sleep(1)
    
END IF

CALL MPI_BARRIER(ierr)

IF(myrank .EQ. 2) THEN
    
    u = -3
    u(2:n-1,2) = (/ 7,7,7 /)     
    u(2:n-1,3) = (/ 8,8,8 /)     
    u(2:n-1,n-1) = (/ 9,9,9 /)  
    

    PRINT*,''
    PRINT*,'Myrank = ',myrank
    PRINT*,''
    PRINT*,'u ='
    DO j = n, 1, -1
    DO i = 1, n
        WRITE(*,FMT='(f5.1)', ADVANCE='NO') u(i,j)
    END DO
    PRINT*,''
    END DO
    
    CALL sleep(1)
    
END IF

CALL MPI_BARRIER(ierr)


IF(myrank .EQ. 3) THEN
    
    u = -4
    u(2:n-1,2) = (/ 10,10,10 /)     
    u(2:n-1,3) = (/ 11,11,11 /)     
    u(2:n-1,n-1) = (/ 12,12,12 /)   
    

    PRINT*,''
    PRINT*,'Myrank = ',myrank
    PRINT*,''
    PRINT*,'u ='
    DO j = n, 1, -1
    DO i = 1, n
        WRITE(*,FMT='(f5.1)', ADVANCE='NO') u(i,j)
    END DO
    PRINT*,''
    END DO

    CALL sleep(1)
    
END IF

CALL MPI_BARRIER(ierr)


!GO TO 111


! tell neighbors that we are ready to receive
CALL recv_signal()

! poll on signal from from neighbors telling us they are ready to receive
CALL poll_signal(A,1)

! send data to neighbors since they are ready
CALL send_msg()

! tell neighbors that we hare done sending
CALL data_sent_signal()

! poll on signal from from neighbors telling us they have sent us data_sent_signal
CALL poll_signal(A,2)


PRINT*,''
PRINT*,'One sided communication complete...Myrank = ',myrank
PRINT*,''


CALL MPI_BARRIER(ierr)


!unpack data from recv buffer

u(1,:) = A(3:3+n-1)  
u(n,:) = A(3+n+2:3+n+2+n-1)
u(:,1) = A(3+2*(n+2):3+2*(n+2)+n-1)
u(:,n) = A(3+3*(n+2):3+3*(n+2)+n-1)


DO k = 0, np-1

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

IF(myrank .EQ. k) THEN

    PRINT*,''
    PRINT*,'Myrank = ',myrank
    PRINT*,''
    PRINT*,'u ='
    DO j = n, 1, -1
    DO i = 1, n
        WRITE(*,FMT='(f5.1)', ADVANCE='NO') u(i,j)
    END DO    
    PRINT*,''
    END DO
    CALL SLEEP(1)
 
END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

END DO



!111 CONTINUE



CALL close_rma_window()

PRINT*,'Done..myrank = ', myrank

DEALLOCATE(u, buffer)




CONTAINS



! send signal to neighbor indicating that we are ready to receive data
SUBROUTINE recv_signal()

 INTEGER :: ierr

    buffer(1) = 1.d0
    
    ! Put in right neighbors window
    disp = 0    
    CALL MPI_PUT(buffer, 1, MPI_DOUBLE_PRECISION, right, &
                 disp, 1, MPI_DOUBLE_PRECISION, win, ierr)
    
    ! Put in left neighbors window
    disp = 0 + n+2    
    CALL MPI_PUT(buffer, 1, MPI_DOUBLE_PRECISION, left, &
                 disp, 1, MPI_DOUBLE_PRECISION, win, ierr)
           

    ! Put in top neighbors window
    disp = 0 + 2*(n+2)    
    CALL MPI_PUT(buffer, 1, MPI_DOUBLE_PRECISION, top, &
                 disp, 1, MPI_DOUBLE_PRECISION, win, ierr)
    
    ! Put in bottom neighbors window
    disp = 0 + 3*(n+2)    
    CALL MPI_PUT(buffer, 1, MPI_DOUBLE_PRECISION, bottom, &
                 disp, 1, MPI_DOUBLE_PRECISION, win, ierr)
    
    CALL MPI_WIN_FLUSH(neighbor, win, ierr)


END SUBROUTINE recv_signal


! sends message to neighbor
SUBROUTINE send_msg()

    INTEGER :: ierr
    

    ! Put in right neighbors window
    disp = 2    
    buffer(:) =  u(n-1,:)
    CALL MPI_PUT(buffer, n, MPI_DOUBLE_PRECISION, right, &
                 disp, n, MPI_DOUBLE_PRECISION, win, ierr)
    
    ! Put in left neighbors window
    disp = 2 + n+2    
    buffer(:) =  u(2,:)
    CALL MPI_PUT(buffer, n, MPI_DOUBLE_PRECISION, left, &
                 disp, n, MPI_DOUBLE_PRECISION, win, ierr)
           

    ! Put in top neighbors window
    disp = 2 + 2*(n+2)    
    buffer(:) =  u(:,n-1)
    CALL MPI_PUT(buffer, n, MPI_DOUBLE_PRECISION, top, &
                 disp, n, MPI_DOUBLE_PRECISION, win, ierr)
    
    ! Put in bottom neighbors window
    disp = 2 + 3*(n+2)    
    buffer(:) =  u(:,2)
    CALL MPI_PUT(buffer, n, MPI_DOUBLE_PRECISION, bottom, &
                 disp, n, MPI_DOUBLE_PRECISION, win, ierr)
                 
    CALL MPI_WIN_FLUSH(neighbor, win, ierr)


END SUBROUTINE send_msg


! send a signal to neighbor indicating that we have sent them data
SUBROUTINE data_sent_signal()

    INTEGER :: ierr

    buffer(1) = 1.d0

    ! First half of window is reserved for data from bottom neighbor
    ! The second half is for top neighbor.
    
     ! Put in right neighbors window
    disp = 1   
    CALL MPI_PUT(buffer, 1, MPI_DOUBLE_PRECISION, right, &
                 disp, 1, MPI_DOUBLE_PRECISION, win, ierr)
    
    ! Put in left neighbors window
    disp = 1 + n+2    
    CALL MPI_PUT(buffer, 1, MPI_DOUBLE_PRECISION, left, &
                 disp, 1, MPI_DOUBLE_PRECISION, win, ierr)
           

    ! Put in top neighbors window
    disp = 1 + 2*(n+2)    
    CALL MPI_PUT(buffer, 1, MPI_DOUBLE_PRECISION, top, &
                 disp, 1, MPI_DOUBLE_PRECISION, win, ierr)
    
    ! Put in bottom neighbors window
    disp = 1 + 3*(n+2)    
    CALL MPI_PUT(buffer, 1, MPI_DOUBLE_PRECISION, bottom, &
                 disp, 1, MPI_DOUBLE_PRECISION, win, ierr)
                 
                 
    CALL MPI_WIN_FLUSH(neighbor, win, ierr)


END SUBROUTINE data_sent_signal


! polls on signal from neighbors
SUBROUTINE poll_signal(rma_window, signal_offset)

    REAL(8), VOLATILE :: rma_window(:) ! must include the volatile attribute here.. to ensure that compiler isn't using a copy of the RMA window that it has pre-loaded into a register
    INTEGER, INTENT(IN) :: signal_offset  ! 1 for recv signal, 2 for sent signal
    INTEGER :: ierr

    
    
    ! poll the signal status
    DO WHILE(rma_window(signal_offset) .NE. 1.d0 .AND. rma_window(signal_offset + n+2) .NE. 1.d0 .AND. &
             rma_window(signal_offset +2*(n+2)) .NE. 1.d0 .AND. rma_window(signal_offset + 3*(n+2)) .NE. 1.d0 )
    
        ! wait till signals from all neighbors have arrived...
    
    END DO

    ! reset the signal status
    rma_window(signal_offset) = 0.d0    
    rma_window(signal_offset + n+2) = 0.d0
    rma_window(signal_offset + 2*(n+2)) = 0.d0
    rma_window(signal_offset + 3*(n+2)) = 0.d0


END SUBROUTINE poll_signal



! start MPI passive RMA epoch
SUBROUTINE open_exposure_epoch()

    INTEGER :: ierr


    CALL MPI_WIN_LOCK(MPI_LOCK_SHARED, neighbor, MPI_MODE_NOCHECK, win, ierr)
    

END SUBROUTINE open_exposure_epoch


! shut down RMA passive epoch and close the window
SUBROUTINE close_rma_window()

    INTEGER ::  ierr

    
    CALL MPI_WIN_UNLOCK(neighbor, win, ierr)
    
    CALL MPI_WIN_FREE(win,ierr)
    
    CALL MPI_FINALIZE(ierr)
    
    NULLIFY(A)

END SUBROUTINE close_rma_window



SUBROUTINE rma_init()

    INTEGER :: i

    CALL MPI_INIT(ierr)

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs(1), ierr)  

    np = numprocs(1)

    !IF(np .NE. 3) THEN
    !    PRINT*,'Need 3 ranks for this test...'
    !    STOP
    !END IF

    ! create a new cartesian communicator for our 2d domain decomposition
    dims(1) = nranks_x 
    dims(2) = nranks_y
    isperiodic(1) = .TRUE. ! periodic boundaries
    isperiodic(2) = .TRUE. ! periodic boundaries
    reorder = .TRUE. ! allow MPI rank reordering 


    CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndim, dims, isperiodic, reorder, comm2d, ierr)

    CALL MPI_COMM_RANK(comm2d, myrank, ierr)

    IF(myrank .EQ. 0) PRINT*,'# processes = ',np


    ! get rank coordinate
    CALL MPI_CART_COORDS(comm2d, myrank, ndim, coord, ierr) 

    ! get neighbor ranks
    CALL MPI_CART_SHIFT(comm2d, 0, 1, bottom, top, ierr)
    CALL MPI_CART_SHIFT(comm2d, 1, 1, left, right, ierr)
    
    DO i = 0, np-1
     
        CALL MPI_BARRIER(comm2d, ierr)
        IF(myrank .EQ. i) THEN
    
        PRINT*,''
        PRINT*,'###################################################################'
        WRITE(*,FMT='("  Myrank = ",i2,", Left = ",i2,", Right = ",i2,", Bottom = ",i2,", Top = ",i2)') myrank, left, right, bottom, top
        PRINT*,'###################################################################'
        PRINT*,''
        
        END IF
        CALL MPI_BARRIER(comm2d, ierr)
    
    END DO
    
    ! allocate work arrays
    ALLOCATE(u(n,n), buffer(n))
    u = 0.d0
    buffer = 0.d0

    ! determine array element size in bytes    
    CALL MPI_SIZEOF(DBLE(1), sizedouble, ierr)
    winsize = (4*(n+2))*sizedouble

    IF(myrank .EQ. 0) PRINT*,'sizedouble, winsize = ',sizedouble, winsize 


    ! create rma local memory allocation and create window
    CALL MPI_WIN_ALLOCATE(winsize, sizedouble, MPI_INFO_NULL, comm2d, rma_cmem, win, ierr)
    CALL C_F_POINTER(rma_cmem, A, (/ 4*(n+2) /)) ! associate C pointer returned by win_allocate to our fortran pointer

    ! Note: The RMA window is broken up into 4 chucks corresponding to ghost cells: x-, x+, y-, y+ (in that order)
    ! The first two slots in each chunk are reserved for the data recv/sent signals from the neighbor corresponding 
    ! to that boundary (hence the size n+2 for each chuck).

    A(:) = 0.d0
    
    PRINT*,'MPI initialization done..'

END SUBROUTINE rma_init



END PROGRAM  test6_onesided
