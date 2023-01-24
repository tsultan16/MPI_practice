! This example program demonstrates halo exchange with 1d domain decomposition
! using one-sided RMA. Each process directly puts data from it's top and 
! bottom interior layers into the neighbor's RMA window. 
! We use passive target RMA and acheive synchronization via MPI_WIN_FLUSH. 

PROGRAM test5_onesided

USE ISO_C_BINDING

IMPLICIT NONE

INCLUDE 'mpif.h'


!*********************************************************************************************************************************

INTEGER :: comm1d, ierr, ndim, dims(2), myrank, numprocs(1), coord(1), req(4), status(MPI_STATUS_SIZE), source, tag
LOGICAL :: isperiodic(1), reorder 
INTEGER, PARAMETER :: n = 5
REAL(8), ALLOCATABLE :: u(:,:), buffer(:)
INTEGER :: i, j, np
INTEGER :: bottom, top, neighbor
! RMA vars
INTEGER :: sizedouble, winsize
INTEGER(KIND=MPI_ADDRESS_KIND) :: disp
INTEGER :: win 
TYPE(C_PTR) :: rma_cmem  ! C pointer to start of RMA window memory block
REAL(8), POINTER :: A(:) ! Fortran pointer that will be associated with the RMA window's C pointer


!*********************************************************************************************************************************



! MPI initialization
CALL MPI_INIT(ierr)

CALL rma_init()

CALL MPI_BARRIER(ierr)

! put some test values inside u
IF(coord(1) .EQ. 0) THEN
    u(:,2) = (/ 1,2,3,4,5 /)     
    u(:,n-1) = (/ 6,7,8,9,10 /)     
    

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
    
END IF

CALL MPI_BARRIER(ierr)

IF(coord(1) .EQ. 1) THEN
    
    u(:,2) = (/ 11,12,13,14,15 /)     
    u(:,n-1) = (/ 16,17,18,19,20 /)     
    

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
    
END IF

CALL MPI_BARRIER(ierr)

IF(coord(1) .EQ. 2) THEN
    
    u(:,2) = (/ 21,22,23,24,25 /)     
    u(:,n-1) = (/ 26,27,28,29,30 /)     
    

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
    
END IF

CALL MPI_BARRIER(ierr)


! open RMA exposure epoch
CALL open_exposure_epoch()



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
PRINT*,'One sided sommunication complete...Myrank = ',myrank
PRINT*,''


CALL MPI_BARRIER(ierr)


!unpack data from recv buffer

u(:,1) = A(3:3+n-1)
u(:,n) = A(3+n+2:3+n+2+n-1)


CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)


IF(myrank .EQ. 2) THEN

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


IF(myrank .EQ. 1) THEN

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

IF(myrank .EQ. 0) THEN

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


CALL close_rma_window()

PRINT*,'Done..myrank = ', myrank

DEALLOCATE(u, buffer)




CONTAINS



! send signal to neighbor indicating that we are ready to receive data
SUBROUTINE recv_signal()

 INTEGER :: ierr

    buffer(1) = 1.d0
    
    
    ! First half of window is reserved for data from bottom neighbor
    ! The second half is for top neighbor.
    
    ! First put signal into bottom neighbors window (second portion)
    disp = 0
    CALL MPI_PUT(buffer, 1, MPI_DOUBLE_PRECISION, top, &
                 disp, 1, MPI_DOUBLE_PRECISION, win, ierr)    
    
    
    ! Now put signal into top neighbors window (first portion)
    disp = 0 + (n+2)
    CALL MPI_PUT(buffer, 1, MPI_DOUBLE_PRECISION, bottom, &
                 disp, 1, MPI_DOUBLE_PRECISION, win, ierr)    
    
    
    CALL MPI_WIN_FLUSH(neighbor, win, ierr)


END SUBROUTINE recv_signal


! sends message to neighbor
SUBROUTINE send_msg()

    INTEGER :: ierr
    

    ! Put in top neighbors window
    disp = 0 + 2    
    buffer(:) =  u(:,n-1)
    CALL MPI_PUT(buffer, n, MPI_DOUBLE_PRECISION, top, &
                 disp, n, MPI_DOUBLE_PRECISION, win, ierr)
    
     ! Put in bottom neighbors window
    disp = 0 + 2 + n+2    
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
    
    ! First put signal into bottom neighbors window (second portion)
    disp = 0 + 1
    CALL MPI_PUT(buffer, 1, MPI_DOUBLE_PRECISION, top, &
                 disp, 1, MPI_DOUBLE_PRECISION, win, ierr)    
    
    
    ! Now put signal into top neighbors window (first portion)
    disp = 0 + (n+2) + 1
    CALL MPI_PUT(buffer, 1, MPI_DOUBLE_PRECISION, bottom, &
                 disp, 1, MPI_DOUBLE_PRECISION, win, ierr)    
                 
                 
    CALL MPI_WIN_FLUSH(neighbor, win, ierr)


END SUBROUTINE data_sent_signal


! polls on signal from neighbor
SUBROUTINE poll_signal(rma_window, signal_offset)

    REAL(8), VOLATILE :: rma_window(:) ! must include the volatile attribute here.. to ensure that compiler isn't using a copy of the RMA window that it has pre-loaded into a register
    INTEGER, INTENT(IN) :: signal_offset  ! 1 for recv signal, 2 for sent signal
    INTEGER :: ierr

    
    
    ! poll the signal status
    DO WHILE(rma_window(signal_offset) .NE. 1.d0 .AND. rma_window(signal_offset + n+2) .NE. 1.d0)
    
        ! wait till signal has arrived...
    
    END DO

    ! reset the signal status
    rma_window(signal_offset) = 0.d0    
    rma_window(signal_offset + n+2) = 0.d0


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

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs(1), ierr)  

    np = numprocs(1)

    IF(np .NE. 3) THEN
        PRINT*,'Need 3 ranks for this test...'
        STOP
    END IF

    ! create a new cartesian communicator for our 1d domain decomposition
    isperiodic(1) = .TRUE. ! periodic boundaries
    reorder = .TRUE. ! allow MPI rank reordering 
    ndim = 1 ! 1d decomposition


    CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndim, numprocs, isperiodic, reorder, comm1d, ierr)

    CALL MPI_COMM_RANK(comm1d, myrank, ierr)

    IF(myrank .EQ. 0) PRINT*,'# processes = ',np


    ! get rank coordinate
    CALL MPI_CART_COORDS(comm1d, myrank, ndim, coord, ierr) 

    ! get neighbor ranks
    CALL MPI_CART_SHIFT(comm1d, 0, 1, bottom, top, ierr)

    
    ! allocate work arrays
    ALLOCATE(u(n,n), buffer(n))
    u = 0.d0
    buffer = 0.d0

    ! determine array element size in bytes    
    CALL MPI_SIZEOF(DBLE(1), sizedouble, ierr)
    winsize = (2*(n+2))*sizedouble

    IF(myrank .EQ. 0) PRINT*,'sizedouble, winsize = ',sizedouble, winsize 


    ! create rma local memory allocation and create window
    CALL MPI_WIN_ALLOCATE(winsize, sizedouble, MPI_INFO_NULL, comm1d, rma_cmem, win, ierr)
    CALL C_F_POINTER(rma_cmem, A, (/ 2*(n+2) /)) ! associate C pointer returned by win_allocate to our fortran pointer

    ! Note: The RMA window has two extra slots for storing send/receive status signals

    A(:) = 0.d0
    
    PRINT*,'MPI initialization done..'

END SUBROUTINE rma_init



END PROGRAM  test5_onesided
