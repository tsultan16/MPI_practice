! This example program demonstrates data  exchange between 2 processes
! using one-sided RMA. 
! We use passive target RMA and acheive synchronization via MPI_WIN_FLUSH. 

PROGRAM test4_onesided

USE ISO_C_BINDING

IMPLICIT NONE

INCLUDE 'mpif.h'


!*********************************************************************************************************************************

INTEGER :: comm1d, ierr, ndim, myrank, numprocs(1), coord(1), req(4), status(MPI_STATUS_SIZE), source, tag
LOGICAL :: isperiodic(1), reorder 
INTEGER, PARAMETER :: n = 5
REAL(8), ALLOCATABLE :: u(:), buffer(:)
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


! put some test values inside u
IF(coord(1) .EQ. 0) THEN
    u(:) = (/ 1,2,3,4,5 /)     

    PRINT*,''
    PRINT*,'Myrank = ',myrank
    PRINT*,''
    PRINT*,'u =',u

    
END IF

CALL MPI_BARRIER(ierr)

IF(coord(1) .EQ. 1) THEN
    
    u(:) = (/ 11,12,13,14,15 /)     
    
    PRINT*,''
    PRINT*,'Myrank = ',myrank
    PRINT*,''
    PRINT*,'u =',u
    
END IF

CALL MPI_BARRIER(ierr)


! open RMA exposure epoch
CALL open_exposure_epoch()

! tell neighbor that we are ready to receive
CALL recv_signal()

! poll on signal from from neighbor telling us they are ready to receive
CALL poll_signal(A,1)

! send data to neighbor since they are ready
CALL send_msg()

! tell neighbor that we hare done sending
CALL data_sent_signal()

! poll on signal from from neighbor telling us they have sent us data_sent_signal
CALL poll_signal(A,2)


PRINT*,''
PRINT*,'One sided sommunication complete...Myrank = ',myrank
PRINT*,''

CALL MPI_BARRIER(ierr)


! unpack data from recv buffer
!u(1:n) = A(3:n+2)



!CALL MPI_BARRIER(ierr)


IF(coord(1) .EQ. 0) THEN

PRINT*,''
PRINT*,'Myrank = ',myrank
PRINT*,''
PRINT*,'RMA WINDOW=',A
PRINT*,''

END IF

CALL MPI_BARRIER(ierr)


IF(coord(1) .EQ. 1) THEN

PRINT*,''
PRINT*,'Myrank = ',myrank
PRINT*,''
PRINT*,'RMA WINDOW=',A
PRINT*,''

END IF

CALL MPI_BARRIER(ierr)

CALL close_rma_window()

PRINT*,'Done..myrank = ', myrank

DEALLOCATE(u, buffer)




CONTAINS



! send signal to neighbor indicating that we are ready to receive data
SUBROUTINE recv_signal()

 INTEGER :: ierr

    buffer(1) = 1.d0
    disp = 0
    
    PRINT*,'Myrank = ',myrank
    PRINT*,'BUFFER=',buffer
    
    ! put a single value (1.0) at the start of neighbors RMA window. This indicates
    ! that we are ready to receive a message from them.
    CALL MPI_PUT(buffer, 1, MPI_DOUBLE_PRECISION, neighbor, &
                 disp, 1, MPI_DOUBLE_PRECISION, win, ierr)    
    
    PRINT*,'PUT COMPLETED..Myrank = ',myrank


    CALL MPI_WIN_FLUSH(neighbor, win, ierr)

    PRINT*,'RECV signal sent. Myrank = ',myrank

END SUBROUTINE recv_signal


! sends message to neighbor
SUBROUTINE send_msg()

    INTEGER :: ierr
    

    disp = 0+2

    PRINT*,''
    PRINT*,'Sending msg...myrank = ',myrank
    PRINT*,''

    CALL MPI_PUT(u, n, MPI_DOUBLE_PRECISION, neighbor, &
                 disp, n, MPI_DOUBLE_PRECISION, win, ierr)
    
    CALL MPI_WIN_FLUSH(neighbor, win, ierr)


END SUBROUTINE send_msg


! send a signal to neighbor indicating that we have sent them data
SUBROUTINE data_sent_signal()

 INTEGER :: ierr


    
    buffer(1) = 1.d0
    disp = 0+1
    
    ! put a single value (1.0) at the start of neighbors RMA window. This indicates
    ! that we are ready to receive a message from them.
    CALL MPI_PUT(buffer, 1, MPI_DOUBLE_PRECISION, neighbor, &
                 disp, 1, MPI_DOUBLE_PRECISION, win, ierr)    
    
    CALL MPI_WIN_FLUSH(neighbor, win, ierr)


    PRINT*,'SEND signal sent. Myrank = ',myrank


END SUBROUTINE data_sent_signal


! polls on signal from neighbor
SUBROUTINE poll_signal(rma_window, signal_offset)

    REAL(8), VOLATILE :: rma_window(:) ! must include the volatile attribute here.. to ensure that compiler isn't using a copy of the RMA window that it has pre-loaded into a register
    INTEGER, INTENT(IN) :: signal_offset  ! 1 for recv signal, 2 for sent signal
    INTEGER :: ierr


    IF(signal_offset .EQ. 1) THEN
        PRINT*,'Polling for RECV signal...Myrank = ',myrank
    ELSE IF(signal_offset .EQ. 2) THEN
        PRINT*,'Polling for SENT signal...Myrank = ',myrank
    END IF
    
    PRINT*,'Myrank, A =',myrank, rma_window
    
    ! poll the signal status
    DO WHILE(rma_window(signal_offset) .NE. 1.d0)
    
        ! wait till signal has arrived...
        PRINT*,'Polling...MYRANK = ',myrank
    
    END DO

    ! reset the signal status
    rma_window(signal_offset) = 0.d0    

    IF(signal_offset .EQ. 1) THEN
        PRINT*,'RECV signal has arrived. Myrank = ',myrank
    ELSE IF(signal_offset .EQ. 2) THEN
        PRINT*,'SENT signal has arrived. Myrank = ',myrank
    END IF


END SUBROUTINE poll_signal



! start MPI passive RMA epoch
SUBROUTINE open_exposure_epoch()

    INTEGER :: ierr



    CALL MPI_WIN_LOCK(MPI_LOCK_SHARED, neighbor, MPI_MODE_NOCHECK, win, ierr)
    
    PRINT*,'Passive RMA epoch now exposed..'


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

    IF(np .NE. 2) THEN
        PRINT*,'Need 2 ranks for this test...'
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
    
    ! find my neighbor's rank
    IF(myrank .EQ. 0) THEN
        neighbor = 1
    ELSE IF(myrank .EQ. 1) THEN
        neighbor = 0
    END IF
    
    
    ! allocate work arrays
    ALLOCATE(u(n), buffer(n))
    u = 0.d0
    buffer = 0.d0

    ! determine array element size in bytes    
    CALL MPI_SIZEOF(DBLE(1), sizedouble, ierr)
    winsize = (n+2)*sizedouble

    IF(myrank .EQ. 0) PRINT*,'sizedouble, winsize = ',sizedouble, winsize 


    ! create rma local memory allocation and create window
    CALL MPI_WIN_ALLOCATE(winsize, sizedouble, MPI_INFO_NULL, comm1d, rma_cmem, win, ierr)
    CALL C_F_POINTER(rma_cmem, A, (/ n+2 /)) ! associate C pointer returned by win_allocate to our fortran pointer

    ! Note: The RMA window has two extra slots for storing send/receive status signals

    A(:) = 0.d0
    
    PRINT*,'MPI initialization done..  A = ',A

END SUBROUTINE rma_init



END PROGRAM  test4_onesided
