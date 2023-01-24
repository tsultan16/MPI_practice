PROGRAM mpiHello2

USE MPI

IMPLICIT NONE

INTEGER :: err_num, this_process, num_process, i
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status


! initialize mpi
CALL MPI_INIT(err_num)  
! get total # of proceses
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_process, err_num)
! get # of this MPI rank
CALL MPI_COMM_RANK(MPI_COMM_WORLD, this_process, err_num)

IF(this_process .EQ. 0) THEN
    ! rank 0 is the head process 
    PRINT*,'Hello! This is process#',this_process,' of ',num_process,' processes.'
   
    ! receive message from all the other processes
    DO i=1, num_process-1
        CALL MPI_RECV(this_process, 1, MPI_INTEGER, i, 1, MPI_COMM_WORLD, status, err_num)
        PRINT*,'Received message from process#',i,' of ',num_process,' processes.'
    END DO

ELSE

     ! send message to head process
     CALL MPI_SEND(this_process, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, err_num)
END IF

! terminate MPI
CALL MPI_FINALIZE(err_num)




END PROGRAM mpiHello2
