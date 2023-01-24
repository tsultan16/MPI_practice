PROGRAM mpiHello

USE MPI

IMPLICIT NONE


INTEGER :: err_num, this_process, num_process

! initialize mpi
CALL MPI_INIT(err_num)  
! get total # of proceses
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_process, err_num)
! get # of this MPI rank
CALL MPI_COMM_RANK(MPI_COMM_WORLD, this_process, err_num)

PRINT*,'Hello! This is process#',this_process,' of ',num_process,' processes.'

! terminate MPI
CALL MPI_FINALIZE(err_num)


END PROGRAM mpiHello
