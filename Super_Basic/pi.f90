INCLUDE  'timing_mod.f90'

! Numerical integration of f(x)=4/(1+x^2) in x=[0,1]
! using midpoint rule.
!
! Parallelized using MPI: Each MPI rank computes part of the total sum
! Then an MPI_REDUCE adds up the partial sums from all the ranks to
! get the final result.

PROGRAM compute_pi

USE timing_mod
USE MPI
IMPLICIT NONE

INTEGER :: i,j,N
REAL*8 :: dx, x, pi_fortran_internal
REAL*8 :: partial_pi, partial_sum, total_pi
INTEGER :: this_process, n_process, err_num 


! initialize mpi
CALL MPI_INIT(err_num)  
! get total # of proceses
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, n_process, err_num)
! get # of this MPI rank
CALL MPI_COMM_RANK(MPI_COMM_WORLD, this_process, err_num)

N = 100000
pi_fortran_internal = 4.0_8*atan(1.0_8)

IF(this_process .EQ. 0) THEN
    PRINT*,'Fortran internal pi = ',pi_fortran_internal
    PRINT*,''
    CALL start_timer()
END IF

DO j = 1,5
    dx = 1.0_8/N
    partial_sum = 0.0_8

    ! have the nth MPI rank accumulate the partial sum from every (n_process)th term
    ! starting from 1+n
    DO i = 1+this_process,N,n_process
        x = dx*(REAL(i,8)-0.5_8)
        partial_sum = partial_sum + f(x)
    END DO
    partial_pi  = dx*partial_sum
    
    ! add all the partial sums to get final result using
    ! MPI_REDUCE(send_buffer, recv_buffer, num_buffer_elements, buffer_data_type, reduce_operation, root_process, communicator, err)
    CALL MPI_REDUCE(partial_pi, total_pi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, err_num)
  
    ! have root process print out the final results
    IF(this_process .EQ. 0) THEN
        PRINT*,'N = ',N
        PRINT*,'pi = ',total_pi
        PRINT*,'pi-pi_fortran =', ABS(total_pi-pi_fortran_internal)
        PRINT*,'Time(sec) = ',time_difference()
        PRINT*,''
    END IF
    N = N*10
END DO

CALL MPI_FINALIZE(err_num)


CONTAINS

FUNCTION f(x) result(fx)

    REAL*8, INTENT(IN) :: x
    REAL*8 :: fx

    fx = 4.0_8/(1.0_8+x*x) 

END FUNCTION f


END PROGRAM compute_pi


