MODULE timing_mod

IMPLICIT NONE

INTEGER, PRIVATE :: dt(8) 
REAL, PRIVATE :: h, m, s, ms, tt, last_tt, initial_tt

CONTAINS

SUBROUTINE start_timer()

    CALL date_and_time(values = dt)
    h = REAL(dt(5))
    m = REAL(dt(6))
    s = REAL(dt(7))
    ms = REAL(dt(8))
    last_tt = 60*(60*h+m)+s+ms/1000.0
    initial_tt = last_tt

END SUBROUTINE start_timer


FUNCTION time_difference() RESULT(tdif)

    REAL :: tdif

    tt = 0.0
    CALL date_and_time(values = dt)
    h = REAL(dt(5))
    m = REAL(dt(6))
    s = REAL(dt(7))
    ms = REAL(dt(8))
    tt = 60*(60*h+m)+s+ms/1000.0
    tdif = tt - last_tt
    last_tt = tt

END FUNCTION time_difference

SUBROUTINE end_timer()

    CALL date_and_time(values = dt)
    h = REAL(dt(5))
    m = REAL(dt(6))
    s = REAL(dt(7))
    ms = REAL(dt(8))
    tt = 60*(60*h+m)+s+ms/1000.0

    PRINT*,'Total time elapsed (sec) = ', tt-initial_tt

END SUBROUTINE end_timer

END MODULE timing_mod
