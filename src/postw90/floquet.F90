module w90_floquet
  
    use w90_constants, only: dp
  
    implicit none
  
    private
  
    public :: floquet_main

    contains

    subroutine floquet_main

        use w90_constants
        use w90_comms
        use w90_io
        use w90_postw90_common
        use w90_parameters
        use w90_get_oper

        !Placeholder.
        open(unit=111,action="write",status="unknown",file="testfloquet.dat")
        write(unit=111,fmt=*) num_wann, "Éxito!"
        close(unit=111)

    end subroutine floquet_main

end module w90_floquet