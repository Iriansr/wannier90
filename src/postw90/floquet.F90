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
        use w90_utility

        !Placeholder.
        complex(kind=dp), dimension(num_wann, num_wann) :: HH_k_t, debug2
        real(kind=dp), dimension(num_wann) :: debug1

        !Test showing that the spectra of both Hamiltonians are not the same even with no t-dep driving.
        !However, they yiedl the same unitary operator on exponentiation.
        call tdep_hamiltonian((/0.2_dp,0.3d0,0.4d0/), 0.0d0, 0.d0, HH_k_t)
            call utility_diagonalize(HH_k_t,num_wann,debug1,debug2)
            print*, debug1!(1)-debug1(2)
        call t_evol_op((/0.2_dp,0.3d0,0.4d0/), 0.0d0, 1.d0, 1000, HH_k_t)
        HH_k_t = utility_logu(HH_k_t,num_wann)/twopi
            call utility_diagonalize(HH_k_t,num_wann,debug1,debug2)
            print*, debug1!(1)-debug1(2)

    end subroutine floquet_main

    subroutine get_q(t, omega, q) !Placeholder.
        use w90_constants, only: dp, twopi
        real(kind=dp), intent(in) :: t, omega
        real(kind=dp), intent(out), dimension(3) :: q
        q = 0.d0!sin(omega*t)
    end subroutine get_q

    subroutine get_dq(t, omega, dq) !Placeholder.
        use w90_constants, only: dp, twopi
        real(kind=dp), intent(in) :: t, omega
        real(kind=dp), intent(out), dimension(3) :: dq
        dq = 0.d0!cos(omega*t)
    end subroutine get_dq

    subroutine t_evol_op(kpt, t0, omega, ntpts, U_k)

        use w90_constants, only: dp, cmplx_0, twopi
        use w90_parameters, only : num_wann
        use w90_utility, only : utility_expsh

        complex(kind=dp), intent(out), dimension(num_wann, num_wann) :: U_k
        real(kind=dp), intent(in), dimension(3) :: kpt
        real(kind=dp), intent(in) :: t0, omega
        integer, intent(in) :: ntpts

        complex(kind=dp), dimension(num_wann, num_wann) :: HH_k_t
        real(kind=dp) :: t
        integer :: i, it

        !Initialization.
        U_k = cmplx_0
        do i = 1, num_wann
            U_k(i,i) = 1.0_dp
        enddo

        do it = 1, ntpts-1
            t = t0 + twopi*real(it-1,dp)/(omega*real(ntpts-1,dp))
            call tdep_hamiltonian(kpt, t, omega, HH_k_t)
            !Get time evolution operator U(t, t + delta t).
            HH_k_t = utility_expsh(twopi*HH_k_t/(omega*real(ntpts-1,dp)),num_wann)
            U_k = matmul(U_k, HH_k_t)
        enddo

    end subroutine t_evol_op

    subroutine tdep_hamiltonian(kpt, t, omega, HH_k_t)

        use w90_constants, only : dp, cmplx_0, cmplx_i
        use w90_postw90_common, only : nrpts, rpt_origin, irvec, pw90common_fourier_R_to_k_new
        use w90_parameters, only : num_wann
        use w90_get_oper, only : HH_R, AA_R, get_HH_R, get_AA_R

        real(kind=dp), intent(in) :: t, omega
        real(kind=dp), intent(in), dimension(3) :: kpt
        complex(kind=dp), intent(out), dimension(num_wann, num_wann) :: HH_k_t

        real(kind=dp), dimension(3) :: q, dq, r_ij
        real(kind=dp), dimension(num_wann, 3) :: wan_centres
        complex(kind=dp), dimension(num_wann, num_wann, nrpts) :: HH_R_t
        integer :: i, j, ir

        call get_HH_R
        call get_AA_R

        call get_q(t, omega, q)
        call get_dq(t, omega, dq)

        !Get Wannier centres.
        do i = 1, num_wann
            wan_centres(i,:) = AA_R(i, i, rpt_origin, :)
        enddo

        do ir = 1, nrpts
            do i = 1, num_wann
                do j = 1, num_wann
                    r_ij = wan_centres(i,:) - wan_centres(j,:)
                    !Compute time-dependent Hamiltonian in the WF basis.
                    if (sqrt(dot_product(r_ij,r_ij))<1.d-3) then
                        HH_R_t(j, i, ir) = HH_R(j, i, ir)
                    else
                        HH_R_t(j, i, ir) = HH_R(j, i, ir) + dot_product(dq, AA_R(j, i, ir, :))
                        HH_R_t(j, i, ir) = HH_R_t(j, i, ir)*exp(cmplx_i*dot_product(q,r_ij))
                    endif
                enddo
            enddo
        enddo

        call pw90common_fourier_R_to_k_new(kpt, HH_R_t, HH_k_t)

    end subroutine tdep_hamiltonian

end module w90_floquet