module w90_floquet
  
    use w90_constants, only: dp
  
    implicit none

    private
  
    public :: floquet_main

    contains

    subroutine floquet_main

        use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
        use w90_postw90_common, only: nrpts, irvec
        use w90_parameters, only: num_kpts, kpt_latt, num_wann, omega_floq, t0, ntpts,&
        floq_calc_st, overwrite_floquet

        use w90_get_oper, only : HH_R, AA_R, get_HH_R, get_AA_R

        complex(kind=dp), dimension(num_wann, num_wann)    :: H_F_k
        complex(kind=dp), dimension(num_wann, num_wann, 3) :: A_F_k

        integer                                            :: ir, ik

        if (overwrite_floquet) then

            open(unit = 343, action = "write", file = "FloquetStatus.dat", status = "unknown")

            !Get the effective Floquet Hamiltonian on the WF basis using the ab-initio kmesh.
            HH_R = cmplx_0
            do ik = 1, num_kpts!Compute for each k.
                call eff_floquet_hamil(kpt_latt(:, ik), omega_floq, H_F_k)
                write(unit = 343, fmt = *), "Calculating H_F(R), progress:", int(100*real(ik,dp)/num_kpts), "%"
                do ir = 1, nrpts!Inregrate over k.
                    HH_R(:, :, ir) = HH_R(:, :, ir) + exp(-cmplx_i*twopi*dot_product(kpt_latt(:, ik), irvec(:, ir)))*H_F_k
                enddo
            enddo
            HH_R = HH_R/real(num_kpts, dp)

            !Get the Berry connection A^F_ij(R) on the Floquet gauge and the WF basis.
            AA_R = cmplx_0
            do ik = 1, num_kpts!Compute for each k.
                call get_berry_connection_on_gauge(kpt_latt(:, ik), HH_R, A_F_k)
                write(unit = 343, fmt = *), "Calculating A_F(R), progress:", int(100*real(ik,dp)/num_kpts), "%"
                do ir = 1, nrpts!Inregrate over k.
                    AA_R(:, :, ir, :) = AA_R(:, :, ir, :) + exp(-cmplx_i*twopi*dot_product(kpt_latt(:, ik), irvec(:, ir)))*A_F_k
                enddo
            enddo
            AA_R = AA_R/real(num_kpts, dp)

            close(unit = 343)

            floq_calc_st = .true. !Signal that effective Floquet quantities have been calculated.

        else

            !!!PLACE TO TEST THINGS

        endif

    end subroutine floquet_main

    subroutine get_q(t, omega, q) !Placeholder for driving force.
        use w90_constants, only: dp, twopi
        real(kind=dp), intent(in) :: t, omega
        real(kind=dp), intent(out), dimension(3) :: q
        q = 1.0_dp!0.0_dp!sin(omega*t)
    end subroutine get_q

    subroutine get_berry_connection_on_gauge(kpt, H_R, A_k)

        !=======================================================!
        !                                                       !
        !Subroutine to obtain the Berry connection A_k on the   !
        !k-point kpt and the gauge that diagonalizes the        !
        !Hamiltonian H_R given in the WF basis.                 !
        !                                                       !
        !=======================================================!

        use w90_constants, only: dp, cmplx_0, cmplx_i
        use w90_parameters, only : num_wann
        use w90_postw90_common, only: pw90common_fourier_R_to_k_new, pw90common_fourier_R_to_k_vec, nrpts
        use w90_utility, only: utility_diagonalize, utility_rotate
        use w90_get_oper, only : AA_R, get_AA_R
        use w90_wan_ham, only: wham_get_D_h_P_value

        real(kind=dp),    intent(in),  dimension(3)                         :: kpt
        complex(kind=dp), intent(in),  dimension(num_wann, num_wann, nrpts) :: H_R
        complex(kind=dp), intent(out), dimension(num_wann, num_wann, 3)     :: A_k

        complex(kind=dp), dimension(num_wann, num_wann)                     :: H_k, U_k
        complex(kind=dp), dimension(num_wann, num_wann, 3)                  :: delH_k, D_k
        real(kind=dp),    dimension(num_wann)                               :: eig

        integer                                                             :: i

        !Get Berry conn. on WF basis and WF gauge.
        call get_AA_R
        call pw90common_fourier_R_to_k_vec(kpt, AA_R, A_k)

        !Get the Hamiltonian and its derivatives on the point kpt and diagonalize it.
        call pw90common_fourier_R_to_k_new(kpt, H_R, H_k, delH_k(:, :, 1), delH_k(:, :, 2), delH_k(:, :, 3))
        call utility_diagonalize(H_k, num_wann, eig, U_k)
        !Compute D_a(k)=U_k^dag.del_a(U_k) (a=x,y,z) using Eq.(24) of WYSV06 and prescription for energy denominator from BK81.
        call wham_get_D_h_P_value(delH_k, U_k, eig, D_k)
        !Compute the Berry connection on the desired gauge (using transformation properties).
        do i = 1, 3
            A_k(:, :, i) = utility_rotate(A_k(:, :, i), U_k, num_wann) + cmplx_i*D_k(:, :, i)
        enddo
    end subroutine get_berry_connection_on_gauge

    subroutine eff_floquet_hamil(kpt, omega, H_F_k, eig_F, UU_F_k)

        !=======================================================!
        !                                                       !
        !Subroutine to obtain the effective Floquet hamiltonian !
        !H_F_k on the k-point kpt and the Wannier gauge.        !
        !Optionally returns the quasienergy spectra, calculated !
        !using convergence prescription, and the unitary matrix !
        !that transforms from Wannier to Floquet gauge.         !
        !                                                       !
        !=======================================================!

        use w90_constants, only: dp, cmplx_0, twopi
        use w90_parameters, only : num_wann, floquet_conv_factor, ntpts, t0
        use w90_utility, only : utility_exphs, utility_logh, utility_diagonalize
        use w90_get_oper, only : HH_R, AA_R, get_HH_R, get_AA_R
        use w90_postw90_common, only : pw90common_fourier_R_to_k_new, pw90common_fourier_R_to_k_vec

        complex(kind=dp), intent(out), dimension(num_wann, num_wann)           :: H_F_k
        complex(kind=dp), intent(out), dimension(num_wann, num_wann), optional :: UU_F_k
        real(kind=dp)   , intent(out), dimension(num_wann), optional           :: eig_F
        real(kind=dp)   , intent(in) , dimension(3)                            :: kpt
        real(kind=dp)   , intent(in)                                           :: omega

        complex(kind=dp), dimension(num_wann, num_wann)                        :: HH_k_t
        real(kind=dp)                                                          :: t
        integer                                                                :: it

        !Initialization.
        H_F_k = cmplx_0
        do it = 1, num_wann
            H_F_k(it,it) = 1.0_dp
        enddo

        do it = 1, ntpts-1
            t = t0 + twopi*real(it-1,dp)/(omega*real(ntpts-1,dp))
            call tdep_hamiltonian(kpt, t, omega, HH_k_t)
            !Get time evolution operator U(t, t + delta t) using convergence prescription,
            !this is equivalent to perform a Suzuki–Trotter expansion.
            HH_k_t = utility_exphs(floquet_conv_factor*HH_k_t/real(ntpts-1,dp),num_wann, .false.)
            !Accumulate for all time slices.
            H_F_k = matmul(H_F_k, HH_k_t)
        enddo

        !Take log with numerical convergence prescription.
        H_F_k = utility_logh(H_F_k,num_wann)/floquet_conv_factor

        !Get the quasienergy spectra eig_F and the unitary matrix that 
        !transforms from Wannier gauge to Floquet gauge.
        if (present(UU_F_k).AND.present(eig_F)) then
            call utility_diagonalize(H_F_k, num_wann, eig_F, UU_F_k)
        endif

    end subroutine eff_floquet_hamil

    subroutine get_Qs(kpt, omega, Qs)

        !=======================================================!
        !                                                       !
        !Subroutine to obtain the Fourier amplitudes of the     !
        !operator P(t), Qs, on the k-point kpt and the Wannier  !
        !gauge.                                                 !                                  
        !                                                       !
        !=======================================================!

        use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
        use w90_parameters, only : num_wann, srange, ntpts, t0
        use w90_utility, only : utility_exphs, utility_logh, utility_diagonalize
        use w90_get_oper, only : HH_R, AA_R, get_HH_R, get_AA_R
        use w90_postw90_common, only : pw90common_fourier_R_to_k_new, pw90common_fourier_R_to_k_vec

        complex(kind=dp), intent(out), dimension(num_wann, num_wann, srange)   :: Qs

        real(kind=dp)   , intent(in) , dimension(3)                            :: kpt
        real(kind=dp)   , intent(in)                                           :: omega

        complex(kind=dp), dimension(num_wann, num_wann)                        :: H_F_K, TEV, AUX
        complex(kind=dp), dimension(num_wann, num_wann, srange)                :: Pt
        real(kind=dp)                                                          :: t
        integer                                                                :: it, is

        !Initialization.
        Qs = cmplx_0
        Pt = cmplx_0
        H_F_K = cmplx_0
        call eff_floquet_hamil(kpt, omega, H_F_k)

        do it = 1, ntpts-1
            t = t0 + twopi*real(it-1,dp)/(omega*real(ntpts-1,dp))
            !Get time evolution operator U(t0, t0 + it* delta t).
            call tdep_hamiltonian(kpt, t, omega, AUX)
            AUX = utility_exphs(cmplx_i*twopi*AUX/(omega*real(ntpts-1,dp)),num_wann, .true.)
            !Accumulate for all time slices.
            TEV = matmul(TEV, AUX)
            !Calcualte P(t) for each time slice.
            Pt(:, :, it) = TEV*utility_exphs(cmplx_i*t*H_F_K,num_wann, .true.)
        enddo

        !Perform a Fourier transform to calculate the Fourier amplitudes of P(t).
        do is = 1, srange
            do it = 1, ntpts-1
                t = t0 + twopi*real(it-1,dp)/(omega*real(ntpts-1,dp))
                Qs(:, :, is) = Qs(:, :, is) + Pt(:, :, it)*exp(cmplx_i*is*omega*t)
            enddo
            Qs(:, :, is) = omega*Qs(:, :, is)/twopi !Divide by period.
        enddo

    end subroutine get_Qs

    subroutine tdep_hamiltonian(kpt, t, omega, HH_k_t)

        !Get the time-dependent Hamiltonian H(t) = H + q(t)*r on the Wannier gauge.

        use w90_constants, only : dp, cmplx_0, cmplx_i
        use w90_postw90_common, only : nrpts, rpt_origin, irvec, pw90common_fourier_R_to_k_new
        use w90_parameters, only : num_wann
        use w90_get_oper, only : HH_R, AA_R, get_HH_R, get_AA_R

        real(kind=dp),    intent(in)                                 :: t, omega
        real(kind=dp),    intent(in),  dimension(3)                  :: kpt
        complex(kind=dp), intent(out), dimension(num_wann, num_wann) :: HH_k_t

        real(kind=dp),    dimension(3)                               :: q, dq, r_ij
        real(kind=dp),    dimension(num_wann, 3)                     :: wan_centres
        complex(kind=dp), dimension(num_wann, num_wann, nrpts)       :: HH_R_t
        integer :: i, j, ir

        call get_HH_R
        call get_AA_R

        call get_q(t, omega, q)

        do ir = 1, nrpts
            do i = 1, num_wann
                do j = 1, num_wann
                    HH_R_t(j, i, ir) = HH_R(j, i, ir) + dot_product(q, AA_R(j, i, ir, :))
                enddo
            enddo
        enddo

        call pw90common_fourier_R_to_k_new(kpt, HH_R_t, HH_k_t)

    end subroutine tdep_hamiltonian

end module w90_floquet