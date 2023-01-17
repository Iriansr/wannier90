module w90_floquet
  
    use w90_constants, only: dp
  
    implicit none

    private
  
    public :: floquet_main, floquet_get_FT_obs_A

    contains

    subroutine floquet_main

        use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
        use w90_postw90_common, only: nrpts, irvec
        use w90_parameters, only: num_kpts, kpt_latt, num_wann, omega_floq, t0, ntpts,&
        floq_calc_st, overwrite_floquet, floquet_conv_factor

        use w90_get_oper, only : HH_R, AA_R, get_HH_R, get_AA_R

        !!!!TODO: CHECK ALL OMEGA DEPENDENCE OF QUANTITIES AND CORRECT THE UNITS USING EV_TO_SECONDS.

        !!!!TESTEBED.
        use w90_utility, only : utility_exphs, utility_logh, utility_diagonalize
        complex(kind=dp), dimension(2, 2)    :: mat, rot
        real(kind=dp), dimension(2) :: eig

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
            mat = cmplx_0
            rot = cmplx_0
            eig = 0.0_dp

            !HERMITIAN
            mat(1,:) = (/cmplx(2.0_dp, 0.0_dp), cmplx(0.0_dp, 1.0_dp)/)
            mat(2,:) = (/cmplx(0.0_dp, -1.0_dp), cmplx(2.0_dp, 0.0_dp)/)
            rot = utility_exphs(mat, 2, .false.)
            do ir = 1, 2
                do ik = 1, 2
                    print*, ir, ik, rot(ir, ik)
                enddo
            enddo

            print*, "log of exp:", mat, utility_logh(utility_exphs(mat, 2, .false.), 2)

            !SKEW-HERMITIAN
            mat(1,:) = cmplx_i*(/cmplx(2.0_dp, 0.0_dp), cmplx(0.0_dp, 1.0_dp)/)
            mat(2,:) = cmplx_i*(/cmplx(0.0_dp, -1.0_dp), cmplx(2.0_dp, 0.0_dp)/)
            rot = utility_exphs(mat, 2, .true.)
            print*, "unitary check", matmul(rot, conjg(transpose(rot)))
            do ir = 1, 2
                do ik = 1, 2
                    print*, ir, ik, rot(ir, ik)
                enddo
            enddo

        endif

    end subroutine floquet_main

    subroutine floquet_get_FT_obs_A(kpt, A, u_target, FTA, error)

        !=======================================================!
        !                                                       !
        !Subroutine to compute the Fourier transform FTA at     !
        !frequency u_target*omega of the expected value of the  !
        !observable A at k-point kpt.                           !
        !The observable A is given by a complex matrix          !
        !containing its matrix elements in the Hamiltonian      !
        !gauge for t=0.                                         !
        !FTA is a complex number.                               !
        !                                                       !
        !WARNING: "A" must transform covariantly under gauge    !
        !transformations, never as a connection, so this        !
        !routine cannot calculate the FT of the pos. operator.  !
        !                                                       !
        !=======================================================!

        use w90_constants, only: dp, cmplx_0, cmplx_i, pi, twopi
        use w90_parameters, only: num_wann, srange, nfermi, kubo_nfreq, kubo_freq_list, &
        fermi_energy_list, kubo_smr_index, kubo_smr_fixed_en_width, floquet_conv_factor, t0, ntpts
        use w90_utility, only: utility_rotate, utility_diagonalize, utility_w0gauss
        use w90_get_oper, only: HH_R, get_HH_R
        use w90_postw90_common, only: pw90common_fourier_R_to_k_new, pw90common_get_occ

        real(kind=dp),    intent(in),  dimension(3)                                  :: kpt
        real(kind=dp),    intent(in)                                                 :: u_target
        complex(kind=dp), intent(in),  dimension(:, :)                               :: A

        complex(kind=dp), intent(out), dimension(kubo_nfreq)                         :: FTA, error

        real(kind=dp),                 dimension(num_wann)                           :: eig_H, occ_H, &
                                                                                     eig_F
        complex(kind=dp),              dimension(num_wann, num_wann)                 :: HH_k, WtoH, occ_mat_H, &
                                                                                     H_F_k, WtoF, FA, occ_mat_F

        complex(kind=dp),              dimension(:, :, :), allocatable   ::    Qs, FQs
        complex(kind=dp),              dimension(:, :, :), allocatable   ::    Pt

        real(kind=dp),                 dimension(kubo_nfreq)                         :: delta, p_value
        complex(kind=dp),              dimension(kubo_nfreq)                         :: RFTA, IFTA
        complex(kind=dp)                                                             :: f_integrand   
        real(kind=dp)                                                                :: t                           

        integer                                                                      :: n, m, l, p, r, s, ifreq, is, it

        !!! This stays here until I figure out how to solve losing variable scope in parallel. 
        kubo_smr_fixed_en_width = 0.025_dp
        ntpts = 100
        srange = 10
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        allocate(Pt(num_wann, num_wann, ntpts))
        allocate(Qs(num_wann, num_wann, -srange:srange))
        allocate(FQs(num_wann, num_wann, -srange:srange))

        FTA = cmplx_0
        error = cmplx_0

        occ_mat_H = cmplx_0
        FA = cmplx_0
        Qs = cmplx_0
        FQs = cmplx_0
        H_F_k = cmplx_0

        !Get undriven Hamiltonian matrix at kpt and Wannier gauge.
        call pw90common_fourier_R_to_k_new(kpt, HH_R, HH_k)
        !Diagonalize and get eigenvalues and Wannier to Hamiltonian gauge transformation.
        call utility_diagonalize(HH_k, num_wann, eig_H, WtoH)
        !Setup occ. factors and obtain density matrix occ_mat_H.
        call pw90common_get_occ(eig_H, occ_H, fermi_energy_list(1))
        do n = 1, num_wann
            occ_mat_H(n, n) = cmplx(occ_H(n), 0.0, dp)
        enddo

        do ifreq = 1, kubo_nfreq !For each frequency:

            !Get the effective floquet Hamiltonian and the Wannier to Floquet gauge transformation.
            call get_tev(kpt, real(kubo_freq_list(ifreq), dp), H_F_K, Pt)
            !Get the Qs-s on the Wannier gauge.

            do is = -srange, srange
                do it = 1, ntpts-1
                    t = t0 + twopi*real(it-1,dp)/(real(kubo_freq_list(ifreq), dp)*real(ntpts-1,dp))
                    !print*, it, t, ntpts, Pt(1, 1, it), Qs(1, 1, is)
                    Qs(:, :, is) = Qs(:, :, is) + Pt(:, :, it)*exp(cmplx_i*is*real(kubo_freq_list(ifreq), dp)*t)/real(ntpts-1,dp)
                enddo
            !TODO: CHECK IF QS NEEDS TO BE NORMALIZED.
            enddo

            call utility_diagonalize(H_F_k, num_wann, eig_F, WtoF)

            !Rotate all quantities to Floquet gauge.
            FA = utility_rotate(A, WtoF, num_wann)
            occ_mat_F = utility_rotate(occ_mat_H, WtoF, num_wann)
            do s = -srange, srange
                FQs(:, :, s) = utility_rotate(Qs(:, :, s), WtoF, num_wann)
            enddo

            do r = -srange, srange
                do s = -srange, srange

                    do l = 1, num_wann
                        do p = 1, num_wann

                            !Compute delta and p-value.
                            delta(ifreq) = utility_w0gauss((s - r + u_target + ((eig_F(l) - eig_F(p))/real(kubo_freq_list(ifreq), dp)))/kubo_smr_fixed_en_width, kubo_smr_index)/kubo_smr_fixed_en_width

                            if ((s - r + u_target + ((eig_F(l) - eig_F(p))/real(kubo_freq_list(ifreq), dp))) .lt. 1.0E-6_dp) then
                                p_value(ifreq) = 0.0_dp
                            else
                                p_value(ifreq) = 1.0_dp/((s - r + u_target + ((eig_F(l) - eig_F(p))/real(kubo_freq_list(ifreq), dp))))
                            endif

                            do n = 1, num_wann
                                do m = 1, num_wann

                                    f_integrand = occ_mat_F(n,m)*FA(l,p)*conjg(FQs(l,m,s))*FQs(p,n,r)
                                    RFTA(ifreq) = RFTA(ifreq) + pi*f_integrand*delta(ifreq)
                                    IFTA(ifreq) = IFTA(ifreq) - f_integrand*p_value(ifreq)

                                enddo!m
                            enddo!n

                        enddo!p
                    enddo!l

                enddo!s
            enddo!r

            FTA(ifreq) = cmplx(real(RFTA(ifreq),dp), real(IFTA(ifreq),dp))
            error(ifreq) = cmplx(aimag(RFTA(ifreq)), aimag(IFTA(ifreq)))

        enddo!ifreq

    end subroutine floquet_get_FT_obs_A

    subroutine get_q(t, omega, q) !Placeholder for driving force.
        use w90_constants, only: dp, twopi
        real(kind=dp), intent(in) :: t, omega
        real(kind=dp), intent(out), dimension(3) :: q
        q(1) = sin(omega*t)
        q(2) = sin(omega*t) !RECALL: CAREFUL WITH EXTREME AMPLITUDES. q has units of eV/A = 10^10 ev/m.
        q(3) = 0.0_dp
        !q = q*1E-40_dp
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
        !Optionally returns the quasienergy spectra and the     !
        !unitary matrix that transforms from Wannier to Floquet !
        !gauge.                                                 !
        !                                                       !
        !=======================================================!

        use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
        use w90_parameters, only : num_wann, floquet_conv_factor, ntpts, t0
        use w90_utility, only : utility_exphs, utility_logh, utility_diagonalize
        use w90_get_oper, only : HH_R, AA_R, get_HH_R, get_AA_R
        use w90_postw90_common, only : pw90common_fourier_R_to_k_new, pw90common_fourier_R_to_k_vec

        complex(kind=dp), intent(out), dimension(num_wann, num_wann)           :: H_F_k
        complex(kind=dp), intent(out), dimension(num_wann, num_wann), optional :: UU_F_k
        real(kind=dp)   , intent(out), dimension(num_wann)          , optional :: eig_F
        real(kind=dp)   , intent(in) , dimension(3)                            :: kpt
        real(kind=dp)   , intent(in)                                           :: omega

        complex(kind=dp)             , dimension(num_wann, num_wann)           :: HH_k_t
        real(kind=dp)                                                          :: t
        integer                                                                :: it

        !!! This stays here until I figure out what causes to set this variables to 0 when this subroutine runs on n-1 processors. 
        floquet_conv_factor = 0.01_dp
        ntpts = 100
        t0 = 0.0_dp
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

    subroutine get_tev(kpt, omega, H_F_K, Pt)

        !========================================================!
        !                                                        !
        !Subroutine to decompose the time-evolution operator     !
        !into the Floquet form                                   !
        !TEV = \sum_{s} Qs exp(-i*s*\omega*t) exp(-i*H_F*t/\hbar)!
        !on the k-point kpt and the Wannier gauge.               !                                  
        !                                                        !
        !========================================================!

        use w90_constants, only: dp, cmplx_0, cmplx_i, twopi
        use w90_parameters, only : num_wann, srange, ntpts, t0
        use w90_utility, only : utility_exphs, utility_logh, utility_diagonalize

        complex(kind=dp), intent(inout), dimension(:, :)                   :: H_F_K
        complex(kind=dp), intent(inout), dimension(:, :, :)  :: Pt

        real(kind=dp)   , intent(in) , dimension(3)                                    :: kpt
        real(kind=dp)   , intent(in)                                                   :: omega

        complex(kind=dp), dimension(:, :, :), allocatable                        :: TEV
        complex(kind=dp), dimension(:, :), allocatable                                :: AUX
        real(kind=dp)                                                                  :: t, dt
        integer                                                                        :: it, is, k

        !NEED TO BE TRIED IN PARALLEL.
        !DUE TO THE BELOW: PASS NTPTS, SRANGE... AS INPUT PARAMETERS RATHER THAN SAVE FROM PARAMETERS.F90

        ntpts = 100
        srange = 10
        t0 = 0.00_dp

        allocate(TEV(num_wann, num_wann, ntpts))
        allocate(AUX(num_wann, num_wann))

        !Initialization.
        H_F_K = cmplx_0 !!!
        !Qs = cmplx_0
        Pt = cmplx_0
        TEV   = cmplx_0 
        do k = 1, num_wann
            TEV(k, k, :) = cmplx(1.0_dp, 0.0_dp)
        enddo  
        AUX   = cmplx_0

        do k = 1, ntpts

            t = t0 + twopi*real(k-1,dp)/(omega*real(ntpts-1,dp))
            !Get time evolution operator U(t0, t).

            do it = 1, ntpts

                dt = (t-t0)/(real(ntpts-1,dp))
                call tdep_hamiltonian(kpt, t0 + real(it-1, dp)*dt, omega, AUX)
                AUX = utility_exphs(-cmplx_i*dt*AUX,num_wann, .true.)
                TEV(:, :, k) = matmul(TEV(:, :, k), AUX)

            enddo
        enddo

        H_F_k = (cmplx_i*omega*utility_logh(TEV(:, :, ntpts),num_wann)/twopi)

        do it = 1, ntpts
            t = t0 + twopi*real(it-1,dp)/(omega*real(ntpts-1,dp))
            Pt(:, :, it) = matmul(TEV(:, :, it),utility_exphs(cmplx_i*H_F_k*t, num_wann, .true.))
        enddo

    end subroutine get_tev

    subroutine tdep_hamiltonian(kpt, t, omega, HH_k_t)

        !Get the time-dependent Hamiltonian H(t) = H + q(t)*r on the Wannier gauge.

        use w90_constants, only : dp, cmplx_0, cmplx_i
        use w90_postw90_common, only : nrpts, pw90common_fourier_R_to_k_new
        use w90_parameters, only : num_wann
        use w90_get_oper, only : HH_R, AA_R, get_HH_R, get_AA_R

        real(kind=dp),    intent(in)                                 :: t, omega
        real(kind=dp),    intent(in),  dimension(3)                  :: kpt
        complex(kind=dp), intent(out), dimension(num_wann, num_wann) :: HH_k_t

        real(kind=dp),    dimension(3)                               :: q
        complex(kind=dp), dimension(num_wann, num_wann, nrpts)       :: HH_R_t
        integer                                                      :: i, j, ir

        HH_k_t = cmplx_0

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