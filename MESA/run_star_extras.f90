! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use auto_diff
      
      implicit none

      real(dp), allocatable :: mdot_hl(:), mdot(:), mdot_edd(:), mdot_hyper(:), fd_hl(:), fd(:), edot_hl(:), edot(:), v(:), Ra(:), Eorb(:), dEorb(:), f(:)
      real(dp), allocatable :: eps_rho(:), mdot_mr15_ratio(:), fd_mr15_ratio(:), mdot_mr15(:), fd_mr15(:)
      real(dp), allocatable :: omega_arr(:), e_arr(:), beta_arr(:)
      real(dp) :: M_ns, eta, op_const, ebind, eorb_change, M_acc, R0, Req, omega, e, Q, Qmax, Qtb, h_val, Rbar, M_crust, beta
      real(dp) :: n_poly, beta_sec, a, b, mom_inert, D
      integer :: azone, stop
      
      ! these routines are called by the standard run_star check_model
      contains

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! OPERATIONAL FUNCTIONS !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !DONE: To find the closest lower value to a target in an array
      subroutine find_closest_value_index(arr, tar, n, closest_index)
         implicit none
         real(dp), intent(in) :: arr(:)      ! Array to search
         real(dp), intent(in) :: tar      ! Target value
         integer, intent(in) :: n           ! Size of the array
         integer, intent(out) :: closest_index  ! Index of the closest value
 
         integer :: l
         real(dp) :: min_diff, diff, temp
         
         ! print *, 'find closest value index called'
         ! Initialize closest_index and min_diff
         closest_index = 1
         min_diff = abs(arr(1) - tar)
 
         ! Loop through the array to find the closest value
         do l = 2, n
            !  temp = arr(l) - tar
             diff = abs(arr(l) - tar)
             if (diff < min_diff) then
                 min_diff = diff
                 closest_index = l
             end if
         end do
      end subroutine find_closest_value_index
      
      !DONE: To find array indices where array values lie in a given range
      subroutine find_inbetween_value_index(arr, v1, v2, n, ind_btw)

         real(dp), intent(in) :: arr(:)
         real(dp), intent(in) :: v1, v2
         integer, intent(in) :: n
         integer, allocatable, intent(out) :: ind_btw(:)
         
         integer :: l

         ! print *, 'find inbetween value index called'
         do l = 1, n
            if (arr(l) >= v1 .and. arr(l) <= v2) then
               call AddToList(ind_btw, l)
            end if
         end do
         ! print *, '*16*'
      end subroutine find_inbetween_value_index
      
      !DONE: To append elements to an allocatable array
      subroutine AddToList(list, element)
         integer :: i, isize
         integer, intent(in) :: element
         integer, dimension(:), allocatable, intent(inout) :: list
         integer, dimension(:), allocatable :: clist

         if(allocated(list)) then
             isize = size(list)
             allocate(clist(isize+1))
             do i=1,isize          
             clist(i) = list(i)
             end do
             clist(isize+1) = element

             deallocate(list)
             call move_alloc(clist, list)

         else
             allocate(list(1))
             list(1) = element
         end if
      end subroutine AddToList
      
      !DONE: To evaluate the function F for the differential equation da/dt = F(a)
      subroutine linspace(n, from, to, array)
         real(dp), intent(in) :: from, to
         integer, intent(in) :: n
         real(dp), allocatable, intent(out) :: array(:)
         real(dp) :: range
         integer :: i
         range = to - from
  
         if (n == 0) return
         allocate(array(n))
         if (n == 1) then
             array(1) = from
             return
         end if
  
         do i=1, n
             array(i) = from + range * (i - 1) / (n - 1)
         end do
      end subroutine
     
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! EVALUATION FUNCTIONS !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !DONE: To evaluate the Hoyle-Lyttleton expressions
      subroutine hl_and_energy(id, ierr)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         
         !subroutine variables
         integer :: k
         
         !getting pointer to star
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! print *, 'hl called'

         !mass is in g, radius is in cm, time is in s
         do k = 1, s% nz
            !hl parameters
            v(k) = SQRT(standard_cgrav*(M_ns + s%m(k))/s%R(k))   !cm/s
            Ra(k) = 2*standard_cgrav*M_ns/(v(k)**2)              !cm
            mdot_hl(k) = pi*(Ra(k)**2)*(s%rho(k))*v(k)           !g/s
            fd_hl(k) = mdot_hl(k)*v(k)                            !dyne
            edot_hl(k) = fd_hl(k)*v(k)                            !erg/s
            
            !energy parameters
            Eorb(k) = standard_cgrav*M_ns*(s%m(k))/(2*s%R(k))                                           !erg
            dEorb(k) = (standard_cgrav*M_ns/(2*s%R(k)))*((s%m(k)/s%R(k)) - 4*pi*(s%R(k)**2)*s%rho(k))   !erg/cm
         end do

      end subroutine hl_and_energy
      
      !DONE: To evaluate the MR15 expressions
      subroutine mr15(id, ierr)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         
         !subroutine variables
         integer :: j
         real(dp) :: f1, f2, f3, mu1, mu2, mu3, mu4
         
         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! print *, 'mr15 called'

         f1 = 1.91791946
         f2 = -1.52814698
         f3 = 0.75992092

         mu1 = -2.14034214
         mu2 = 1.94694764
         mu3 = 1.19007536
         mu4 = 1.05762477

         do j = 1, s%nz
            eps_rho(j) = Ra(j)/s%scale_height(j)
            fd_mr15_ratio(j) = f1 + f2*eps_rho(j) + f3*(eps_rho(j)**2)
            mdot_mr15_ratio(j) = (10)**(mu1 + mu2/(1 + mu3*eps_rho(j) + mu4*(eps_rho(j)**2)))
            mdot_mr15(j) = mdot_mr15_ratio(j)*mdot_hl(j)
            fd_mr15(j) = fd_mr15_ratio(j)*fd_hl(j)
         end do
      end subroutine mr15
      
      !DONE: To evaluate the Eddington and Hypereddington accretion rates
      subroutine edd_and_hyper(n)         
         !subroutine variables
         integer :: k, n

         ! print *, 'edd and hyper called'

         do k = 1, n
            mdot_edd(k) = 3.5*(1d-8)*(M_ns/(1.33*Msun))*(0.34/op_const)*Msun/secyer
            mdot_hyper(k) = 8.9*(1d-5)*((op_const/0.34)**(-0.73))*Msun/secyer
         end do
      end subroutine edd_and_hyper
      
      !DONE: To evaluate the Holgado accretion parametrisation
      subroutine hol_acc_par(n)     
         !subroutine variables
         integer :: k, n

         ! print *, 'holgado accretion parametrisation called'
         
         do k = 1, n
            if (eta*mdot_mr15(k) >= mdot_hyper(k) .or. eta*mdot_mr15(k) <= mdot_edd(k)) then
               mdot(k) = eta*mdot_mr15(k)
            else
               mdot(k) = mdot_edd(k)
            end if
            fd(k) = eta*fd_mr15(k)
         end do   
         
      end subroutine hol_acc_par

      !DONE: to evaluate the function F for the differential equation da/dt = F(a)
      subroutine orbital_evolution_function(id, ierr, t, a, fa) 
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         
         !subroutine variables
         real(dp), intent(in) :: a, t
         real(dp), intent(out) :: fa
         integer :: k, azone_temp
         
         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! print *, 'orbital evolution function called'

         !finding the zone where the function is to be evaluated
         call find_closest_value_index(s%R(1:s%nz), a, s%nz, azone_temp)

         !evaluating the function
         fa = f(azone_temp)

      end subroutine orbital_evolution_function

      
      !DONE: fourth order runge kutta solver
      !currently, this considers the orbital_evolution_function as the function to be solved
      subroutine rk4(id, ierr, a, t, dt, aout)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         
         !subroutine variables
         real(dp), intent(in) :: a, t, dt
         real(dp), intent(out) :: aout
         real(dp) :: k1, k2, k3, k4
         procedure(orbital_evolution_function), pointer :: G
         
         !getting pointer to star
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! print *, 'rk4 called'

         !setting the function pointer
         G => orbital_evolution_function
         
         !evaluating the k values
         call G(id, ierr, t, a, k1)
         call G(id, ierr, t + dt/2, a + k1/2, k2)
         call G(id, ierr, t + dt/2, a + k2/2, k3)
         call G(id, ierr, t + dt, a + k3, k4)
         
         !evaluating the next value of a
         aout = a + dt*(k1 + 2*k2 + 2*k3 + k4)/6
      end subroutine rk4
      
      !DONE: To find the equatorial radius and beta secular
      subroutine Req_and_beta(e1, Req1, Rbar1, beta1)
         real(dp), intent(in) :: e1
         real(dp), intent(out) :: Req1, Rbar1, beta1
         
         beta1 = 3*(1 - ((e1*sqrt(1-e1**2))/(asin(e1))))/(2*e1**2) - 1
         Rbar1 = R0*((asin(e1) * ((1-e1**2)**(1/6)) * (1-beta1)) / e1)**(-1 * n_poly / (3-n_poly))
         Req1 = Rbar1/((1-e1**2)**(1/6))
      end subroutine Req_and_beta

      ! subroutine omega_and_beta_vs_e(M, e_in, omega_out, beta_out)
      !    real(dp), intent(in) :: M
      !    real(dp), intent(in), allocatable :: e_in(:)
      !    real(dp), intent(out), allocatable :: omega_out(:), beta_out(:)
      !    real(dp) :: rho_bar, qn
      !    integer :: j

      !    rho_bar = 3*M/(4*pi*(R0**3))
      !    qn = (1-n_poly/5)

      !    allocate(omega_out(size(e_in)), beta_out(size(e_in)))

      !    do j = 1, size(e_in)
      !       omega_out(j) = sqrt(2*pi*standard_cgrav*rho_bar*( (sqrt(1-e_in(j)**2)*(3-2*e_in(j)**2)*asin(e_in(j))/(e_in(j)**3)) - 3*(1-e_in(j)**2)/(e_in(j)**2) )/qn)
      !       beta_out(j) = 3*(1 - ((e_in(j)*sqrt(1-e_in(j)**2))/(asin(e_in(j)))))/(2*e_in(j)**2) - 1         
      !    end do      
      ! end subroutine omega_and_beta_vs_e
      
      !DONE: To evaluate the omega function and solve for the inverse
      subroutine omega_function(e1, M, omega1)
         !subroutine variables
         real(dp), intent(in) :: e1, M
         real(dp), intent(out) :: omega1
         real(dp) :: rho_bar, qn

         rho_bar = 3*M/(4*pi*(R0**3))
         qn = (1-n_poly/5)
         
         omega1 = sqrt(2*pi*standard_cgrav*rho_bar*( (sqrt(1-e1**2)*(3-2*e1**2)*asin(e1)/(e1**3)) - 3*(1-e1**2)/(e1**2) )/qn)
      end subroutine omega_function
      
      !DONE: To solve for the inverse of the omega function
      subroutine omega_func_solve_inverse(e_start, omega_val, tol, M, e_out)
         !subroutine variables
         real(dp), intent(in) :: e_start, omega_val, tol, M
         real(dp), intent(out) :: e_out

         !local variables
         real(dp) :: e1, omega1
         integer :: iterations

         iterations = 0

         e1 = e_start
         call omega_function(e1, M, omega1)

         do while (abs(omega1 - omega_val) >= tol .and. iterations < 100)
            e1 = e1*omega_val/omega1
            call omega_function(e1, M, omega1)
            iterations = iterations + 1
         end do

         e_out = e1
      end subroutine omega_func_solve_inverse
      
      !DONE: To evaluate the beta function
      subroutine beta_func(e1, beta1)
         !subroutine variables
         real(dp), intent(in) :: e1
         real(dp), intent(out) :: beta1
         
         beta1 = 3*(1 - ((e1*sqrt(1-e1**2))/(asin(e1))))/(2*e1**2) - 1
      end subroutine beta_func

      !DONE: to evaluate the spin evolution and quadrupole moment evolution
      subroutine omega_and_q(id, ierr, M_next, M_acc_next, a_next_zone, Qmax_next, Qtb_next, Q_next, omega_next, Req_next, Rbar_next, beta_next, e_next, aa_next, b_next, mom_inert_next)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables
         real(dp), intent(in) :: M_next, M_acc_next
         integer, intent(in) :: a_next_zone
         real(dp), intent(out) :: Qmax_next, Qtb_next, Q_next, omega_next, Req_next, Rbar_next, beta_next, e_next
         real(dp), intent(out) :: aa_next, b_next, mom_inert_next
         real(dp) :: omega_func
         integer :: i, omega_next_ind

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         Qmax_next = 4*1d39*(M_acc_next/M_crust)
                  
         omega_func = ((mdot(azone)*sqrt(standard_cgrav*M_ns*Req)) - ((32*standard_cgrav*(omega**5)*(Q**2))/(5*(clight**5))))/mom_inert
         
         if (s%x_ctrl(18) /= s%x_ctrl(19)) then
            omega_next = omega + omega_func*s%dt_next
         else 
            omega_next = omega
         end if

         call omega_func_solve_inverse(e, omega_next, 1d-9, M_next, e_next)

         ! call find_closest_value_index(omega_arr, omega_next, size(omega_arr), omega_next_ind)
         ! e_next = e_arr(omega_next_ind)

         aa_next = R0*(1 + e_next/2)
         b_next = R0*(1 - e_next/2)
         mom_inert_next = M_next*(aa_next**2 + b_next**2)/5

         call Req_and_beta(e_next, Req_next, Rbar_next, beta_next)

         Qtb_next = sqrt((5*(clight**5)*mdot(a_next_zone)*sqrt(standard_cgrav*M_next*Req_next))/(32*standard_cgrav*(omega_next**5)))
         if (Qtb_next > Qmax_next) then
            Q_next = Qmax_next*exp(-1*s%time*mdot(a_next_zone)/(M_crust))
         else
            Q_next = Qtb_next*exp(-1*s%time*mdot(a_next_zone)/(M_crust))
         end if

      end subroutine omega_and_q
      
      !DONE: To do the orbital evolution and inject energy into the star
      subroutine inject_energy(id, ierr)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables
         real(dp) :: M_acc_next, M_ns_next, a_curr, a_next, r1, r2, efactor, temp
         real(dp) :: Qmax_next, Qtb_next, Q_next, omega_next, Req_next, Rbar_next, beta_next, e_next, h_next
         real(dp) :: aa_next, b_next, mom_inert_next
         ! real(dp), allocatable :: rand(:)
         ! real(dp) :: 
         integer, allocatable :: rzones(:)
         integer :: i, j, nr, a_next_zone! omega_next_ind, beta_sec_ind

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! INITIALISING AND ALLOCATING
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         s%x_ctrl(18) = s%model_number            !current model number for checking retries

         ! print *, 'inject energy called'

         !Initialising
         op_const = s%x_ctrl(6)                   !opacity constant (cm^2/g)
         eta = s%x_ctrl(7)                        !efficiency factor
         efactor = s%x_ctrl(8)                    !multiplication factor for injected energy
         M_crust = s%x_ctrl(20)*Msun              !mass of the crust (gm)
         n_poly = s%x_ctrl(24)                    !polytropic index
         beta_sec = s%x_ctrl(25)                  !beta secular
         D = s%x_ctrl(50)*1d3*pc                  !distance to the source (cm)
         R0 = s%x_ctrl(22)                        !equatorial radius (cm)

         !allocating relevant variables
         if (allocated(fd)) then 
            deallocate(fd, mdot, f, edot)
            deallocate(mdot_hl, fd_hl, edot_hl, v, Ra, Eorb, dEorb)
            deallocate(mdot_mr15, fd_mr15, mdot_mr15_ratio, fd_mr15_ratio, eps_rho)
            deallocate(mdot_edd, mdot_hyper)!, omega_arr, e_arr, beta_arr)!, rand)
         end if
         allocate(fd(s%nz), mdot(s%nz), f(s%nz), edot(s%nz))
         allocate(mdot_hl(s% nz), fd_hl(s% nz), edot_hl(s% nz), v(s% nz), Ra(s% nz), Eorb(s%nz), dEorb(s%nz))
         allocate(mdot_mr15(s%nz), fd_mr15(s%nz), mdot_mr15_ratio(s%nz), fd_mr15_ratio(s%nz), eps_rho(s%nz))
         allocate(mdot_edd(s%nz), mdot_hyper(s%nz))

         print *, 'model = ', s%model_number
         !initialising the orbital seperation
         if (s%model_number == 1) then 
            s%x_ctrl(2) = s%R(1)                          !cm
            a_curr = s%x_ctrl(2)                          !cm

            s%x_ctrl(5) = 0.0D0                           !total injected energy

            M_ns = s%x_ctrl(1)*Msun                  !gm
            s%x_ctrl(12) = M_ns                      !gm
            
            s%x_ctrl(16) = 0                         !gm
            M_acc = s%x_ctrl(16)                     !gm
            
            omega = 2*pi*s%x_ctrl(21)                !spin frequency (Hz)
            e = s%x_ctrl(23)                         !ellipticity
            Qmax = 4*1d39*(M_acc/M_crust)            !gm cm^2/s^2
            
            s%x_ctrl(42) = R0*(1+e/2)                !a (cm)
            a = s%x_ctrl(42)                         !cm
            s%x_ctrl(44) = R0*(1-e/2)                !b (cm)
            b = s%x_ctrl(44)                         !cm
            s%x_ctrl(46) = M_ns*(a**2 + b**2)/5      !moment of inertia (gm cm^2)
            mom_inert = s%x_ctrl(46)                 !gm cm^2

            s%x_ctrl(19) = 2                         !to check for retries
         ! else if (s%model_number == 3351) then (Only used when restarting - specific to the model restarting from, so needs to be changed accordingly)
         !    print *, 'blaq'
         !    a_curr  =  6.7202675798626758d12
         !    M_ns  =  2.6481245677598265d33
         !    M_acc  =  3.5394397314246373d30
         !    Qmax  =  1.4240282282171865d38
         !    Qtb  =  4.0017361680840237d40
         !    Q  =  9.6660149013829673d37
         !    omega  =  1.6148748790971920d3
         !    Req  =  9.9978230993989506d5
         !    Rbar  =  9.9978230993989506d5
         !    beta  =  4.4764110268387114d-3
         !    e  =  1.8177567016667656d-1
         !    a  =  1.0908878350833382d6
         !    b  =  9.0911216491666168d5
         !    mom_inert  =  1.0679998647068851d45

         !    h_val = 2*standard_cgrav*(omega**2)*Q/(D*(clight**4))
            
         !    if (s%x_ctrl(19) /= 3351) then
         !       s%x_ctrl(19) = 3352                     !to check for retries
         !    end if
         else 
            !setting current values to the ones outputted from last timestep

            a_curr = s%x_ctrl(3)                     !cm
            M_ns = s%x_ctrl(13)                      !gm
            M_acc = s%x_ctrl(17)                     !gm
            Qmax = s%x_ctrl(27)                      !gm cm^2/s^2
            Qtb = s%x_ctrl(29)                       !gm cm^2/s^2
            Q = s%x_ctrl(31)                         !gm cm^2/s^2
            omega = s%x_ctrl(33)                     !spin frequency (Hz)
            Req = s%x_ctrl(35)                       !equatorial radius (cm)
            Rbar = s%x_ctrl(37)                      !mean radius (cm)
            beta = s%x_ctrl(39)                      !beta secular
            e = s%x_ctrl(41)                         !ellipticity
            a = s%x_ctrl(43)                         !a (cm)
            b = s%x_ctrl(45)                         !b (cm)
            mom_inert = s%x_ctrl(47)                 !moment of inertia (gm cm^2)
            h_val = s%x_ctrl(49)                     !h
         endif

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! EVALUATING CURRENT PARAMETERS
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !evaluating the HL, Eddington, Hypereddington and MR15 model parameters and the energy variables
         call hl_and_energy(id, ierr)
         call mr15(id, ierr)
         call edd_and_hyper(s%nz)

         !Performing the accretion rate parametrisation as specified in holgado et al. 2018
         call hol_acc_par(s%nz)

         !evaluating edot and the function G for the differential equation
         do j = 1, s%nz
            ! mdot(j) = mdot_edd(j)
            ! fd(j) = mdot_edd(j)*v(j)
            edot(j) = fd(j)*v(j)
            f(j) = -fd(j)*v(j)*(1/dEorb(j))
         end do

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! FINDING ZONES
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         !zone array - to find the zone closest to a_curr
         call find_closest_value_index(s%R(1:s%nz), a_curr, s%nz, azone)
         print *, 'azone at a_curr = ', azone

         !finding the radius range where energy is to be injected
         r1 = a_curr - Ra(azone)
         r2 = a_curr + Ra(azone)
         ! print *, 'Enenrgy is injected from r1 = ', r1/Rsun, ' to r2 = ', r2/Rsun, 'Rsun' 
         ! print *, 'delta(r) = ', 2*Ra(azone)/Rsun, 'Rsun'

         call find_inbetween_value_index(s%R(1:s%nz), r1, r2, s%nz, rzones)
         nr = size(rzones)
         
         !setting initial signal value
         if (s%model_number == 1) then
            call Req_and_beta(e, Req, Rbar, beta)
            Qtb = sqrt((5*(clight**5)*mdot(azone)*sqrt(standard_cgrav*M_ns*Req))/(32*standard_cgrav*(omega**5)))
            Q = Qtb

            s%x_ctrl(48) = 2*standard_cgrav*(omega**2)*Q/(D*(clight**4))
            h_val = s%x_ctrl(48)
         end if

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! INJECTING ENERGY
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !injecting energy
         do j = 1, nr
            i = rzones(j)
            s% extra_heat(i)%val = efactor*edot(azone)/((s%dm(i)) * nr)
         end do

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! EVALUATING THE STOPPING CONDITION
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !calculating the binding energy and the change in orbital energy
         ebind = 0
         do i = 1, azone
            temp = (s%energy(i) - standard_cgrav*s%m(i)/s%R(i))*s%dm(i)
            ebind = ebind + temp
         end do
         eorb_change = Eorb(azone) - Eorb(1)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! FINDING VARIABLES FOR NEXT TIMESTEP
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !finding the next orbital seperation
         call rk4(id, ierr, a_curr, s%time, s%dt_next, a_next)
         call find_closest_value_index(s%R(1:s%nz), a_next, s%nz, a_next_zone)
         
         print *, 'x_ctrl(18)', s%x_ctrl(18)
         print *, 'x_ctrl(19)', s%x_ctrl(19)

         if (s%x_ctrl(18) /= s%x_ctrl(19)) then
            !finding the next neutron star mass
            M_ns_next = M_ns + mdot(azone)*s%dt_next

            !finding the next total accreted mass
            M_acc_next = M_acc + mdot(azone)*s%dt_next

         else !this is satisfied when there is a retry
            !finding the next neutron star mass
            M_ns_next = M_ns

            !finding the next total accreted mass
            M_acc_next = M_acc

         end if
     
         !finding the next quadrupole moment and spin frequency
         call omega_and_q(id, ierr, M_ns_next, M_acc_next, a_next_zone, Qmax_next, Qtb_next, Q_next, omega_next, Req_next, Rbar_next, beta_next, e_next, aa_next, b_next, mom_inert_next)
         
         !finding the next signal strength
         h_next = 2*standard_cgrav*(omega_next**2)*Q_next/(D*(clight**4))
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! SAVING RELEVANT VARIABLES
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !save other relevant variables
         s%x_ctrl(2) = a_curr                     !a_curr (cm)
         s%x_ctrl(3) = a_next                     !a_next (cm)

         s%x_ctrl(4) = edot(azone)*s%dt_next      !injected energy per timestep (ergs)
         s%x_ctrl(5) = s%x_ctrl(5) + s%x_ctrl(4)  !total injected energy (ergs)

         s%x_ctrl(9) = azone                      !azone
         s%x_ctrl(10) = r1                        !r1 (cm)
         s%x_ctrl(11) = r2                        !r2 (cm)

         s%x_ctrl(12) = M_ns                      !M_ns current (gm)
         s%x_ctrl(13) = M_ns_next                 !M_ns next (gm)

         s%x_ctrl(16) = M_acc                     !M_acc (gm)
         s%x_ctrl(17) = M_acc_next                !M_acc next (gm)

         s%x_ctrl(14) = -ebind                    !-binding energy (ergs)
         s%x_ctrl(15) = eorb_change               !change in orbital energy (ergs)
         
         s%x_ctrl(19) = s%model_number            !current model number

         s%x_ctrl(26) = Qmax                      !Qmax gm cm^2/s^2
         s%x_ctrl(27) = Qmax_next                 !Qmax next gm cm^2/s^2

         s%x_ctrl(28) = Qtb                       !Qtb gm cm^2/s^2
         s%x_ctrl(29) = Qtb_next                  !Qtb next gm cm^2/s^2

         s%x_ctrl(30) = Q                         !Q gm cm^2/s^2
         s%x_ctrl(31) = Q_next                    !Q next gm cm^2/s^2

         s%x_ctrl(32) = omega                      !omega (Hz)
         s%x_ctrl(33) = omega_next                 !omega next (Hz)

         s%x_ctrl(34) = Req                        !Req (cm)
         s%x_ctrl(35) = Req_next                   !Req next (cm)

         s%x_ctrl(36) = Rbar                       !Rbar (cm)
         s%x_ctrl(37) = Rbar_next                  !Rbar next (cm)

         s%x_ctrl(38) = beta                       !beta
         s%x_ctrl(39) = beta_next                  !beta next

         s%x_ctrl(40) = e                          !e
         s%x_ctrl(41) = e_next                     !e next

         s%x_ctrl(42) = a                          !a (cm)
         s%x_ctrl(43) = aa_next                     !a next (cm)

         s%x_ctrl(44) = b                          !b (cm)
         s%x_ctrl(45) = b_next                     !b next (cm)

         s%x_ctrl(46) = mom_inert                  !moment of inertia (gm cm^2)
         s%x_ctrl(47) = mom_inert_next             !next moment of inertia (gm cm^2)

         s%x_ctrl(48) = h_val                          !h
         s%x_ctrl(49) = h_next                     !h next
         
         print *, 'a_curr = ', a_curr/Rsun
         
         print *, 'Injected energy in dt = ', edot(azone)*s%dt_next

         print *, '-ebind = ', -ebind
         print *, 'eorb_change = ', eorb_change

         print *, 'Q_curr = ', Q

         print *, 'omega_curr = ', omega

         print *, 'e_curr = ', e

         print *, 'M_acc', M_acc

         print *, 'h = ', h_val

         print *, '##############################################################################'
         print *, 'ONE TIMESTEP DONE'
         print *, '##############################################################################'


      end subroutine inject_energy

      subroutine spin_evol(id, ierr)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


      end subroutine spin_evol
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.+
         ! otherwise we use a null_ version which does nothing (except warn).

         s% other_energy => inject_energy
         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

      end subroutine extras_controls
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_startup
      

      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0

      end function extras_start_step


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if
         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 22
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.

         !column 1: orbital seperation
         names(1) = 'a_curr'
         vals(1) = s%x_ctrl(2)

         !column 2: injected energy per timestep
         names(2) = 'Injected_E_per_timestep'
         vals(2) = s%x_ctrl(4)

         !column 3: total injected energy
         names(3) = 'Total_Injected_E'
         vals(3) = s%x_ctrl(5)

         !column 4: azone
         names(4) = 'azone'
         vals(4) = azone

         !column 5: r1
         names(5) = 'r1'
         vals(5) = s%x_ctrl(10)

         !column 6: r2
         names(6) = 'r2'
         vals(6) = s%x_ctrl(11)

         !column 7: M_ns
         names(7) = 'M_ns'
         vals(7) = s%x_ctrl(12)

         !column 8: M_acc
         names(8) = 'M_acc'
         vals(8) = s%x_ctrl(16)

         !column 9: -binding energy
         names(9) = '-ebind'
         vals(9) = s%x_ctrl(14)

         !column 10: change in orbital energy
         names(10) = 'delta_Eorb'
         vals(10) = s%x_ctrl(15)

         !column 11: Qmax
         names(11) = 'Qmax'
         vals(11) = s%x_ctrl(26)

         !column 12: Qtb
         names(12) = 'Qtb'
         vals(12) = s%x_ctrl(28)

         !column 13: Q
         names(13) = 'Q'
         vals(13) = s%x_ctrl(30)

         !column 14: omega
         names(14) = 'omega'
         vals(14) = s%x_ctrl(32)

         !column 15: Req
         names(15) = 'Req'
         vals(15) = s%x_ctrl(34)

         !column 16: Rbar
         names(16) = 'Rbar'
         vals(16) = s%x_ctrl(36)

         !column 17: beta
         names(17) = 'beta'
         vals(17) = s%x_ctrl(38)

         !column 18: e
         names(18) = 'e'
         vals(18) = s%x_ctrl(40)

         !column 19: a
         names(19) = 'a'
         vals(19) = s%x_ctrl(42)

         !column 20: b
         names(20) = 'b'
         vals(20) = s%x_ctrl(44)

         !column 21: moment of inertia
         names(21) = 'I'
         vals(21) = s%x_ctrl(46)

         !column 22: h
         names(22) = 'h'
         vals(22) = s%x_ctrl(48)

         ierr = 0
      
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 14
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'

         print *, 'data_for_extra_profile_columns called'

         names(1) = 'mdot_hl'
         do k = 1, s%nz
           vals(k,1) = mdot_hl(k)
         end do

         names(2) = 'mdot_mr15'
         do k = 1, s%nz
           vals(k,2) = mdot_mr15(k)
         end do

         names(3) = 'mdot'
         do k = 1, s%nz
           vals(k,3) = mdot(k)
         end do

         names(4) = 'fd_hl'
         do k = 1, s%nz
           vals(k,4) = fd_hl(k)
         end do

         names(5) = 'fd_mr15'
         do k = 1, s%nz
           vals(k,5) = fd_mr15(k)
         end do

         names(6) = 'fd'
         do k = 1, s%nz
           vals(k,6) = fd(k)
         end do

         names(7) = 'edot_hl'
         do k = 1, s%nz
           vals(k,7) = edot_hl(k)
         end do

         names(8) = 'edot'
         do k = 1, s%nz
           vals(k,8) = edot(k)
         end do

         names(9) = 'eps_rho'
         do k = 1, s%nz
           vals(k,9) = eps_rho(k)
         end do

         names(10) = 'v'
         do k = 1, s%nz
           vals(k,10) = v(k)
         end do

         names(11) = 'Ra'
         do k = 1, s%nz
           vals(k,11) = Ra(k)
         end do

         names(12) = 'Eorb'
         do k = 1, s%nz
           vals(k,12) = Eorb(k)
         end do

         names(13) = 'dEorb'
         do k = 1, s%nz
           vals(k,13) = dEorb(k)
         end do

         names(14) = 'f'
         do k = 1, s%nz
           vals(k,14) = f(k)
         end do

         ierr = 0
         
      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! OTHER FUNCTIONS !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !DONE: to evaluate the orbital velocity of the neutron star
      ! subroutine v_Ra(id, ierr, M_ns, v, Ra)
      !    !star variables
      !    integer, intent(in) :: id
      !    integer, intent(out) :: ierr
      !    type (star_info), pointer :: s
         
      !    !subroutine variables
      !    real(dp), intent(in) :: M_ns
      !    real(dp), intent(out), allocatable :: v(:), Ra(:)
      !    integer :: k
         
      !    !getting pointer to star
      !    ierr = 0
      !    call star_ptr(id, s, ierr)
      !    if (ierr /= 0) return
         
      !    allocate(v(s% nz), Ra(s%nz))

      !    !mass is in g, radius is in cm, time is in s

      !    do k = 1, s% nz
      !       v(k) = SQRT(standard_cgrav*(M_ns + s%m(k))/s%R(k))   !cm/s
      !       Ra(k) = 2*standard_cgrav*M_ns/(v(k)**2)
      !    end do

      ! end subroutine v_Ra

      
      !injecting energy in the specified zones
      ! subroutine other_energy_zones(id, ierr)
      !    use const_def, only: Rsun
      !    integer, intent(in) :: id
      !    integer, intent(out) :: ierr
      !    type (star_info), pointer :: s
         
      !    integer :: k, k1, kn, i, nz, rsize
      !    integer, allocatable :: zones(:), rzone(:)
      !    real(dp) :: const_e
         
      !    ierr = 0
      !    call star_ptr(id, s, ierr)
      !    if (ierr /= 0) return
         
      !    nz = s% nz                     !number of zones
      !    allocate(zones(nz))
      !    zones = [(i, i=1, nz, 1)]      !array of each zone from 1 to nz
         
      !    k1 = 500
      !    kn = 550
      !    const_e = 10**(9)
         
      !    do k = k1, kn
      !      s% extra_heat(k)%val = const_e
      !    end do
         
      !    print *, 'Adding energy ', const_e, ' in zones', k1, ' to ', kn ! r range 10 to 50 solar radius'

      !    return
      !    ! note that extra_heat is type(auto_diff_real_star_order1) so includes partials.
      ! end subroutine other_energy_zones
      
      
      !injecting energy in the specified radius
      ! subroutine other_energy_radius(id, ierr)
      !    use const_def, only: Rsun
      !    integer, intent(in) :: id
      !    integer, intent(out) :: ierr
      !    type (star_info), pointer :: s
         
      !    integer :: i, j, nz, rsize
      !    integer, allocatable :: zones(:), rzones(:)
      !    real(dp) :: E_extra_tot                 !ergs
      !    real(dp) :: r1, r2                      !Rsun
         
      !    !calling star pointer
      !    ierr = 0
      !    call star_ptr(id, s, ierr)
      !    if (ierr /= 0) return
         
         
      !    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !    ! SPECIFYING THE TOTAL EXTRA *ENERGY* (NOT EXTRA_HEAT)
      !    ! AND THE RADIUS RANGE WHERE IT IS TO BE ADDED
      !    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
      !    E_extra_tot = 1.0D44                        !ergs
      !    !print *, E_extra_tot
      !    r1 = 60.0                                   !Rsun
      !    r2 = 100.0                                  !Rsun
         
         
      !    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !    ! FINDING ZONES FROM r1 TO r2
      !    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
      !    !defining an array integers from 1 to nz
      !    nz = s% nz                                  !number of zones
      !    allocate(zones(nz))
      !    zones = [(i, i=1, nz, 1)]
         
      !    !array of zones that are in 
      !    !the range or radius r1 to r2 solar radius
      !    call find_inbetween_value_index(s%R, r1, r2, s%nz, rzones)  
         
         
      !    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !    ! INJECTING EXTRA HEAT
      !    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
      !    !looping over rzones and injecting energy at those zones
      !    do j = 1, size(rzones)
      !      i = rzones(j)
      !      s% extra_heat(i)%val = E_extra_tot/( (s%dt) * (s%dm(i)) * (size(rzones)) )
      !    end do
       
      !    !print *, 'Adding energy ', E_extra_tot, ' between radius', r1, ' to ', rn 
      !    print *, 'Extra heat at r1 = ', s% extra_heat(rzones(1))%val
      !    print *, 'denom(1) = ', ( (s%dt) * (s%dm(rzones(1))) * (size(rzones)) )


      !    return
      !    ! note that extra_heat is type(auto_diff_real_star_order1) so includes partials.
      ! end subroutine other_energy_radius

      end module run_star_extras
      
