module global_variables
  implicit none

! mathematical constants
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! model paramters
  real(8) :: Egap

  complex(8) :: zrho_dm(2,2)
  complex(8) :: zH_mat(2,2)

! relaxation paramters
  real(8) :: T1_relax, T2_relax

! parameters for time-propagation
  integer :: nt, nt_cycle
  real(8) :: dt, Tprop
  real(8),allocatable :: tt(:)
  

! laser paraemter
  real(8) :: E0, omega0
  real(8),allocatable :: Et(:),Et_dt2(:)
  



end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none

  call input
  call initialize


  call time_propagation



end program main
!-------------------------------------------------------------------------------
subroutine input
  use global_variables
  implicit none


  Egap = 1d0
  T1_relax = 1d0
  T2_relax = 1d0

  omega0 = Egap
  E0 = 1d-2


  Tprop = 60d0*2d0*pi/omega0
  dt = 0.1d0




  write(*,"(A,2x,e26.16e3)")"input   dt=",dt
  nt_cycle = (2d0*pi/omega0)/dt + 1
  dt = (2d0*pi/omega0)/nt_cycle

  write(*,"(A,2x,e26.16e3)")"refined dt=",dt
  nt = Tprop/dt + 1



end subroutine input
!-------------------------------------------------------------------------------
subroutine initialize
  use global_variables
  implicit none
  integer :: it

  zrho_dm = 0d0
  zrho_dm(2,2) = 1d0


  allocate(tt(0:nt+1))
  do it = 0, nt+1
    tt(it) = it*dt
  end do


  call initialize_laser

end subroutine initialize
!-------------------------------------------------------------------------------
subroutine initialize_laser
  use global_variables
  implicit none
  integer :: it
  real(8) :: xx

  allocate(Et(0:nt+1),Et_dt2(0:nt+1))

  do it = 0, nt+1
    xx = tt(it)
    Et(it) = E0*sin(omega0*xx)

    xx = tt(it) + dt*0.5d0
    Et_dt2(it) = E0*sin(omega0*xx)


  end do
  

end subroutine initialize_laser
!-------------------------------------------------------------------------------
subroutine time_propagation
  use global_variables
  implicit none
  integer :: it
  real(8) :: pop, dip


  open(21,file='quantities_t.out')
  write(21,"(999e26.16e3)")tt(0),Et(0),real(zrho_dm(1,1)),real(zrho_dm(2,2)),&
                          &2d0*real(zi*zrho_dm(1,2))


  do it = 0,nt

    call dt_evolve(it)
    write(21,"(999e26.16e3)")tt(it+1),Et(it+1), &
                           & real(zrho_dm(1,1)),real(zrho_dm(2,2)),&
                           & 2d0*real(zi*zrho_dm(1,2))

  end do

  close(21)

end subroutine time_propagation
!-------------------------------------------------------------------------------
subroutine dt_evolve(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  complex(8) :: zrho_t(2,2)
  complex(8) :: zHam_mat(2,2)
  complex(8) :: zLrho_RK(2,2,4)
  real(8) :: Et_tmp


  Et_tmp = Et(it)

  zHam_mat(1,1) =  0.5d0*Egap
  zHam_mat(2,1) =  zI*Et_tmp
  zHam_mat(1,2) = -zI*Et_tmp
  zHam_mat(2,2) = -0.5d0*Egap

!RK 1, t
  zrho_t = zrho_dm
  call Lrho_op(zrho_t,zHam_mat, zLrho_RK(:,:,1))


  Et_tmp = Et_dt2(it)
  zHam_mat(1,1) =  0.5d0*Egap
  zHam_mat(2,1) =  zI*Et_tmp
  zHam_mat(1,2) = -zI*Et_tmp
  zHam_mat(2,2) = -0.5d0*Egap

!RK 2, t+dt/2
  zrho_t = zrho_dm + dt*0.5d0*zLrho_RK(:,:,1)
  call Lrho_op(zrho_t,zHam_mat, zLrho_RK(:,:,2))

!RK 3, t+dt/2
  zrho_t = zrho_dm + dt*0.5d0*zLrho_RK(:,:,2)
  call Lrho_op(zrho_t,zHam_mat, zLrho_RK(:,:,3))


  Et_tmp = Et(it+1)
  zHam_mat(1,1) =  0.5d0*Egap
  zHam_mat(2,1) =  zI*Et_tmp
  zHam_mat(1,2) = -zI*Et_tmp
  zHam_mat(2,2) = -0.5d0*Egap

!RK 4, t+dt
  zrho_t = zrho_dm + dt*zLrho_RK(:,:,3)
  call Lrho_op(zrho_t,zHam_mat, zLrho_RK(:,:,4))


  zrho_dm = zrho_dm + dt/6d0*(zLrho_RK(:,:,1) &
                         +2d0*zLrho_RK(:,:,2) &
                         +2d0*zLrho_RK(:,:,3) &
                         +    zLrho_RK(:,:,4))



end subroutine dt_evolve
!-------------------------------------------------------------------------------
subroutine Lrho_op(zrho_in, zHam_in, zLrho_out)
  use global_variables
  implicit none
  logical,parameter :: if_hermite = .true.
  complex(8),intent(in)  :: zrho_in(2,2), zHam_in(2,2)
  complex(8),intent(out) :: zLrho_out(2,2)


  zLrho_out = -zI*(matmul(zHam_in,zrho_in) - matmul(zrho_in,zHam_in))

  if(if_hermite)then
  
    zLrho_out(1,1) = real(zLrho_out(1,1) -zrho_in(1,1)/T1_relax)
    zLrho_out(2,1) = zLrho_out(2,1) -zrho_in(2,1)/T2_relax
    zLrho_out(1,2) = conjg(zLrho_out(2,1))
    zLrho_out(2,2) = real(zLrho_out(2,2) +zrho_in(1,1)/T1_relax)
  else
    zLrho_out(1,1) = zLrho_out(1,1) -zrho_in(1,1)/T1_relax
    zLrho_out(2,1) = zLrho_out(2,1) -zrho_in(2,1)/T2_relax
    zLrho_out(1,2) = zLrho_out(1,2) -zrho_in(1,2)/T2_relax
    zLrho_out(2,2) = zLrho_out(2,2) +zrho_in(1,1)/T1_relax
  end if

  
end subroutine Lrho_op
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
