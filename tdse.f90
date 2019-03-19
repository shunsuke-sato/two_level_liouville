module global_variables
  implicit none

! mathematical constants
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! model paramters
  real(8) :: Egap

  complex(8) :: zH_mat(2,2)

  complex(8) :: zpsi(2)

! parameters for time-propagation
  integer :: nt, nt_cycle
  real(8) :: dt, Tprop
  real(8),allocatable :: tt(:)
  

! laser paraemter
  real(8) :: E0, omega0, T0
  real(8),allocatable :: Et(:)
  

! Floquet
  integer,parameter :: nmax_floquet = 32
  complex(8),allocatable :: zpsi_F(:,:)
  real(8) :: eps_F(2)

! Photoelectron spectroscopy
  logical,parameter :: if_calc_PES = .true.
  integer :: it_PES_ini, it_PES_fin
  integer :: NE_PES
  real(8),allocatable :: eps_PES(:)
  real(8),allocatable :: pop_PES(:)
  complex(8),allocatable :: zdip_PES(:)
  real(8),allocatable :: Et_env_PES(:)
  real(8) :: omega_PES, omega_range_PES
  real(8) :: Tpulse_PES

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

  omega0 = Egap
  E0 = 0.1d0
  T0 = 10d0*2d0*pi/omega0


  Tprop = 60d0*2d0*pi/omega0
  dt = 0.005d0


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


  zpsi = 0d0
  zpsi(2) = 1d0

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

  allocate(Et(0:nt+1))

  do it = 0, nt+1
    xx = tt(it)
    if(xx < T0)then
      Et(it) = E0*sin(omega0*xx)*sin(0.5d0*pi*xx/T0)**2
    else
      Et(it) = E0*sin(omega0*xx)
    end if

  end do
  

end subroutine initialize_laser
!-------------------------------------------------------------------------------
subroutine time_propagation
  use global_variables
  implicit none
  integer :: it
  real(8) :: pop, dip


  open(21,file='quantities_t.out')
  write(21,"(999e26.16e3)")tt(0),Et(0),abs(zpsi(1))**2&
                                      ,abs(zpsi(2))**2&
                                      ,2d0*real(conjg(zpsi(1))*zpsi(2))

  it = 0

  do it = 0,nt

    call dt_evolve(it)

    write(21,"(999e26.16e3)")tt(it+1),Et(it+1),abs(zpsi(1))**2&
                                              ,abs(zpsi(2))**2&
                                              ,2d0*real(conjg(zpsi(1))*zpsi(2))

  end do

  close(21)

end subroutine time_propagation
!-------------------------------------------------------------------------------
subroutine dt_evolve(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  complex(8) :: zpsi_t(2)
  complex(8) :: zHam_mat(2,2)
  complex(8) :: zeigv(2,2)
  real(8) :: lambda(2)
  real(8) :: Et_tmp


! t -> t+dt/2

  Et_tmp = Et(it)

  zHam_mat(1,1) =  0.5d0*Egap
  zHam_mat(2,1) =  Et_tmp
  zHam_mat(1,2) =  Et_tmp
  zHam_mat(2,2) = -0.5d0*Egap

  call diag_2x2(zHam_mat,zeigv,lambda)

  zpsi_t(1) = sum(conjg(zeigv(:,1))*zpsi(:))
  zpsi_t(2) = sum(conjg(zeigv(:,2))*zpsi(:))

  zpsi(:) = zpsi_t(1)*exp(zI*lambda(1)*0.5d0*dt)*zeigv(:,1) &
           +zpsi_t(2)*exp(zI*lambda(2)*0.5d0*dt)*zeigv(:,2)


! t+dt/2 -> t+dt



  Et_tmp = Et(it+1)
  zHam_mat(1,1) =  0.5d0*Egap
  zHam_mat(2,1) =  Et_tmp
  zHam_mat(1,2) =  Et_tmp
  zHam_mat(2,2) = -0.5d0*Egap

  call diag_2x2(zHam_mat,zeigv,lambda)

  zpsi_t(1) = sum(conjg(zeigv(:,1))*zpsi(:))
  zpsi_t(2) = sum(conjg(zeigv(:,2))*zpsi(:))

  zpsi(:) = zpsi_t(1)*exp(zI*lambda(1)*0.5d0*dt)*zeigv(:,1) &
           +zpsi_t(2)*exp(zI*lambda(2)*0.5d0*dt)*zeigv(:,2)




end subroutine dt_evolve
!-------------------------------------------------------------------------------
subroutine diag_2x2(zmat,zvec,lambda)
  implicit none
  complex(8),intent(in)  :: zmat(2,2)
  complex(8),intent(out) :: zvec(2,2)
  real(8),intent(out)     :: lambda(2)
  real(8) :: a,b,ss
  complex(8) :: zc
  complex(8) :: zx,zy
  
  a = real(zmat(1,1))
  b = real(zmat(2,2))
  zc = zmat(2,1)

  if(a>b)then
     lambda(1) = 0.5d0*( (a+b)+sqrt((a-b)**2+4d0*abs(zc)**2) )
     lambda(2) = 0.5d0*( (a+b)-sqrt((a-b)**2+4d0*abs(zc)**2) )
  else
     lambda(1) = 0.5d0*( (a+b)-sqrt((a-b)**2+4d0*abs(zc)**2) )
     lambda(2) = 0.5d0*( (a+b)+sqrt((a-b)**2+4d0*abs(zc)**2) )
  end if
  zy = zc/(lambda(1)-b)
  zx = conjg(zc)/(lambda(2)-a)

! vector 1
  ss = 1d0/sqrt((1d0+abs(zy)**2))
  zvec(1,1) = ss
  zvec(2,1) = ss*zy
  
! vector 2
  ss = 1d0/sqrt((1d0+abs(zx)**2))
  zvec(1,2) = ss*zx
  zvec(2,2) = ss

  
end subroutine diag_2x2
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
