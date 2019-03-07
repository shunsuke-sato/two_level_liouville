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
  

! Floquet
  integer,parameter :: nmax_floquet = 8
  complex(8),allocatable :: zpsi_F(:,:)
  real(8) :: eps_F(2)

! Photoelectron spectroscopy
  logical,parameter :: if_calc_PES = .false.
  integer :: it_PES_ini, it_PES_fin
  integer :: NE_PES
  real(8),allocatable :: eps_PES(:)
  real(8),allocatable :: pop_PES(:)
  complex(8),allocatable :: zdip_PES(:)

end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none

  call input
  call initialize

  call calc_floquet_state
  call time_propagation



end program main
!-------------------------------------------------------------------------------
subroutine input
  use global_variables
  implicit none


  Egap = 1d0
  T1_relax = 30d0
  T2_relax = 30d0

  omega0 = Egap
  E0 = 1d-6


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
  real(8) :: s_floquet

  open(21,file='quantities_t.out')
  write(21,"(999e26.16e3)")tt(0),Et(0),real(zrho_dm(1,1)),real(zrho_dm(2,2)),&
                          &2d0*real(zi*zrho_dm(1,2))

  it = 0
  open(22,file='fidelity_t.out')
  call calc_floquet_fidelity_inst(it,s_floquet)
  write(22,"(999e26.16e3)")tt(0),s_floquet


  do it = 0,nt

    call dt_evolve(it)
    write(21,"(999e26.16e3)")tt(it+1),Et(it+1), &
                           & real(zrho_dm(1,1)),real(zrho_dm(2,2)),&
                           & 2d0*real(zrho_dm(1,2))

  call calc_floquet_fidelity_inst(it,s_floquet)
  write(22,"(999e26.16e3)")tt(it),s_floquet

  end do

  close(21)
  close(22)

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
  zHam_mat(2,1) =  Et_tmp
  zHam_mat(1,2) =  Et_tmp
  zHam_mat(2,2) = -0.5d0*Egap

!RK 1, t
  zrho_t = zrho_dm
  call Lrho_op(zrho_t,zHam_mat, zLrho_RK(:,:,1))


  Et_tmp = Et_dt2(it)
  zHam_mat(1,1) =  0.5d0*Egap
  zHam_mat(2,1) =  Et_tmp
  zHam_mat(1,2) =  Et_tmp
  zHam_mat(2,2) = -0.5d0*Egap

!RK 2, t+dt/2
  zrho_t = zrho_dm + dt*0.5d0*zLrho_RK(:,:,1)
  call Lrho_op(zrho_t,zHam_mat, zLrho_RK(:,:,2))

!RK 3, t+dt/2
  zrho_t = zrho_dm + dt*0.5d0*zLrho_RK(:,:,2)
  call Lrho_op(zrho_t,zHam_mat, zLrho_RK(:,:,3))


  Et_tmp = Et(it+1)
  zHam_mat(1,1) =  0.5d0*Egap
  zHam_mat(2,1) =  Et_tmp
  zHam_mat(1,2) =  Et_tmp
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
subroutine calc_floquet_state
  use global_variables
  implicit none
  integer :: nmat
  complex(8),allocatable :: zHam_F(:,:)
  complex(8) :: zHam00(2,2),zHam01(2,2),zHam10(2,2)
  integer :: ifloquet, jfloquet
  integer :: i,j

! for LAPACK    =====
  integer :: lwork,info
  complex(8),allocatable :: work(:)
  real(8),allocatable    :: rwork(:),w(:)

  nmat = (2*nmax_floquet + 1)*2

  lwork = min(64, nmat)*max(1,2*nmat-1)
  if(lwork < 1024)lwork = 1024
  allocate(w(nmat))
  allocate(work(lwork), rwork(max(1, 3*nmat-2)))
! for LAPACK    =====


  nmat = (2*nmax_floquet + 1)*2
  allocate(zHam_F(nmat,nmat))

  zHam00 = 0d0
  zHam00(1,1) =  0.5d0*Egap
  zHam00(2,2) = -0.5d0*Egap

  zHam01 = 0d0
  zHam01(2,1) = -0.5d0*zI*E0
  zHam01(1,2) = -0.5d0*zI*E0

  zHam10 = 0d0
  zHam10(2,1) =  0.5d0*zI*E0
  zHam10(1,2) =  0.5d0*zI*E0


  zHam_F = 0d0
  do ifloquet = -nmax_floquet,nmax_floquet

! (m,m)
    jfloquet = ifloquet
    i = 2*(ifloquet+nmax_floquet)+1
    j = 2*(jfloquet+nmax_floquet)+1

    zHam_F(i:i+1,j:j+1) = zHam00(1:2,1:2)
    zHam_F(i,i) = zHam_F(i,i) - ifloquet*omega0
    zHam_F(i+1,i+1) = zHam_F(i+1,i+1) - ifloquet*omega0

! (m,m+1)
    jfloquet = ifloquet + 1
    if(jfloquet <= nmax_floquet)then
      i = 2*(ifloquet+nmax_floquet)+1
      j = 2*(jfloquet+nmax_floquet)+1

      zHam_F(i:i+1,j:j+1) = zHam01(1:2,1:2)
    end if

! (m,m+1)
    jfloquet = ifloquet - 1
    if(jfloquet >= -nmax_floquet)then
      i = 2*(ifloquet+nmax_floquet)+1
      j = 2*(jfloquet+nmax_floquet)+1

      zHam_F(i:i+1,j:j+1) = zHam10(1:2,1:2)
    end if

  end do

  call zheev('V', 'U', nmat, zHam_F, nmat, w, work, lwork, rwork, info)


  allocate(zpsi_F(nmat,2))
  i = 2*(0+nmax_floquet)+1
  zpsi_F(:,1:2) = zHam_F(:,i:i+1)
  eps_F(1:2) = w(i:i+1)


end subroutine calc_floquet_state
!-------------------------------------------------------------------------------
subroutine calc_instantaneous_floquet_state(it, zpsi_F_inst)
  use global_variables
  implicit none
  integer,intent(in) :: it
  complex(8),intent(out) :: zpsi_F_inst(2,2)
  integer :: ifloquet, i
  real(8) :: xx

  xx = tt(it)

  zpsi_F_inst = 0d0
  do ifloquet = -nmax_floquet, nmax_floquet
    i = 2*(ifloquet+nmax_floquet)+1

    zpsi_F_inst(:,:) = zpsi_F_inst(:,:) &
      + exp(-zI*ifloquet*omega0*xx)*zpsi_F(i:i+1,1:2)
  end do

  zpsi_F_inst(:,1) = exp(-zI*eps_F(1)*xx)*zpsi_F_inst(:,1)
  zpsi_F_inst(:,2) = exp(-zI*eps_F(2)*xx)*zpsi_F_inst(:,2)



end subroutine calc_instantaneous_floquet_state
!-------------------------------------------------------------------------------
subroutine calc_floquet_fidelity_inst(it,s_floquet)
  use global_variables
  implicit none
  integer,intent(in) :: it
  real(8),intent(out) :: s_floquet
  complex(8) :: zpsi_F_inst(2,2)
  complex(8) :: zpsi_NO(2,2)
  real(8) :: smat(2,2)
  real(8) :: lambda(2)

  call calc_instantaneous_floquet_state(it, zpsi_F_inst)
  call diag_2x2(zrho_dm,zpsi_NO,lambda)

  smat(1,1) = abs(sum(conjg(zpsi_NO(:,1))*zpsi_F_inst(:,1)))**2
  smat(2,1) = abs(sum(conjg(zpsi_NO(:,2))*zpsi_F_inst(:,1)))**2
  smat(1,2) = abs(sum(conjg(zpsi_NO(:,1))*zpsi_F_inst(:,2)))**2
  smat(2,2) = abs(sum(conjg(zpsi_NO(:,2))*zpsi_F_inst(:,2)))**2

  s_floquet = smat(1,1)*smat(2,2)-smat(1,2)*smat(2,1)
  s_floquet = abs(s_floquet)


end subroutine calc_floquet_fidelity_inst
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
subroutine initialize_PES
  use global_variables
  implicit none

end subroutine initialize_PES
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
