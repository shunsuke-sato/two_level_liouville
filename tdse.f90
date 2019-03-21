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
  logical,parameter :: if_floquet_state = .true.
  integer,parameter :: nmax_floquet = 32
  complex(8),allocatable :: zpsi_F(:,:)
  real(8) :: eps_F(2)

! Photoelectron spectroscopy
  logical,parameter :: if_calc_PES = .true.
  integer :: it_PES_ini, it_PES_fin
  integer :: NE_PES
  real(8),allocatable :: eps_PES(:)
  complex(8),allocatable :: zCt_PES(:)
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

  if(if_calc_pes)call initialize_PES

  call time_propagation

end program main
!-------------------------------------------------------------------------------
subroutine input
  use global_variables
  implicit none


  Egap = 1d0

  omega0 = Egap
  E0 = 0.3d0
  T0 = 20d0*2d0*pi/omega0


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

  if(if_floquet_state) call calc_floquet_state
  call initialize_laser

end subroutine initialize
!-------------------------------------------------------------------------------
subroutine initialize_laser
  use global_variables
  implicit none
  integer :: it
  real(8) :: xx

  allocate(Et(0:nt+1))
  Et = 0d0

  if(if_floquet_state)then
    do it = 0, nt+1
      xx = tt(it)
      Et(it) = E0*sin(omega0*xx)
    end do

  else

    do it = 0, nt+1
      xx = tt(it)
      if(xx < T0)then
        Et(it) = E0*sin(omega0*xx)*sin(0.5d0*pi*xx/T0)**2
      else
        Et(it) = E0*sin(omega0*xx)
      end if
      
    end do
  end if

end subroutine initialize_laser
!-------------------------------------------------------------------------------
subroutine time_propagation
  use global_variables
  implicit none
  integer :: it, ipes
  real(8) :: pop, dip


  open(21,file='quantities_t.out')
  write(21,"(999e26.16e3)")tt(0),Et(0),abs(zpsi(1))**2&
                                      ,abs(zpsi(2))**2&
                                      ,2d0*real(conjg(zpsi(1))*zpsi(2))

  it = 0

  do it = 0,nt

    if(if_calc_pes .and. it>=it_PES_ini .and. it<=it_PES_fin)call dt_evolve_PES_1st_half(it)
    call dt_evolve(it)
    if(if_calc_pes .and. it>=it_PES_ini .and. it<=it_PES_fin)call dt_evolve_PES_2nd_half(it)

    write(21,"(999e26.16e3)")tt(it+1),Et(it+1),abs(zpsi(1))**2&
                                              ,abs(zpsi(2))**2&
                                              ,2d0*real(conjg(zpsi(1))*zpsi(2))

  end do

  close(21)

  if(if_calc_pes)then
    open(20,file='pes_spectrum.out')
    do ipes = 0,NE_PES
      write(20,"(999e26.16e3)")eps_PES(ipes)-omega_PES,abs(zCt_PES(ipes))**2
    end do
    close(20)
  end if

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

  zpsi(:) = zpsi_t(1)*exp(-zI*lambda(1)*0.5d0*dt)*zeigv(:,1) &
           +zpsi_t(2)*exp(-zI*lambda(2)*0.5d0*dt)*zeigv(:,2)


! t+dt/2 -> t+dt



  Et_tmp = Et(it+1)
  zHam_mat(1,1) =  0.5d0*Egap
  zHam_mat(2,1) =  Et_tmp
  zHam_mat(1,2) =  Et_tmp
  zHam_mat(2,2) = -0.5d0*Egap

  call diag_2x2(zHam_mat,zeigv,lambda)

  zpsi_t(1) = sum(conjg(zeigv(:,1))*zpsi(:))
  zpsi_t(2) = sum(conjg(zeigv(:,2))*zpsi(:))

  zpsi(:) = zpsi_t(1)*exp(-zI*lambda(1)*0.5d0*dt)*zeigv(:,1) &
           +zpsi_t(2)*exp(-zI*lambda(2)*0.5d0*dt)*zeigv(:,2)




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
subroutine initialize_PES
  use global_variables
  implicit none
  real(8) :: wi, wf, dw
  integer :: iw, it
  real(8) :: xx

  if(.not. if_calc_PES) stop 'Error in initialize_PES.'

  NE_PES = 512
  omega_PES = Egap * 100d0
  omega_range_PES = 3d0*Egap

  Tpulse_PES = 20d0*2d0*pi/omega0

  allocate(eps_PES(0:NE_PES), zCt_PES(0:NE_PES))
  allocate(Et_env_PES(0:nt+1))
  zCt_PES = 0d0

  wi = omega_PES - omega_range_PES
  wf = omega_PES + omega_range_PES
  dw = (wf-wi)/NE_PES

  do iw = 0, NE_PES
    
    eps_PES(iw) = wi + dw*iw
    
  end do

  it_PES_ini = -1
  it_PES_fin = -1


  Et_env_PES = 0d0
  do it = 0,nt+1
    if(tt(it) >= Tprop - Tpulse_PES)then
      if(it_PES_ini == -1)it_PES_ini = it
      xx = tt(it) - (Tprop - Tpulse_PES)
      Et_env_PES(it) = sin(pi*xx/tpulse_PES)**2
    else if(tt(it) >= Tprop)then
      if(it_PES_fin == -1)it_PES_fin = it
    end if
  end do
  if(it_PES_fin == -1)it_PES_fin = nt+1
  write(*,"(A,2x,I7)")"it_pes_ini=",it_pes_ini
  write(*,"(A,2x,I7)")"it_pes_fin=",it_pes_fin

  open(20,file='Et_env_PES.out')
  do it = 0,nt+1
    write(20,"(999e26.16e3)")tt(it),Et_env_PES(it)
  end do
  close(20)
  

end subroutine initialize_PES
!-------------------------------------------------------------------------------
subroutine dt_evolve_PES_1st_half(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer :: ipes
  complex(8) :: zs


  do ipes = 0, NE_PES
    zs = 0.5d0*dt*0.5d0*Et_env_PES(it)*zpsi(2) ! |g>

    zCt_PES(ipes) = zCt_PES(ipes) + zs*exp(zI*(eps_PES(ipes) - omega_PES)*tt(it))
  end do
  

end subroutine dt_evolve_PES_1st_half
!-------------------------------------------------------------------------------
subroutine dt_evolve_PES_2nd_half(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer :: ipes
  complex(8) :: zs


  do ipes = 0, NE_PES
    zs = 0.5d0*dt*0.5d0*Et_env_PES(it+1)*zpsi(2) ! |g>

    zCt_PES(ipes) = zCt_PES(ipes) + zs*exp(zI*(eps_PES(ipes) - omega_PES)*tt(it))
  end do
  
end subroutine dt_evolve_PES_2nd_half
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

  zpsi = 0d0
  do ifloquet = -nmax_floquet, nmax_floquet
    i = 2*(ifloquet+nmax_floquet)+1

    zpsi(:) = zpsi(:) &
      + zpsi_F(i:i+1,2) ! 1 or 2
  end do


end subroutine calc_floquet_state
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
