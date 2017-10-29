!=======================================================================
!   This program solves the Euler equations (2D) with the Lax Method
!=======================================================================
!un cambio diferente en el programa
!otro cambio
!mas cambios
!=======================================================================
!   This module contains global variables
      module globals
      implicit none

!------------------------------------------------------------------------------
!   nx = number of points used to discretize X
!   ny = number of points used to discretize Y
!   neq = number of equations

!   xmax = extent of X
!   ymax = extent of Y
!   gamma = 4/3 (1.3333) or 7/5 (1.4) or 5/3 (1.6666)
!
!   tmax = max time for the simulation
!   dtprint = outputs every dtprint
!
!   rhoin, pin, vxin, vyin = density, pressure, velocity at x<rad, y<rad
!   rhoout, pout, vxout, vyout = density, pressure, velocity at x>rad, y>rad
!
!   Co = courant number
!
!   **we can modify these values**
!------------------------------------------------------------------------------
      integer, parameter :: nx=10, ny=10, neq=4
      real, parameter :: xmax=1., dx=xmax/float(nx)
      real, parameter :: ymax=1., dy=ymax/float(ny)
      real, parameter :: gamma=1.6666666

      real, parameter :: tmax= 1.0        ! maximum integration time
      real, parameter :: dtprint=0.05       ! interval between outputs

      real, parameter :: rhoin = 1.0
      real, parameter :: rhoout = 1.0
      real, parameter :: pin = 1.0
      real, parameter :: pout = 1.0
      real, parameter :: vxin = 0.0
      real, parameter :: vyin = 0.0
      real, parameter :: vxout = 0.0
      real, parameter :: vyout = 0.0
      real, parameter :: xc = 0.5
      real, parameter :: yc = 0.5

      real, parameter :: Co = 0.7

!--------------------------------------
! Boundary Conditions
! 0.0 = open boundary conditions
! 1.0 = Jet injection
!--------------------------------------
      real, parameter :: bound=0.0

!------------------------------------------------------------------------------
!   This is a vector that contains u(x,y) (do not touch)
!------------------------------------------------------------------------------
      real :: u (neq,0:nx+1,0:ny+1)
      real :: qu(neq)
      real :: qp(neq)

        end module globals
!------------------------------------------------------------------------------
! End module globals
!==============================================================================

!--------

!==============================================================================
!   This module is the main program
!------------------------------------------------------------------------------
      program euler2D_lax
      use globals
      implicit none
      real      :: time, dt         ! t, $\Delta t$
      real      :: tprint           ! time of next output
      integer   :: itprint          ! number of current output

!------------------------------------------------------------------------------
! This subroutine generates the initial conditions main loop,
! iterate until maximum time is reached
!------------------------------------------------------------------------------
      call initconds(time, tprint, itprint)     ! conditions at t=0
      do while (time.lt.tmax)                   ! stops when t = tmax
        if(time.ge.tprint) then                 ! prints data in terminal
        !  write (*,*) itprint+1,time,tmax,dt
          call output (itprint)                 ! print data output
          tprint=tprint+dtprint
          itprint=itprint+1
        end if

!------------------------------------------------------------------------------
! Calculates and uses the CFL dt criterium
!------------------------------------------------------------------------------
        call courant(dt)
!------------------------------------------------------------------------------
! Calculates U(time) and F => U(time+dt)
!------------------------------------------------------------------------------
        call ulax(dt,time)
        time=time+dt

      end do

      stop
      end program euler2D_lax
!------------------------------------------------------------------------------
! end of main program module
!==============================================================================

!--------

!==============================================================================
! In this module we set the initial condition
!------------------------------------------------------------------------------
      subroutine initconds(time,tprint,itprint)
      use globals
      implicit none
      real, intent(out) :: time, tprint
      integer, intent (out) :: itprint
      integer ::i,j
      real :: x,y, rad, lorin,lorout, hin,hout

!------------------------------------------------------------------------------
! For the 2D circular blast:
! u(1,i,j) = rho(i,j)
! u(2,i,j) = vx(i,j)
! u(3,i,j) = vy(i,j)
! u(4,i,j) = etot(i,j) = eint + ekin = P/(gamma-1)
!------------------------------------------------------------------------------
      do i=0,nx+1
        do j=0,ny+1
          x=float(i)*dx          ! obtain the position $x_i$
          y=float(j)*dy          ! obtain the position $y_j$
          rad=sqrt((x-xc)**2+(y-yc)**2)
!lorin=1.d0
!hin=1.d0
!lorout=1.d0
!hout=1.d0
          if (rad < 0.0) then

            lorin=1.d0!(1-(vxin**2+vyin**2))**(-1/2)
            hin=1.d0!1.d0+gamma/(gamma-1.d0)*rhoin


            u(1,i,j)=rhoin*lorin
            u(2,i,j)=rhoin*vxin*lorin**2*hin
            u(3,i,j)=rhoin*vyin*lorin**2*hin
            u(4,i,j)=pin/(gamma-1.) + 0.5*u(2,i,j)*u(2,i,j)/u(1,i,j) + 0.5*u(3,i,j)*u(3,i,j)/u(1,i,j)
          else

            lorout=1.d0!(1-(vxout**2+vyout**2))**(-1/2)
            hout=1.d0!1.d0+gamma/(gamma-1.d0)*rhoout


            u(1,i,j)=rhoout*lorout
            u(2,i,j)=rhoout*vxout*lorout**2*hout
            u(3,i,j)=rhoout*vyout*lorout**2*hout
            u(4,i,j)=pout/(gamma-1.) + 0.5*u(2,i,j)*u(2,i,j)/u(1,i,j) + 0.5*u(3,i,j)*u(3,i,j)/u(1,i,j)

         end if

!write (*,*) u(1,1,1),p(1)

qu(1) = u(1,i,j)
qu(2) = u(2,i,j)
qu(3) = u(3,i,j)
qu(4) = u(4,i,j)
 print*, "entra" , rhoout , vxout , vyout , pout
call uprim(qu,qp)
         print*, "sale" , qp(1) , qp(2) , qp(3) , qp(4)          

!r(1) 
        end do
       


        !call uprim
        !write(*,*)
      end do
!write (*,*) u(1,i,j),u(2,i,j),u(3,i,j),u(4,i,j)
!------------------------------------------------------------------------------
! end of the 2D circular blast initial condition
! reset the counters and time to 0
!------------------------------------------------------------------------------
      time=0
      tprint=0
      itprint=0

      return
      end subroutine initconds
!------------------------------------------------------------------------------
! end of the init condition module
!==============================================================================

!--------

!==============================================================================
! output files module
!------------------------------------------------------------------------------
      subroutine output(itprint)
      use globals
      implicit none
      integer, intent(in) :: itprint
      character (len=20) filerho, fileE, fileVel
      integer :: i,j
      real :: rho,vx,vy,P

!------------------------------------------------------------------------------
! the name of the output files are:
!------------------------------------------------------------------------------
      write(filerho,'(a,i2.2,a)') '2D_lax',itprint,'.txt'
      open(unit=10,file=filerho,status='unknown')

!------------------------------------------------------------------------------
! writes the output datafiles (ascii)
! writes x, y, rho, vx, vy and P
!------------------------------------------------------------------------------
      do j=1,ny
        do i=1,nx
          rho=u(1,i,j)
          vx=u(2,i,j)/rho
          vy=u(3,i,j)/rho
          P=(u(4,i,j)-0.5*rho*(vx**2)*(vy**2))*(gamma-1.)
          write(10,'(7es12.5)') float(i)*dx,float(j)*dy,rho, vx, vy, P, u(4,i,j)
        end do
        write(10,*)
      end do
      close(10)

      return
      end subroutine output
!------------------------------------------------------------------------------
! end of the output files module
!==============================================================================

!--------

!==============================================================================
! CFL criterium module
!------------------------------------------------------------------------------
      subroutine courant(dt)
      use globals
      implicit none
      real, intent(out) ::dt
      real :: rho, vx, vy, P, cs
      integer :: i,j

!------------------------------------------------------------------------------
! Calculate the CFL criterium
!------------------------------------------------------------------------------
      dt=1E30
      do i=0,nx+1
        do j=0,ny+1
          rho=u(1,i,j)
          vx=u(2,i,j)/rho
          vy=u(3,i,j)/rho
          P=(u(4,i,j)-0.5*rho*(vx**2+vy**2))*(gamma-1.)
          cs=sqrt(gamma*P/rho)
          dt=min( dt,Co*dx/(abs(vx)+cs) )
          dt=min( dt,Co*dy/(abs(vy)+cs) )

        end do
      end do

      return
      end subroutine courant
!------------------------------------------------------------------------------
! end of the CFL criterium module
!==============================================================================

!--------

!==============================================================================
!Lax integration method module to calculate the conservative variables
!------------------------------------------------------------------------------
      subroutine ulax(dt,time)
      use globals
      implicit none
      real, intent(in) :: dt, time
      real :: up(neq,0:nx+1,0:ny+1), f(neq,0:nx+1,0:ny+1), g(neq,0:nx+1,0:ny+1)
      real :: dtx, dty
      integer :: i,j

!------------------------------------------------------------------------------
! obtain the fluxes
!------------------------------------------------------------------------------
      call fluxes(nx,ny,neq,gamma,u,f,g,bound)

!------------------------------------------------------------------------------
!   Here is the Lax method, notice that the values at the extremes can not be
!   calculated, we need to enter then as a boundary conditions
!------------------------------------------------------------------------------
      dtx=dt/dx
      dty=dt/dy

      do i=1,nx
        do j=1,ny
          up(:,i,j)=0.25*( u(:,i-1,j)+u(:,i+1,j)+u(:,i,j-1)+u(:,i,j+1) ) &
                         -dtx*0.5*(f(:,i+1,j)-f(:,i-1,j) ) &
                         -dty*0.5*(g(:,i,j+1)-g(:,i,j-1) )
        end do
      end do

!------------------------------------------------------------------------------
!   Boundary conditions to the U^n+1
!------------------------------------------------------------------------------
      call boundaries(nx,ny,neq,up,time,gamma,bound)
!------------------------------------------------------------------------------
!   copy the up to the u
!------------------------------------------------------------------------------
      u(:,:,:)=up(:,:,:)

      return
      end subroutine ulax
!------------------------------------------------------------------------------
!   end of the lax integration method module
!==============================================================================

!--------

!==============================================================================
! Calculation of the fluxes module
!------------------------------------------------------------------------------
      subroutine fluxes(nx,ny,neq,gamma,u,f,g,bound)
      implicit none
      integer, intent(in) :: nx,ny,neq
      real, intent(in) :: gamma
      real, intent(in) :: bound
      real, intent(in) :: u(neq,0:nx+1,0:ny+1)
      real, intent(out) :: f(neq,0:nx+1,0:ny+1), g(neq,0:nx+1,0:ny+1)
      integer :: i, j
      real :: rho, vx, vy, P,lor, h
	lor=1.d0
	h=1.d0

      do i=0,nx+1
        do j=0,ny+1
          rho=u(1,i,j)
          vx=u(2,i,j)/rho
          vy=u(3,i,j)/rho
          P=(u(4,i,j)-0.5*rho*(vx**2+vy**2))*(gamma-1.)

          f(1,i,j)=rho*vx*lor
          f(2,i,j)=rho*vx*vx*lor**2*h+P
          f(3,i,j)=rho*vx*vy*lor**2*h
          f(4,i,j)=vx*(u(4,i,j)+P)
	 ! f(4,i,j)=rho*vx*lor**2*h

	  g(1,i,j)=rho*vy*lor
          g(2,i,j)=rho*vx*vy*lor**2*h
          g(3,i,j)=rho*vy*vy*lor**2*h+P
          g(4,i,j)=vy*(u(4,i,j)+P)
	! g(4,i,j)=rho*vy*lor**2*h


        end do
      end do

      return
      end subroutine fluxes
!------------------------------------------------------------------------------
! end of the fluxes module
!==============================================================================

!--------

!==============================================================================
! Boundary conditions module
!------------------------------------------------------------------------------
      subroutine boundaries(nx,ny,neq,u,time,gamma,bound)
      implicit none
      integer, intent(in) :: nx,ny,neq
      real, intent(in) :: bound
      real, intent(out) ::u(neq,0:nx+1,0:ny+1)
      real, intent(in) :: time,gamma
      integer :: j,i
      real :: rhoj , vj, Pj, Tauj , rhoi, vi, Pi, Taui, theta

!------------------------------------------------------------------------------
! Periodic or open boundaries, or your own boundary condition
!------------------------------------------------------------------------------

        u(:,0   ,:)=u(:,1 ,:)
        u(:,nx+1,:)=u(:,nx,:)
        u(:,:,0   )=u(:,:,1 )
        u(:,:,ny+1)=u(:,:,ny)

!---------------------------------------

        if(bound.eq.0.)then
          u(:,0   ,:)=u(:,1 ,:)
          u(:,nx+1,:)=u(:,nx,:)
          u(:,:,0   )=u(:,:,1 )
          u(:,:,ny+1)=u(:,:,ny)
        end if

!---------------------------------------

        if(bound.eq.1.)then
           u(1,0   ,:)=u(1,1 ,:)
           u(2,0   ,:)=-u(2,1 ,:)  !-u(2,1 ,:) si se quiere ref en x=0
           u(3,0   ,:)=u(3,1 ,:)
           u(4,0   ,:)=u(4,1 ,:)
         end if

!---------------------------------------

        if(bound.eq.2.)then
          u(1,nx+1,:)=u(1,nx,:)
          u(2,nx+1,:)=-u(2,nx,:)  !-u(2,nx,:) si se quiere ref en x=nx
          u(3,nx+1,:)=u(3,nx,:)
          u(4,nx+1,:)=u(4,nx,:)
         end if

!---------------------------------------

        if(bound.eq.3.)then
          u(1,:,0   )=u(1,:,1 )
          u(2,:,0   )=u(2,:,1 )
          u(3,:,0   )=-u(3,:,1 )  !-u(3,:,1 ) si se quiere ref en y=0
          u(4,:,0   )=u(4,:,1 )
         end if

!---------------------------------------

        if(bound.eq.4.)then
          u(1,:,ny+1)=u(1,:,ny)
          u(2,:,ny+1)=u(2,:,ny)
          u(3,:,ny+1)=-u(3,:,ny)  !-u(3,:,ny) si se quiere ref en y=ny
          u(4,:,ny+1)=u(4,:,ny)
         end if

!---------------------------------------

         if(bound.eq.5.)then
           u(1,0   ,:)=u(1,1 ,:)
           u(2,0   ,:)=-u(2,1 ,:)  !-u(2,1 ,:) si se quiere ref en x=0
           u(3,0   ,:)=u(3,1 ,:)
           u(4,0   ,:)=u(4,1 ,:)

          u(1,nx+1,:)=u(1,nx,:)
          u(2,nx+1,:)=-u(2,nx,:)  !-u(2,nx,:) si se quiere ref en x=nx
          u(3,nx+1,:)=u(3,nx,:)
          u(4,nx+1,:)=u(4,nx,:)
         end if

!---------------------------------------

        if(bound.eq.6.)then
          u(1,:,0   )=u(1,:,1 )
          u(2,:,0   )=u(2,:,1 )
          u(3,:,0   )=-u(3,:,1 )  !-u(3,:,1 ) si se quiere ref en y=0
          u(4,:,0   )=u(4,:,1 )

          u(1,:,ny+1)=u(1,:,ny)
          u(2,:,ny+1)=u(2,:,ny)
          u(3,:,ny+1)=-u(3,:,ny)  !-u(3,:,ny) si se quiere ref en y=ny
          u(4,:,ny+1)=u(4,:,ny)
         end if

!---------------------------------------

         if(bound.eq.7.)then
          u(1,0   ,:)=u(1,1 ,:)
          u(2,0   ,:)=-u(2,1,:)  !-u(2,1 ,:) si se quiere ref en x=0
          u(3,0   ,:)=u(3,1 ,:)
          u(4,0   ,:)=u(4,1 ,:)

          u(1,nx+1,:)=u(1,nx,:)
          u(2,nx+1,:)=-u(2,nx,:)  !-u(2,nx,:) si se quiere ref en x=nx
          u(3,nx+1,:)=u(3,nx,:)
          u(4,nx+1,:)=u(4,nx,:)

          u(1,:,0   )=u(1,:,1 )
          u(2,:,0   )=u(2,:,1 )
          u(3,:,0   )=-u(3,:,1 )  !-u(3,:,1 ) si se quiere ref en y=0
          u(4,:,0   )=u(4,:,1 )

          u(1,:,ny+1)=u(1,:,ny)
          u(2,:,ny+1)=u(2,:,ny)
          u(3,:,ny+1)=-u(3,:,ny)  !-u(3,:,ny) si se quiere ref en y=ny
          u(4,:,ny+1)=u(4,:,ny)
         end if

!---------------------------------------
!  Jet injection
!---------------------------------------
         if(bound.eq.10.)then
           rhoj=9.5
           Tauj=100.2
           vj=0.5 !+ 0.3*sin(2*3.1416*time/Tauj)
           Pj=0.01     !0.1
           theta=50

!          do j=0,ny
!             if (abs(j-ny/2) <= ny/20)  then
!               u(1,0,j)=rhoj
!               u(2,0,j)=rhoj*vj
!               if((j-ny/2).ge.0)then
!                 u(3,0,j)=0.5*rhoj*vj
!               else
!                 u(3,0,j)=-0.5*rhoj*vj
!               end if
!               u(4,0,j)=0.5*rhoj*vj**2 + Pj/(gamma-1.)
!             end if
!           end do

          do j=0,nx
             if (abs(j-nx/2) <= nx/20)  then
               u(1,j,0)=rhoj
               if((j-nx/2).ge.0)then
                 u(2,j,0)= 0.5*rhoj*vj     !*COS(45*3.14*time/180)
               else
                 u(2,j,0)=-0.5*rhoj*vj   !*COS(45*3.14*time/180)
               end if
               u(3,j,0)=1.0*rhoj*vj   !*sin(45*3.14*time/180)
               u(4,j,0)=0.5*rhoj*vj**2 + Pj/(gamma-1.)
             end if
           end do


         end if
!---------------------------------------
!  End jet injection
!---------------------------------------

  return
end subroutine boundaries
!=======================================================================




!subrutina de desacoplamiento
subroutine uprim(q,p)

 ! use parameters, only : neq, gamma
  use globals
  implicit none
  real :: q(neq), p(neq)
  real :: w, m2, lor
  real :: u2, alpha, chi
  integer :: k
  m2 = sum(q(2:3)**2)               ! v^2

  call newrap(q, w, m2)

  alpha = m2 / w**2   ! alpha < 1 !
  u2  = alpha/(1d0-alpha)

  lor = sqrt(1d0 + u2)

  !  velocities
  p(2:3) = q(2:3) / w

   !determination of the mass density
  p(1) = q(1)/lor

   !thermal pressure
  chi = (w - q(1)*(1d0+u2/(lor+1d0)))/(1d0+u2)

  p(4)   = (gamma - 1d0)/gamma * chi

  p(4) = max(p(4),1d-10*p(4))

  if(p(1).lt.0d0) then
    print*,'negative density in uprim'
    stop
  end if

end subroutine uprim
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
subroutine newrap(q, w, m2)

  !use parameters, only : neq, gamma
  use globals
  implicit none

  real, intent(in) :: q(neq)
  real, intent(in)  :: m2
  real, intent(out) :: w
  real, parameter :: eps = 1d-10
  real :: a, b, c, mu, alpha, u2, lor, chi, dpdchi, dpdrho
  real :: w0, dpdw, f, dfdw, pg, dv2dw, dchidw, drhodw !q(neq)
  integer :: i,j,k

  if(q(4)**2.lt.m2+q(1)**2.or.q(4).le.0d0) then
    print*,'error in newrap'
    stop
  end if

  a = 3d0;   b = 2d0 * (-q(4));   c = m2
  if(b**2-a*c.lt.0d0) then
     print*,'b**2-a*c<0'
     stop
  end if
  w = ( - b + sqrt(b**2-a*c) ) / a   ! initial guess for w = rho * h * lor**2

  w0 = w
  mu = 1d0
100 continue
  do k = 1, 40

    alpha = m2 / w**2   ! alpha < 1 !

    u2  = alpha/(1d0-alpha)

    if(u2.lt.0d0) then
      print*,'u2<0cc'
      print*,u
      stop
    end if

    lor = sqrt(1d0 + u2)

    chi = (w - q(1)*(1d0+u2/(lor+1d0)))/(1d0+u2)

    ! ideal gas case
    pg     = (gamma - 1d0)/gamma * chi
    dpdchi = (gamma - 1d0)/gamma
    dpdrho = 0d0

    f = w - pg - q(4)  ! f(w) = 0

    if(abs(f).lt.eps) return

    dv2dw  = lor/w**3*m2
    dchidw = 1d0/lor**2 + dv2dw*(q(1)+2d0*lor*chi)

    drhodw = dv2dw*q(1)

    dpdw = dpdchi*dchidw + dpdrho*drhodw

    dfdw   = 1d0 - dpdw    ! df/dw

    w  = w0 - mu * f / dfdw                       ! Newton-raphson iteration

    if(abs(w-w0).lt.eps) return

    w0 = w

  end do

  if(mu.gt.0.1) then
    mu = mu/2d0
    goto 100
  end if

end subroutine newrap
!------------------------------------------------------------------------------!
