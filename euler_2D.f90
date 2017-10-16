!=======================================================================
!   This program solves the Euler equations (2D) with the Lax Method
!=======================================================================
!un cambio diferente en el programa
!otro cambio
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
      integer, parameter :: nx=100, ny=100, neq=4
      real, parameter :: xmax=1., dx=xmax/float(nx)
      real, parameter :: ymax=1., dy=ymax/float(ny)
      real, parameter :: gamma=1.6666666

      real, parameter :: tmax= 30.0        ! maximum integration time
      real, parameter :: dtprint=1.0       ! interval between outputs

      real, parameter :: rhoin = 10.0
      real, parameter :: rhoout = 0.001
      real, parameter :: pin = 100.0
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
      real, parameter :: bound=10.0

!------------------------------------------------------------------------------
!   This is a vector that contains u(x,y) (do not touch)
!------------------------------------------------------------------------------
      real :: u (neq,0:nx+1,0:ny+1)

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
          write (*,*) itprint+1,time,tmax,dt
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
      real :: x,y, rad

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

          if (rad < 0.0) then
            u(1,i,j)=rhoin
            u(2,i,j)=rhoin*vxin
            u(3,i,j)=rhoin*vyin
            u(4,i,j)=pin/(gamma-1.) + 0.5*u(2,i,j)*u(2,i,j)/u(1,i,j) + 0.5/u(1,i,j)*u(3,i,j)*u(3,i,j)
          else
            u(1,i,j)=rhoout
            u(2,i,j)=rhoout*vxout
            u(3,i,j)=rhoout*vyout
            u(4,i,j)=pout/(gamma-1.) + 0.5/u(1,i,j)*u(2,i,j)*u(2,i,j) + 0.5/u(1,i,j)*u(3,i,j)*u(3,i,j)

          end if

        end do
      end do

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
      real :: rho, vx, vy, P

      do i=0,nx+1
        do j=0,ny+1
          rho=u(1,i,j)
          vx=u(2,i,j)/rho
          vy=u(3,i,j)/rho
          P=(u(4,i,j)-0.5*rho*(vx**2+vy**2))*(gamma-1.)

          f(1,i,j)=rho*vx
          f(2,i,j)=rho*vx*vx+P
          f(3,i,j)=rho*vx*vy
          f(4,i,j)=vx*(u(4,i,j)+P)

          g(1,i,j)=rho*vy
          g(2,i,j)=rho*vx*vy
          g(3,i,j)=rho*vy*vy+P
          g(4,i,j)=vy*(u(4,i,j)+P)

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
      integer :: j , i
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
