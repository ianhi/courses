c units; space hlen=sqrt(D/2*omega), time 1/omega
c V(x) = 1/2 m*omega*omegax*x;  
c -L/2<x<L/2   L=sizex
c with dif=0 rho(x) is steady, A and B oscillate standing
c --------------------------------------------------------
       parameter (nx=1001)
       dimension Aold(nx),A(nx),Anew(nx)
       dimension Bold(nx),B(nx),Bnew(nx)
       dimension x(nx),rho(nx),pot(nx),qpot(nx)
       dimension hc(nx),hw(nx),he(nx)
       dimension jw(nx),je(nx)
c ------------------------------------
       open(9,file='diagno.out')
       open(21,file='rhoAB.out')
       open(31,file='movie.out')
       open(33,file='qmovie.out')
c ------------------ -------------------------
       pi     = 4.0*atan(1.0)
       nout   = 1000 
       write(6,*) '0,1,>1 = free particle,harmonic potential, square barrier'
       read(5,*) lop

       dif    = 1.0     ! hbar/m in natural units hbar=m=1
       omh    = 1.0     ! harmonic frequency
c the harmonic length hl = sqrt(dif/omh)
       hlen   = sqrt(0.5*dif/omh) 
       hvel   = sqrt(0.5*dif*omh)

       sigma  = hlen
       vel    = 1.0*hvel    ! initial speed of the wavepacket in units hlen/htim

       htim   = 2.0*pi/omh

       sizex  = 20.0*hlen   

       sizet  = 200*htim
       dx     = sizex/float(nx-1)
       dt     = 0.01/omh                    ! 1/100 oscillation time            
       nt     = sizet/dt

       if(lop.eq.0)  then
        vamp = 0.0    ! free particle
       else
        vamp = 1.0
       endif
         
       write(6,*) 'dx,dt prior to CFL',dx,dt

       dtd    = 1.0
       if(dif.ne.0)  dtd = 0.5*dx*dx/dif
       dta    = 1.0
       if(vel.ne.0)  dta = dx/abs(vel)

       cfl    = 0.5
       dt     = cfl*min(dta,dtd)

       delta  = dif*dt/(dx*dx)
       alfa   = vel*dt/dx
       omega  = vamp*dt

       if(alfa.gt.1.0)  stop 'dt too large, alfa'
       if(delta.gt.1.0) stop 'dt too large, delta'
       if(omega.gt.1.0) stop 'dt too large, omega'

       write(6,*) 'sizex,sigma',sizex,sigma
       write(6,*) 'dx,dt',dx,dt
       write(6,*) 'alfa,delta,omega',alfa,delta,omega
c init
       potmax = 0.5*omh*omh*(0.25*sizex*sizex)
       if(lop.gt.1) potmax=vamp

       j0 = (nx-1)/2
       do j=1,nx
        x(j) = -0.5*sizex+dx*(float(j-1))
        xj   = x(j)
        rho(j)  = exp(-0.25*xj*xj/sigma/sigma)
        Aold(j) = rho(j)*cos(vel*xj)
        Bold(j) = rho(j)*sin(vel*xj)
        rhoA    = Aold(j)*Aold(j)
        rhoB    = Bold(j)*Bold(j)

        if(lop.eq.0.or.lop.eq.1) then           ! free or harmonic
          pot(j) = vamp*0.5*omh*omh*xj*xj
        else
         pot(j)=0.0                           ! square barrier
         j1 = 3*nx/4
         j2 = j1+nx/10
         if(j.gt.j1.and.j.lt.j2) pot(j)=vamp
        endif
        write(21,*) xj,Aold(j),Bold(j),pot(j)
        write(31,*) xj,rhoA,rhoB,rho(j),pot(j)/potmax
       end do
       write(31,'(bn)')
       write(31,'(bn)')
c -------------------------------------
       write(6,*) 'time begins'
c Euler start-up
        do j=1,nx
c periodic bc
         je(j) = j+1
         if(je(j).gt.nx) je(j)=1
         jw(j) = j-1
         if(jw(j).lt.1)  jw(j)=nx
         omega = pot(j)*dt
c discretized Hamiltonian: once and for all for simplicity
         hc(j) =  delta+omega
         hw(j) = -delta/2.
         he(j) = -delta/2.
c Euler start-up
         hB    = hc(j)*Bold(j)+hw(j)*Bold(jw(j))+he(j)*Bold(je(j))
         A(j)  = Aold(j)+hB
         hA    = hc(j)*Aold(j)+hw(j)*Aold(jw(j))+he(j)*Aold(je(j))
         B(j)  = Bold(j)-hA
        end do
c now unleash the time-engine -------      
       do it=2,nt
        do j=1,nx
         hB = hc(j)*B(j)+hw(j)*B(jw(j))+he(j)*B(je(j))
         Anew(j) = Aold(j)+2.0*hB
        end do
        do j=1,nx
         hA = hc(j)*A(j)+hw(j)*A(jw(j))+he(j)*A(je(j))
         Bnew(j) = Bold(j)-2.0*hA

        end do
c prepare new timestep
        do j=1,nx
         Aold(j) = A(j) 
         Bold(j) = B(j) 
         A(j)    = Anew(j) 
         B(j)    = Bnew(j) 

         rho(j) = sqrt(A(j)*A(j)+B(j)*B(j))
        end do
c printout
        if(mod(it,nout).eq.1) then
         tmass = 0.0
         ave   = 0.0
         var   = 0.0
         epot  = 0.0
         ekin  = 0.0
         do j=1,nx
          rhoA = A(j)*A(j)
          rhoB = B(j)*B(j)
          if(rho(j).gt.1.e-3) then
           qpot(j)= sqrt(rho(je(j)))-2.0*sqrt(rho(j))+sqrt(rho(jw(j)))
           qpot(j)=-qpot(j)/(dx*dx*sqrt(rho(j)))
          else
           qpot(j)=0.0
          endif

          ue=atan(B(je(j))/A(je(j)))
          uw=atan(B(jw(j))/A(jw(j)))
          uvel=0.5*(ue-uw)/dx
          tmass  = tmass+ rho(j)*dx
          ave    = ave  + x(j)*rho(j)*dx
          var    = var  + x(j)*x(j)*rho(j)*dx

          dA   = 0.5*(A(je(j))-A(jw(j)))/dx
          dB   = 0.5*(B(je(j))-B(jw(j)))/dx
          ekin = ekin + 0.5*(dA*dA+dB*dB)*dx
          epot = epot + pot(j)*rho(j)*dx

          write(31,*) x(j),rhoA,rhoB,rho(j),pot(j)/potmax,uvel
          write(33,*) x(j),rho(j),qpot(j)
         end do
         ave = ave/tmass
         rms = sqrt(var/tmass -ave*ave)
         write(6,*) 'total mass',it,tmass,ave,rms,epot,ekin
         write(9,*) it,tmass,ave,rms,epot,ekin
         do iout=31,33
          write(iout, '(bn)') 
          write(iout, '(bn)') 
         end do

        endif
c ================================================
       end do
       stop
       end
