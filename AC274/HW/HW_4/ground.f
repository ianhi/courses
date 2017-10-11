c   Solve for the ground state of the 1d harmonic oscillator with Schroedinger in imaginary-time
c   centered differences in space 
c   periodic boundary conditions (thru buffers at j=0 and j=NX+1)
c /////////////////////////////////////////////////////////////////////
       parameter (NX=101,NM=4)
       dimension x(NX),f(0:NX+1),fnew(0:NX+1)
       dimension pot(0:NX+1)
       dimension fmom(0:NM),hmom(0:NM)
c --------------------------------------
       open(31,file='fmom.out')
       open(32,file='hmom.out')
       open(33,file='probe.out')     ! wath the local deacy to ground state
       open(36,file='movie.out')
c baseline parameters and input data
       Nsteps = 10000
       nout   = 200

       sizex  = 10.0
       dx     = sizex/float(nx-1)
       dif    = 1.0
       vel    = 0.0
       vamp   = 1.0   ! omega square
c -----------------------------------------------------
c CFL numbers fix the timestep
       cfl   = 0.10
       dtd   = 1.0
       if(dif.ne.0.0) dtd = dx*dx/dif
       dta   = 1.0
       if(vel.ne.0.0) dta = dx/vel
       dtc   = 1.0
       if(vamp.ne.0.0) dtc = 1.0/(vamp*sizex*sizex)
c pick up the safest dt           (play with them)
       dt    = cfl*min(min(dta,dtd),dtc)
       write(6,*) dx,dtd,dta,dtc,dt

       alfa   = vel*dt/dx
       delta  = dif*dt/(dx*dx)
c initial conditions: unnormalized  (play wity them)
       wid   = 1.0
       do j=1,NX+1
        x(j)   = dx*float(j-1-nx/2)
        xj     = x(j)/wid
c mixture of ground and excited state
        f(j)   = (1.0+0.5*xj)*exp(-0.5*xj*xj)
        pot(j) = 0.5*vamp*x(j)*x(j)
        potn = pot(j)/(0.25*vamp*sizex*sizex)
        write(36,*) x(j),f(j),potn
       end do
       write(36,'(bn)')
       write(36,'(bn)')
c time-stepper ---------------------------------
       do it=1,Nsteps
        time = (it-1)*dt
        f(0)    = f(nx)
        f(nx+1) = f(1)
        do j=1,NX
c build the tridiagonal coefficients
         a = delta+0.5*alfa
         b = delta-0.5*alfa
         c = 1.0-a-b-pot(j)*dt
         fnew(j) = a*f(j-1)+c*f(j)+b*f(j+1)
        end do
c prepare for new timestep
        do j=1,NX
         f(j) = fnew(j)
        end do
c -----------------------------------------------------------------------
c diagnostics: standard moments, tail moment and entropies
        if(mod(it,nout).eq.1) then
         do m=0,NM
          fmom(m) = 0.0
          hmom(m) = 0.0
         end do
         do j=1,NX
          fj = f(j)
          fmom(0) = fmom(0) + fj*dx
          fmom(1) = fmom(1) + fj*x(j)*dx
          fmom(2) = fmom(2) + fj*x(j)*x(j)*dx
          hmom(2) = hmom(2) + fj*fj*dx
c output for gnuplot movies
          write(36,*) x(j),f(j),pot(j)/pot(nx)
         end do
c blank lines for gnuplot movies
         write(36,'(bn)')
         write(36,'(bn)')
c monitor the moments in time and local probes for expaonential decay
         fmean = fmom(1)/fmom(0)
         fvar  = fmom(2)/fmom(0)-fmean*fmean
         write(31,*), it,fmom(0),fmean,fvar,xmom/fmom(0)
         write(32,*), it,hmom(2),hmom(3),hmom(4)
c teh probes: in the long-time they should decay like exp(-omega_0*t)
         write(33,*), time,f(nx/4),f(nx/2),f(3*nx/4)
        endif
        write(6,*) 'total mass and entropy',it,fmom(0),hmom(2)
c end evolution -----------------------------------------------------
       end do

       stop
       end
