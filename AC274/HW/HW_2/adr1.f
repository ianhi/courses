c   Solve the 1d ADR equation using
c   centered differences in space and forward euler in time.
c   periodic boundary conditions (thru buffers at j=0 and j=NX+1)
c /////////////////////////////////////////////////////////////////////
       parameter (NX=100,NM=4)
       dimension x(NX),f(0:NX+1),fnew(0:NX+1)
       dimension fmom(0:NM),hent(0:NM)
c --------------------------------------
       open(31,file='moment.out')
       open(32,file='entropy.out')
       open(36,file='movie.out')
c baseline parameters and input data
       Nsteps = 10000
       nout   = 100

       sizex  = 1.0
       dx     = sizex/float(nx-1)

       dif    = 0.1
       vel    = 0.1
       cha    = 10.0
       chb    = 1.0
c chemical reaction is df/dt = cha*f-chb*f^2
c with cha and chb> 0 converges to cap=cha/chb
       Pe  = vel*sizex/dif          ! Peclet number
       Da  = cha*sizex*sizex/dif    ! Damkohler number (diffusive)
       write(6,*) 'Peclet and Damkohler',Pe,Da
c -----------------------------------------------------
c CFL numbers fix the timestep
       cfl   = 0.25
       dtd   = 1.0
       if(dif.ne.0.0) dtd = dx*dx/dif
       dta   = 1.0
       if(vel.ne.0.0) dta = dx/vel
       dtc   = 1.0
       if(cha.ne.0.0) dtc = 1.0/cha
c pick up the safest dt           (play with them)
       dt    = cfl*min(min(dta,dtd),dtc)
       write(6,*) dx,dtd,dta,dtc,dt

       alfa   = vel*dt/dx
       delta  = dif*dt/(dx*dx)
c initial conditions: unnormalized  (play wity them)
       wid   = 1
       sigma = wid*dx            !  play with width
       do j=1,NX
         x(j) = dx*float(j-1)
         xj   = (x(j)-sizex/2)/sigma
         f(j) = exp(-0.5*xj*xj)
       end do
c time-stepper ---------------------------------
       do it=1,Nsteps
        f(0)    = f(nx)
        f(nx+1) = f(1)
        do j=1,NX
c build the tridiagonal coefficients
         a = delta+0.5*alfa
         b = delta-0.5*alfa
c insert space dependent chemistry if you wish
         cshape = +1.0        
         chr = cshape*(cha-chb*f(j))*dt

         c = 1.0-a-b+chr
         fnew(j) = a*f(j-1)+c*f(j)+b*f(j+1)
        end do
c prepare for new timestep
        do j=1,NX
         f(j) = fnew(j)
        end do
c diagnostics: standard moments, tail moment and entropies
        if(mod(it,nout).eq.1) then
         do m=0,NM
          fmom(m) = 0.0
          hent(m) = 0.0
         end do
         xmom = 0.0
         do j=1,NX
          fj = f(j)*dx
          fmom(0) = fmom(0) + fj
          fmom(1) = fmom(1) + fj*x(j)
          fmom(2) = fmom(2) + fj*x(j)*x(j)

          hent(2) = hent(2) - fj*fj
          hent(3) = hent(3) - fj*fj*fj
          hent(4) = hent(4) - fj*fj*fj*fj
c probe the tail, should grow exponentially in time
          if(x(j).gt.0.99*sizex) xmom=xmom+fj
c output for gnuplot movies
          write(36,*) it,j,f(j)
         end do
c blank lines for gnuplot movies
         write(36,'(bn)')
         write(36,'(bn)')
c monitor the moments in time
         fmean = fmom(1)/fmom(0)
         fvar  = fmom(2)/fmom(0)-fmean*fmean
         write(31,*), it,fmom(0),fmean,fvar,xmom/fmom(0)
         write(32,*), it,hent(2),hent(3),hent(4)
        endif
        write(6,*) 'total mass at time step',it,fmom(0)
c end evolution -----------------------------------------------------
       end do

       stop
       end
