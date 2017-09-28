c /////////////////////////////////////////////////////////////////////
c   Solve the 1d continuity equation using Flux-Corrected Transport SHASTA 
c   periodic boundary conditions (thru buffers at j=0 and j=NX+1)
c   S. Succi, AC274, Fall semester 2017
c ///////////////////////////////////////////////////////////////////
       parameter (NX=101,NM=4)
       dimension x(NX),u(-1:NX+2),f(-1:NX+2),fnew(-1:NX+2)
       dimension fmom(0:NM),hmom(0:NM)
c --------------------------------------
       open(31,file='fmom.out')
       open(32,file='hmom.out')
       open(36,file='movie.out')
c ---------------------------------------------
c baseline parameters and input data
c  0: no antidif, no limiter, 1: antidif only, 2: antidif+limiter    
       write(6,*) 'enter shasta and uni/nuni flow option'
       read(5,*)  lop,nuni
       write(6,*) lop,nuni

       Nsteps = 1001
       nout   = 100    ! output every nout steps
       sizex  = 1.0
       dif    = 0.0
       adif   =-0.125
       vel    = 0.5
       dx     = sizex/float(nx-1)
c CFL numbers fix the timestep
       cfl   = 0.3
       dtd   = 1.0
       if(dif.ne.0.0) dtd = 0.5*dx*dx/dif
       if(vel.ne.0.0) dta = dx/vel

       dt = cfl*min(dta,dtd)
       write(6,*) 'dx,dta,dtd,dt',dx,dta,dtd,dt
       dnum   = dx*dx/dt
       delta  = dif/dnum

c initial conditions: unnormalized  (play wity them)
       sigma = 5.0*dx
       xm = 0.5*sizex-sigma
       xp = 0.5*sizex+sigma
       do j=1,NX
        x(j) = dx*float(j-1)
        f(j) = 0.0
        if(x(j).gt.xm.and.x(j).lt.xp) f(j) = 1.0

        xj   = (x(j)-sizex/2)/(sizex/2)
c gaussian
c        r2 =0.5* ((x(j)-0.5*sizex)/sigma)**2 
c        f(j) = exp(-r2)

        u(j) = vel
        if(nuni.ne.0) u(j) = vel*(1-xj*xj)
       end do
c periodicity, note extra buffers because shasta requires double-neighbor info
       u(-1)   = u(nx-1)
       u(0)    = u(nx)
       u(nx+1) = u(1)
       u(nx+2) = u(2)
c time-stepper ---------------------------------
       do it=1,Nsteps
        time = dt*it
        f(-1)   = f(nx-1)
        f(0)    = f(nx)
        f(nx+1) = f(1)
        f(nx+2) = f(2)
        do j=1,NX
c lagrangian move plus interpolation back to the grid
         alfal  = u(j-1)*dt/dx
         alfac  = u(j)*dt/dx
         alfar  = u(j+1)*dt/dx
         ql = (0.5+alfac)/(1.0-(alfal-alfac))
         qr = (0.5-alfac)/(1.0+(alfar-alfac))
         a = 0.5*ql*ql
         b = 0.5*qr*qr
         c = ql+qr-a-b
c pure transport, no flux correction, on a uniform grid this should be Lax-Wendroff plus 1/8 diffusion
         fnew(j) = a*f(j-1)+c*f(j)+b*f(j+1)

c anti-diffusion options
         if (lop.eq.1) then
          fluxr = adif*(f(j+1)-f(j))
          fluxl = adif*(f(j-1)-f(j))
          fnew(j) = fnew(j) + fluxr + fluxl
          write(6,*) 'lop',lop, fluxl,fluxr
         endif
c --------------------------------------------------------
         if (lop.eq.2) then
c ***************************** WARNING !!!!
c  This option is based on eq. 23 of Boris-Book paper
c  I did not re-derive myself from scratch. Boris-Book only gives the 
c rightward flux (fluxr) expression, so there might be an 
c ambiguity/error in the way I constructed the left flux (fluxl)
c Their notation leaves several ambiguities on the signs of the fluxes
c I have no guarantee that the implementation is correct
c In fact it does not seem to be because i see no or less
c wiggles but also some mass non-conservation.
c ****************************************
c flux limiter: compact form
          dl1 = f(j)-f(j-1)
          dl2 = f(j-1)-f(j-2)
          dr1 = f(j+1)-f(j)
          dr2 = f(j+2)-f(j+1)
c sign switches, danger on the tails
          sr1 = 0.0
          if(abs(dr1).gt.0.0) sr1=dr1/abs(dr1)
          sl1 = 0.0
          if(abs(dl1).gt.0.0) sl1=dl1/abs(dl1)
c min(x,y,z) = min(x,min(y,z))
          difa  = abs(adif)
          fminr = min( min(dl1*sl1,difa*abs(dr1)),dr2*sr1)
          fluxr = sr1*max(0.0,fminr)

          fminl = min( min(dl2*sl2,difa*abs(dl1)),dr1*sl1)
          fluxl = sl1*max(0.0,fminl)

          fnew(j) = fnew(j) + fluxl - fluxr

          write(6,*) 'j,dl1,dl2,dr1,dr2,sl1,sr1',j,dl1,dl2,dr1,dr2,sl1,sr1
          write(6,*) 'fnew(j),fluxl,fluxr',fnew(j),fluxl,fluxr

         endif
        
        end do             ! end loop in space
c ---------------------------------------------
c prepare for new timestep
        do j=1,NX
         f(j) = fnew(j)
        end do
c diagnostics: standard moments, tail moment and entropies
        if(mod(it,nout).eq.1) then
         do m=0,NM
          fmom(m)=0.0
          hmom(m)=0.0
         end do
         xmom = 0.0
         do j=1,NX
          fj = f(j)
          fmom(0) = fmom(0) + fj*dx
          fmom(1) = fmom(1) + fj*x(j)*dx
          fmom(2) = fmom(2) + fj*x(j)*x(j)*dx

          hmom(2) = hmom(2) - fj*fj*dx
          hmom(3) = hmom(3) - fj*fj*fj*dx
          hmom(4) = hmom(4) - fj*fj*fj*fj*dx
c probe the tail, should grow exponentially in time
          if(x(j).gt.0.99*sizex) xmom=xmom+fj*dx
c output for gnuplot movies
          write(36,*) x(j),f(j),u(j)
         end do
c blank lines for gnuplot movies
         write(36,'(bn)')
         write(36,'(bn)')
c monitor the moments in time
         fmean = fmom(1)/fmom(0)
         fvar  = fmom(2)/fmom(0)-fmean*fmean
         write(31,*), it,fmom(0),fmean,fvar,xmom/fmom(0)
         write(32,*), it,hmom(2),hmom(3),hmom(4)
        endif
        write(6,*) 'time step: mass and mean position',it,fmom(0),fmean
c end evolution -----------------------------------------------------
       end do

       stop
       end
