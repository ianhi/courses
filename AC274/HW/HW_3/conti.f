c /////////////////////////////////////////////////////////////////////
c Solve 1d Continuity Equation with "low-order' FD schemes:
c 1. centered FD :  unstable
c 2. Upwind:     : No wiggles, but overdiffusive
c 3. Lax-Wendroff: Gibbs wiggles on the trail side
c S. Succi, AC274, Fall semester 2017
c ----------------------------------------------------------
       parameter (NX=100,NM=4)
       dimension x(NX),u(0:NX+1),f(0:NX+1),fnew(0:NX+1)
       dimension fmom(0:NM),hmom(0:NM)
c --------------------------------------
       open(31,file='fmom.out')
       open(32,file='hmom.out')
       open(36,file='movie.out')
c --------------------------------------
c baseline parameters and input data
c   ifd =1,2,3: centered,upwind,lax-wendroff 
c   nuni=0,1 uniform/parabolic velocity profile

       write(6,*) 'enter ifd and nuni'
       read(5,*) ifd,nuni
       write(6,*) ifd,nuni

       Nsteps = 500
       nout   = 25    ! output every nout steps
       sizex  = 1.0
       dif    = 0.0
       vel    = 0.1
       dx     = sizex/float(nx-1)
c CFL numbers fix the timestep
       cfl   = 0.1
       dtd   = 1.0
       if(dif.ne.0.0) dtd = 0.5*dx*dx/dif
       if(vel.ne.0.0) dta = dx/vel
       dt = cfl*min(dta,dtd)
       write(6,*) 'dx and dt',dx,dt
       dnum   = dx*dx/dt
       delta  = dif/dnum
c initial conditions: unnormalized  
       sigma = 8.0*dx            
       xm = 0.5*(sizex-sigma)
       xp = 0.5*(sizex+sigma)
       do j=1,NX
         x(j) = dx*float(j-1)
        f(j) = 0.0
        if(x(j).gt.xm.and.x(j).lt.xp) f(j) = 1.0   !box
        xj   = (x(j)-sizex/2)/(sizex/2)
c imposed parabolic flow-field 
        u(j) = vel               ! uniform
        if(nuni.gt.0) u(j) = vel*(1-xj*xj)
        write(36,*) x(j),f(j),u(j)
       end do
       write(36,'(bn)')
       u(0)    = u(nx)
       u(nx+1) = u(1)
c time-stepper ---------------------------------
       do it=1,Nsteps
        time = dt*it
        f(0)    = f(nx)
        f(nx+1) = f(1)
        do j=1,NX
         alfal  = u(j-1)*dt/dx
         alfac  = u(j  )*dt/dx
         alfar  = u(j+1)*dt/dx
c build the tridiagonal coefficients 
         if(ifd.eq.1) then             ! Centered
          a = delta+0.5*alfal
          b = delta-0.5*alfar
          c = 1.0-2*delta
         endif
c  -----------------------------
         if(ifd.eq.2) then             ! Upwind
          if(u(j).gt.0) then
           a = alfal+delta      ! f(j)-f(j-1)
           b = delta
           c = 1.-alfac-2.0*delta
          else                  ! f(j+1)-f(j)
           a = delta
           b = delta-alfar
           c = 1.0+alfac-2.0*delta
          endif
         endif
c ---------------------------------
         if(ifd.eq.3) then             ! Lax-Wendroff
c centered conservative (unstable)
          a =  0.5*alfal+0.5*alfac*alfac + delta
          b = -0.5*alfar+0.5*alfac*alfac + delta
          c = 1.0-alfac*alfac-2.0*delta 
         endif
c tridiagonal update
         fnew(j) = a*f(j-1)+c*f(j)+b*f(j+1)
        end do
c ---------------------------------------------
c prepare for new timestep
        do j=1,NX
         f(j) = fnew(j)
        end do
c diagnostics: standard moments, tail moment and entropies
        if(mod(it,nout).eq.0) then
         do m=0,NM
          fmom(m) = 0.0
          hmom(m) = 0.0
         end do
         xmom = 0.0
         do j=1,NX
          fj      = f(j)
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
         write(31,*), time,fmom(0),fmean,fvar,xmom/fmom(0)
         write(32,*), time,hmom(2),hmom(3),hmom(4)
        endif
        write(6,*) 'completed time step',it,fmom(0)
c end evolution -----------------------------------------------------
       end do

       stop
       end
