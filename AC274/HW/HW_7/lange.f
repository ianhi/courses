       program lange
c      solve the langevin equation 
c      dv/dt=-gamma(v)*v + a(x) + noise 
c      AC274 Fall 2017
c =============================================================================
       implicit double precision (a-h,o-z)
       parameter(nsteps=50000)
       dimension vv(nsteps),xx(nsteps)
       dimension vvcf(0:nsteps/2),xxcf(0:nsteps/2),xvcf(0:nsteps/2)
c ---------------------------------------
        open(7,file  ='lange.inp')
        open(9,file  ='traj.out')
        open(11,file ='cfun.out')
        open(21,file ='sfun.out')
        open(51,file ='bin.out')
        open(52,file ='bbin.out')
c ----------------------------------------
c number of realizations, initial conditions, random seed
        read(7,*) nreal
        
        a0     = 30.0         ! external field  
        wk0    = 0.5          ! spatial periodicity
        b      = 0.3          ! exponential scale
        gamma  = 1.0        
        vth    = 0.1
        rnl    = 0.0          ! non linear damping rate

        do ir=1,nreal
        read(7,*) x0,iseed
        dt    = 0.01
        nout  = 50
c initial conditions
        x = x0
        v = 0.0
c time evolution
        vrave=0.
        vrvar=0.
        do it=1,nsteps
         rv = (2.0*ranpang(iseed)-1.0)
c single well, linear
         if(rnl.eq.0.0) then
          z  = gamma*dt
          gd = exp(-z)
          gc = (1.0-exp(-z))/z
          gr = sqrt(0.5*(1.0-exp(-2.*z))/z)
c gamma-integrator
          xk = 6.28*wk0*x
          a  = -a0*sin(xk)*cos(xk)*exp(b*abs(x))
          vc = a*dt
          vr = rv*vth*sqrt(3.0)

          vrave=vrave+vr
          vrvar=vrvar+vr*vr
cd          write(6,*) 'gd,gc,gr',gd,gc,gr
cd          write(6,*) 'v,vc,vr,x',v,vc,vr,x
          vnew = gd*v + gc*vc + gr*vr
         else
c double well: v=0 unstable, v=\pm 1 stable
c plain euler because we don't have the analytical solution of the purely damped problem 
          z = gamma*dt
c note the sign change because v=0 must be unstable
          xk = 6.28*k0*x
          a  = -a0*sin(xk)*cos(xk)*exp(b*x)
          vnew = v*(1.0+z*(1.0-rnl*v*v)) + a*dt + vth*rv*sqrt(dt)
         endif
         xnew = x + 0.5*(v+vnew)*dt
         if(mod(it,nout).eq.1) then
          write(9,*) it,xnew,vnew
         endif
         vv(it) = vnew
         xx(it) = xnew
c new step
         v = vnew
         x = xnew
        end do       ! end time evolution
        write(6,*) 'random v:average and variance'
     $             ,vrave/nsteps, (vrave/nsteps)**2-vrvar/nsteps
        end do       ! end realizations
c ---------------
c data analysis
c ---------------
        call bin(vv,nsteps)     ! simple binning
        write(6,*) 'out of binning'
        call bbin(vv,nsteps)    ! piecewise linear binning
        write(6,*) 'out of b-binning'
c        call cfun(xx,vv,nsteps) ! correlation functions
        write(6,*) 'out of correlation function'
c        call sfun(xx,vv,nsteps) ! structure functions
        write(6,*) 'out of structure function'

        stop 
        end
c ====================================
        subroutine cfun(xx,vv,N)
c ====================================
        implicit double precision (a-h,o-z)
        dimension xx(N),vv(N)
        dimension xxcf(0:N/4),vvcf(0:N/4),xvcf(0:N/4)
c  ---------------------------------------------
        do j=0,N/4-1
         vvcf(j)=0.0
         xxcf(j)=0.0
         xvcf(j)=0.0
         do i=1,3*N/4
          vvcf(j)=vvcf(j)+vv(i)*vv(i+j)
          xxcf(j)=xxcf(j)+xx(i)*xx(i+j)
          xvcf(j)=xvcf(j)+0.5*(xx(i)*vv(i+j)+vv(i)*xx(i+j))
         end do 
         cfvv = vvcf(j)/vvcf(0)
         cfxx = xxcf(j)/xxcf(0)
         cfxv = xvcf(j)/xvcf(0)
         write(11,*) j,cfvv,cfxx,cfxv
        end do 

        return
        end
c ====================================
        subroutine sfun(xx,vv,N)
c ====================================
        implicit double precision (a-h,o-z)
        parameter (np=6)
        dimension xx(N),vv(N)
        dimension sf(N/4,np)
c  ---------------------------------------------
        do j=1,N/4
         do i=1,3*N/4
          do ip=1,np
           sf(j,ip) = 0.0
          end do
          dx = abs(vv(i+j)-vv(i))
          do ip=1,np
           sf(j,ip) = sf(j,ip)+ dx**ip
          end do
         end do
         write(21,*) j,(sf(j,ip),ip=1,np)
        end do
        
        stop
        end
c ====================================
        subroutine bin(f,N)
c ====================================
        implicit double precision (a-h,o-z)
        parameter (Nbin=101)
        dimension f(N),b(Nbin),icount(Nbin)
c  ---------------------------------------------
        fmax = -1.e14
        fmin = - fmax
        do j=1,N
         if(f(j).gt.fmax) fmax=f(j)
         if(f(j).lt.fmin) fmin=f(j)
        end do
        write(6,*) 'min-max', fmin,fmax

        do k=1,Nbin
         b(k) = fmin+(k-1)*(fmax-fmin)/float(Nbin-1)
         icount(k) = 0
        end do

        do j=1,N
         do k=1,Nbin-1
          if (f(j).gt.b(k).and.f(j).lt.b(k+1)) then
           icount(k) = icount(k)+1
          endif
         end do
        end do

        do k=1,Nbin
         write(51,*) b(k),icount(k)
        end do
        
        return
        end
c ====================================
          subroutine bbin(f,N)
c better bin with pwl charge assignement
c ====================================
          implicit double precision (a-h,o-z)
          parameter (Nbin=101)
          dimension f(N),b(Nbin),fcount(Nbin)
c -------------------------------------------
          fmax = -1.e14
          fmin = - fmax
          do j=1,N
            if(f(j).gt.fmax) fmax=f(j)
            if(f(j).lt.fmin) fmin=f(j)
          end do

          h =(fmax-fmin)/float(Nbin-1)
          do i=1,Nbin
           b(i) = fmin+(i-1)*h
           fcount(i) = 0.0
          end do
c inner bins
          do j=1,N
           do i=2,Nbin-2
            if (f(j).gt.b(i).and.f(j).lt.b(i+1)) then
             x  = (f(j)-b(i))/h
             wl = 0.5*(1-x)*(1-x)       
             wr = x*x
             wc = 1.0 - wl - wr
             fcount(i-1) = fcount(i-1)+wl
             fcount(i)   = fcount(i)+wc
             fcount(i+1) = fcount(i+1)+wr
            endif
           end do
          end do
c boundary bins
          if(f(j).gt.b(1).and.f(j).lt.b(2)) then
            fcount(1) = fcount(1)+0.5 !to be improved
            fcount(2) = fcount(2)+0.5
           endif 
           if(f(j).gt.b(Nbin-1).and.f(j).lt.b(Nbin)) then
            fcount(Nbin-1) = fcount(Nbin-1)+0.5    ! to be improved
            fcount(Nbin)   = fcount(Nbin)+0.5
           endif 
           do i=1,Nbin
            write(52,*) b(i), fcount(i)
           end do
          
           return
           end
c ============================          
         function ranpang(iseed)
c Pang, p.47         
        implicit double precision (a-h,o-z)
c ============================          
         i2e30 = 2**30
         ia=16807
         ic=i2e30-1    ! ic=2**31-1, but 2**31 is a overflow
         ic=ic+i2e30
          
         iq=ic/ia
         ir=mod(ic,ia)

         ih=iseed/iq
         il=mod(iseed,iq)

         it=ia*il-ir*ih
         if(it.gt.0) then
           iseed=it
         else
           iseed=it+ic
         endif

         ranpang=iseed/float(ic)

         return
         end
