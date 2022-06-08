c
c     test driver that solves dielectric interface problem
c     for point source in upper half space (zk0) 
c     and one in lower half space with dielectric (zk1) ->
c     artificial data 
c
c     us_0 = (-1/b0) S_k0[Gamma] sigma + (1/b0) D_k0[Gamma] mu 
c           
c     us_1 = (-1/b1) S_k1[Gamma] sigma + (1/b1) D_k1[Gamma] mu 
c
c     [au] = [a0 us_0 - a1 us_1] = [a0 uex0 - a1 uex1] 
c     [bu_n] = [b0 u_n_0 - b1 u_n_1] = [b0 uex0_n - b1 uex1_n] 
c
c     [  XMAT   ]  [ sig ]     [au] on Omega
c     [  XMAT   ]  [ mu  ]     [bu_n] on Omega 
c
c
c                        /\
c     |---------|-------/  \------|---------|---
c    -L-a      -a                 a         a+L
c
c     call real chunk code for [-a,a].
c     add complexified chunks to chunk list  after complex conversion of
c     chunks -> srcinfo
c
c
      implicit real *8 (a-h,o-z)

      parameter (maxc = 10000)
 
      complex *16 zks(0:10),alphas(0:10),betas(0:10),eye,zk0,zk1
      complex *16 zkl,zkr,al,ar,bl,br
      complex *16, allocatable :: xmat(:,:)
      real *8 pars(1000)
      real *8, allocatable :: chunks(:,:,:),ders(:,:,:),ders2(:,:,:)
      real *8, allocatable :: hs(:)
      real *8, allocatable :: whtsr(:,:)
      complex *16, allocatable :: whts(:,:)
      integer, allocatable :: adjs(:,:)
      integer, allocatable :: iregionl(:)
      integer, allocatable :: iregionr(:)
      real *8, allocatable :: srccoefs(:,:,:),srcvals(:,:,:)
      complex *16, allocatable :: zsrccoefs(:,:,:),zsrcvals(:,:,:)

      complex *16 pars1(10),pars2(10)

      complex *16 xyin(2),xyout(2),xsr(2),xsl(2)
ccc      real *8 srcinfoin(8)
ccc      real *8 srcinfoout(8)
      complex *16 grad(100)

      complex *16, allocatable :: zrhs(:,:,:)
      complex *16, allocatable :: zsoln(:,:,:)
      real *8 errs(1000),src_normal(2),targ_normal(2)
      complex *16 targ(2)
      complex *16 zz,pot,potex,zpars(6),zid
      complex *16 a0,a1,a2,b0,b1,b2
      complex *16 pots,potd,gx,gy
      complex *16 zpar3(3),zq

      external fcurve,fcurvel,fcurver,h2d_combc
      data eye/(0.0d0,1.0d0)/

      call prini(6,13)

      done = 1
      pi = 4*atan(done)

      zks(0) = 4.0d0 + 0.0d0*eye
      zks(1) = 5.0d0 + 0.0d0*eye
      zk0 = zks(0)
      zk1 = zks(1)

      eps = 1.0d-12

      k = 16

      nover = 2
      npw=10
      irefl = 0
      irefr = 0
      rl = 20.0d0
      pars(1) = 7.0d0
ccc      pars(1) = 0.0d0
      pars(2) = rl
      pars(3) = 3.0d0
      pars(4) = 0.0d0

      ier = 0
      ifclosed = 0
      chsmall = 1.0d0

      xyin(1) = 0.0d0
      xyin(2) = 1.2d0
      xyout(1) = 0.0d0
      xyout(2) = -1.2d0

      allocate(chunks(2,k,maxc),ders(2,k,maxc),ders2(2,k,maxc))
      allocate(adjs(2,maxc),hs(maxc))
      allocate(iregionl(maxc),iregionr(maxc))
c
c     routine fcurve1 defines smooth curve
c
      nch = 0
c
c     left end
c
      xa = -1.0d0
      ta = -rl
      tb = 0
      pars(4) = xa
      call chunkfunc(eps,ifclosed,chsmall,ta,tb,fcurvel,pars,nover,k,
     1       nch1,chunks,adjs,ders,
     2       ders2,hs)
c
c     right end
c
      xb = 1.0d0
      ta = 0
      tb = rl
      pars(4) = xb
      call chunkfunc(eps,ifclosed,chsmall,ta,tb,fcurver,pars,nover,k,
     1       nch2,chunks(1,1,nch1+1),adjs(1,nch1+1),ders(1,1,nch1+1),
     2       ders2(1,1,nch1+1),hs(nch1+1))
      nch = nch1+nch2
c
c     middle region
c
      ta = xa
      tb = xb
      call chunkfunc(eps,ifclosed,chsmall,ta,tb,fcurve,pars,nover,k,
     1       nch3,chunks(1,1,nch+1),adjs(1,nch+1),ders(1,1,nch+1),
     2       ders2(1,1,nch+1),hs(nch+1))
      nch = nch+nch3


       call prinf('nch1=*',nch1,1)
       call prinf('nch2=*',nch2,1)
       call prinf('nch3=*',nch3,1)
      open (unit=17,status='unknown',file='pts.m')
      write(17,*) ' xpts = ['
      do i = 1,nch
      do j = 1,k
         write(17,*) chunks(1,j,i), chunks(2,j,i)
      enddo
      enddo
      write(17,*) ' ];'
      write(17,*) ' plot(xpts(:,1),xpts(:,2),''.'')'
      ntot = k*nch
      nsys = ntot+nquad
      allocate(xmat(2*nsys,2*nsys))
      
      t1 = second()


      allocate(whtsr(k,nch))
      allocate(srcvals(8,k,nch),srccoefs(6,k,nch),whts(k,nch))
      allocate(zsrcvals(8,k,nch),zsrccoefs(6,k,nch))

      call chunks_to_srcinfoc(k,nch1+nch2,chunks,ders,ders2,hs,srcvals,
     1  srccoefs,whts)

      call chunks_to_srcinfo(k,nch3,chunks(1,1,nch1+nch2+1),
     1     ders(1,1,nch1+nch2+1),ders2(1,1,nch1+nch2+1),hs(nch1+nch2+1),
     2     srcvals(1,1,nch1+nch2+1),srccoefs(1,1,nch1+nch2+1),
     3     whtsr(1,nch1+nch2+1))

      n = k*nch
c
      a0 = 1.0d0
      b0 = 1.0d0
      a1 = 1.0d0
      b1 = 1.0d0
      zpars(1) = zks(0)
      zpars(2) = a0
      zpars(3) = b0
      zpars(4) = zks(1)
      zpars(5) = a1
      zpars(6) = b1
      zq = 0.5*(zpars(2)/zpars(3)+zpars(5)/zpars(6))
      call prin2(' zq is *',zq,2)
c
c  move ends to cmplx format
c
      do i=1,nch1+nch2
        do j=1,k
          zsrcvals(1,j,i) = srcvals(1,j,i) + eye*srcvals(2,j,i)
          zsrcvals(2,j,i) = 0.0d0
          zsrcvals(3,j,i) = srcvals(3,j,i) + eye*srcvals(4,j,i)
          zsrcvals(4,j,i) = 0.0d0
          zsrcvals(5,j,i) = srcvals(5,j,i) + eye*srcvals(6,j,i)
          zsrcvals(6,j,i) = 0.0d0
          srcvals(7,j,i) = 0.0d0
          srcvals(8,j,i) = -1.0d0
          zsrcvals(7,j,i) = srcvals(7,j,i)
          zsrcvals(8,j,i) = srcvals(8,j,i)
c      
          zsrccoefs(1,j,i) = srccoefs(1,j,i) + eye*srccoefs(2,j,i)
          zsrccoefs(2,j,i) = 0.0d0
          zsrccoefs(3,j,i) = srccoefs(3,j,i) + eye*srccoefs(4,j,i)
          zsrccoefs(4,j,i) = 0.0d0
          zsrccoefs(5,j,i) = srccoefs(5,j,i) + eye*srccoefs(6,j,i)
          zsrccoefs(6,j,i) = 0.0d0
        enddo
      enddo
c
      do i=nch1+nch2+1,nch1+nch2+nch3
        do j=1,k
          zsrcvals(1,j,i) = srcvals(1,j,i) 
          zsrcvals(2,j,i) = srcvals(2,j,i)
          zsrcvals(3,j,i) = srcvals(3,j,i) 
          zsrcvals(4,j,i) = srcvals(4,j,i)
          zsrcvals(5,j,i) = srcvals(5,j,i) 
          zsrcvals(6,j,i) = srcvals(6,j,i)
          zsrcvals(7,j,i) = srcvals(7,j,i)
          zsrcvals(8,j,i) = srcvals(8,j,i)
c      
          zsrccoefs(1,j,i) = srccoefs(1,j,i) 
          zsrccoefs(2,j,i) = srccoefs(2,j,i)
          zsrccoefs(3,j,i) = srccoefs(3,j,i) 
          zsrccoefs(4,j,i) = srccoefs(4,j,i)
          zsrccoefs(5,j,i) = srccoefs(5,j,i) 
          zsrccoefs(6,j,i) = srccoefs(6,j,i)
          whts(j,i) = whtsr(j,i)
        enddo
      enddo
c


      call cpu_time(t1)
C$       t1 = omp_get_wtime()

      call helm2d_dielec_matc(k,nch,n,zsrcvals,zsrccoefs,zpars,
     1     iregionl,iregionr,xmat)

      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      
      call prin2('matrix generation time=*',t2-t1,1)

      allocate(zrhs(k,nch,2),zsoln(k,nch,2))

      do ich=1,nch
        do j=1,k
           call h2d_slpc(xyin,8,zsrcvals(1,j,ich),0,dpars,1,zk0,0,
     1       ipars,zz)
           zrhs(j,ich,1) = zz*zpars(2)/zq
           call h2d_slpc(xyout,8,zsrcvals(1,j,ich),0,dpars,1,zk1,0,
     1       ipars,zz)
          zrhs(j,ich,1) = zrhs(j,ich,1) - zz*zpars(5)/zq
        enddo
      enddo

      do ich=1,nch
        do j=1,k
           call h2d_sprimec(xyin,8,zsrcvals(1,j,ich),0,dpars,1,zk0,0,
     1       ipars,zz)
           zrhs(j,ich,2) = zz*zpars(3)
           call h2d_sprimec(xyout,8,zsrcvals(1,j,ich),0,dpars,1,zk1,0,
     1       ipars,zz)
          zrhs(j,ich,2) = zrhs(j,ich,2) - zz*zpars(6)
        enddo
      enddo

ccc      zid = b2/2
      zid = 1.0d0
      numit = 200
      niter = 0

      eps = 1.0d-14

      call prin2(' zrhs1 is *',zrhs(1,1,1),k*nch*2)
      call prin2(' zrhs2 is *',zrhs(1,1,2),k*nch*2)
      call zgmres_solver(2*n,xmat,zid,zrhs,numit,eps,niter,errs,rres,
     1   zsoln)
      call prin2(' zsoln1 is *',zsoln(1,1,1),k*nch*2)
      call prin2(' zsoln2 is *',zsoln(1,1,2),k*nch*2)
      
      open (unit=18,status='unknown',file='densall.m')
      write(18,*) ' dens = ['
      do i = 1,nch
      do j = 1,k
         write(18,*) abs(zsoln(j,i,1)),abs(zsoln(j,i,2))
      enddo
      enddo
      write(18,*) ' ];'
      write(18,*) ' plot(xpts(:,1),dens(:,1),''.'')'
      write(18,*) ' hold on; plot(xpts(:,1),dens(:,2),''.'')'
c
c
c    test soln at an exterior point
c
      targ(1) = 0.31d0
      targ(2) = -3000.33d0

      call h2d_slpc(xyin,2,targ,0,dpars,1,zk0,0,ipars,potex)

      zpar3(1) = zk0
      zpar3(2) = -1.0d0/zpars(3)
      zpar3(3) = 1.0d0/zpars(3)
      
      pot = 0
      pots = 0
      potd = 0
      do ich=1,nch
        do j=1,k
          zz = 0
          call h2d_dlpc(zsrcvals(1,j,ich),2,targ,0,dpars,1,zk0,0,
     1      ipars,zz)
          pot = pot + zsoln(j,ich,1)*zz*whts(j,ich)/zpars(3)
ccc             pot = pot + zrhs(j,ich,1)*zz*whts(j,ich)/zpars(3)
             potd = potd + zsoln(j,ich,1)*zz*whts(j,ich)/zpars(3)
          call h2d_slpc(zsrcvals(1,j,ich),2,targ,0,dpars,1,zk0,0,
     1      ipars,zz)
          pot = pot - zsoln(j,ich,2)*zz*whts(j,ich)/zpars(3)
ccc          pot = pot - zrhs(j,ich,2)*zz*whts(j,ich)/zpars(3)
          pots = pots - zsoln(j,ich,2)*zz*whts(j,ich)/zpars(3)
cc          call prin2('zz=*',zz,2)
        enddo
      enddo

      call prin2('pots=*',pots,2)
      call prin2('potd=*',potd,2)
      call prin2('pot=*',pot,2)
      call prin2('potex=*',potex,2)
      call prin2('abs error=*',abs(potex-pot),1)
      call prin2('rel error=*',abs(potex-pot)/abs(potex),1)
 
      return
      end
c-----------------------------------------------------------
c
c
c
c
c
c

        subroutine fcurve(t,pars,x,y,dxdt,dydt,dxdt2,dydt2)
        implicit real *8 (a-h,o-z)
        real *8 pars(1)

        pi = 4.0d0*datan(1.0d0)
c
        cc = 1.0d0/pars(1)
        ff = 1.0d0
c
        x = t
        y = ff*exp(-x**2/cc**2)

        dxdt= 1.0d0
        dydt= -2*x*y/cc**2

        dxdt2=0.0d0
        dydt2= -2*y/cc**2 + 4*x*x*y/cc**4
c
        return
        end
c
c
c
c
c

        subroutine fcurveall(t,pars,x,y,dxdt,dydt,dxdt2,dydt2)
        implicit real *8 (a-h,o-z)
        real *8 pars(1)

        pi = 4.0d0*datan(1.0d0)
c
        cc = pars(1)
        rl = pars(2)
        a = pars(3)
c
        x = t
        y = -cc*(erfc((t+rl)/a)-erfc((rl-t)/a))

        dxdt= 1.0d0
        dydt=cc*(exp(-(t+rl)**2/a**2)+exp(-(rl-t)**2/a**2))*2/sqrt(pi)/a

        dxdt2=0.0d0
        dydt2= -(cc*4/(sqrt(pi)*a**3))*
     1        ((t+rl)*exp(-(t+rl)**2/a**2)+(t-rl)*exp(-(rl-t)**2/a**2))
c
        return
        end
c
c
c

        subroutine fcurvel(t,pars,x,y,dxdt,dydt,dxdt2,dydt2)
        implicit real *8 (a-h,o-z)
        real *8 pars(1)

        pi = 4.0d0*datan(1.0d0)
c
        cc = pars(1)
        rl = pars(2)
        a = pars(3)
        aend = pars(4)
c
        x = t + aend
        y = -cc*erfc((t+rl)/a)
c
        dxdt= 1.0d0
        dydt=cc*(exp(-(t+rl)**2/a**2))*2/sqrt(pi)/a
c
        dxdt2=0.0d0
        dydt2= -(cc*4/(sqrt(pi)*a**3))*
     1        ((t+rl)*exp(-(t+rl)**2/a**2))
c
        return
        end
c
c
c
c
        subroutine fcurver(t,pars,x,y,dxdt,dydt,dxdt2,dydt2)
        implicit real *8 (a-h,o-z)
        real *8 pars(1)

        pi = 4.0d0*datan(1.0d0)
c
        cc = pars(1)
        rl = pars(2)
        a = pars(3)
        bend = pars(4)
c
        x = t + bend
        y = -cc*erfc((t+rl)/a)
        y = cc*erfc((rl-t)/a)

        dxdt= 1.0d0
        dydt=cc*exp(-(rl-t)**2/a**2)*2/sqrt(pi)/a

        dxdt2=0.0d0
        dydt2= -(cc*4/(sqrt(pi)*a**3))*
     1        ((t-rl)*exp(-(rl-t)**2/a**2))
c
c
        return
        end
c
c
c
c

c

c
        subroutine chunks_to_srcinfoc(k,nch,chunks,ders,ders2,hs,
     1     srcinfo,srccoefs,whts)
c
c  This subroutine returns the chunks info into the new format 
c  and stores it in the srcinfo, srccoefs arrays.
c
c
c
c      srcinfo(8,k,nch)
c        srcinfo(1:2,npts) x(t),y(t)
c        srcinfo(3:4,npts) dxdt(t),dydt(t)
c        srcinfo(5:6,npts) d2xdt2(t),d2ydt2(t)
c        srcinfo(7:8,npts) rnx(t),rny(t) = dydt/dsdt,-dxdt/dsdt
c
c      srccoefs(6,k,nch)
c        srcinfo(1:2,k,nch) xcoefs,ycoefs
c        srcinfo(3:4,k,nch) dxdtcoefs, dydtcoefs
c        srcinfo(5:6,k,nch) d2xdt2coefs,d2ydt2coefs
c
c  Note that in this representation, the chunks $\gamma_{\ell}$
c  are maps from [-1,1] as opposed to the original parameterization
c  
c  So the derivatives are scaled by hs(ich) and the second derivatives
c  are scaled by hs(ich)**2
c
        implicit real *8(a-h,o-z)
        integer k,nch
        real *8 chunks(2,k,nch),ders(2,k,nch),ders2(2,k,nch)
        real *8 srcinfo(8,k,nch),srccoefs(6,k,nch)
        real *8 hs(nch)
        complex *16 whts(k,nch),dsdt,eye
        real *8, allocatable :: ts(:),wts(:),umat(:,:),vmat(:,:)
        data eye/(0.0d0,1.0d0)/

        itype = 2
        allocate(ts(k),wts(k),umat(k,k),vmat(k,k))
        call legeexps(itype,k,ts,umat,vmat,wts)

        do ich=1,nch
          do j=1,k 
            srcinfo(1,j,ich) = chunks(1,j,ich)
            srcinfo(2,j,ich) = chunks(2,j,ich)
            srcinfo(3,j,ich) = ders(1,j,ich)*hs(ich)
            srcinfo(4,j,ich) = ders(2,j,ich)*hs(ich)
            srcinfo(5,j,ich) = ders2(1,j,ich)*hs(ich)**2
            srcinfo(6,j,ich) = ders2(2,j,ich)*hs(ich)**2
            dsdt = ders(1,j,ich) + eye*ders(2,j,ich)
            srcinfo(7,j,ich) = ders(2,j,ich)/dsdt
            srcinfo(8,j,ich) = -ders(1,j,ich)/dsdt
            whts(j,ich) = dsdt*hs(ich)*wts(j)
          enddo

          do j=1,k
            do l=1,6
              srccoefs(l,j,ich) = 0
            enddo
            do l=1,k
              do m=1,6
                srccoefs(m,j,ich) = srccoefs(m,j,ich) + 
     1              umat(j,l)*srcinfo(m,l,ich)
              enddo
            enddo
          enddo
        enddo
        
        return
        end



