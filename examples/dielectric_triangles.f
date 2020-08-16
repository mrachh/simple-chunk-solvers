c
c     test driver that solves dielectric interface problem
c     on two triangles 
c
c              (0,1)
c               /\
c            2 /  \ 1
c             /    \    k_0
c            / k_1  \
c     (-1,0)/______5_\  (1,0)
c           \        /
c            \ k_2  /
c           3 \    / 4
c              \  /
c               \/
c              (0,-1)
c
      implicit real *8 (a-h,o-z)

      parameter (maxc = 10000)
 
      complex *16 zks(0:10),alphas(0:10),betas(0:10),eye
      complex *16 zkl,zkr,al,ar,bl,br
      complex *16, allocatable :: xmat(:,:)
      real *8 pars(1000)
      real *8, allocatable :: chunks(:,:,:),ders(:,:,:),ders2(:,:,:)
      real *8, allocatable :: hs(:)
      real *8, allocatable :: whts(:,:)
      integer, allocatable :: adjs(:,:)
      integer, allocatable :: iregionl(:)
      integer, allocatable :: iregionr(:)
      real *8, allocatable :: srccoefs(:,:,:),srcvals(:,:,:)

      complex *16 pars1(10),pars2(10)

      real *8 xyin1(2),xyin2(2),xyout(2),xsr(2),xsl(2)
      real *8 srcinfoin(8)
      real *8 srcinfoout(8)
      complex *16 grad(100)

      complex *16, allocatable :: zrhs(:,:,:)
      complex *16, allocatable :: zsoln(:,:,:)
      real *8 errs(1000),targ(2),src_normal(2),targ_normal(2)
      complex *16 zz,pot,potex,zpars(6),zid
      complex *16 a0,a1,a2,b0,b1,b2
      complex *16 zpar3(3),zq

      external fcurve1,h2d_comb
      data eye/(0.0d0,1.0d0)/

      call prini(6,13)

      done = 1
      pi = 4*atan(done)

      zks(0) = 1.2d0 + 0.0d0*eye
      zks(1) = 1.4d0 + 0.0d0*eye
      zks(2) = 1.6d0 + 0.0d0*eye

      eps = 1.0d-7
      ifclosed = 1
      ier = 0

      lused = 0
      k = 16

      nover = 1
      npw = 10
      irefl = 1
      irefr = 1
      ta = 0
      tb = 1.0d0
      ier = 0
      ifclosed = 0
      chsmall = 1.0d-4

      xyin1(1) = 0.3d0
      xyin1(2) = 0.2d0
      xyin2(1) = 0.2d0
      xyin2(2) = -0.3d0
      xyout(1) = 6.3d0
      xyout(2) = 6.2d0

      allocate(chunks(2,k,maxc),ders(2,k,maxc),ders2(2,k,maxc))
      allocate(adjs(2,maxc),hs(maxc))
      allocate(iregionl(maxc),iregionr(maxc))
c
c     chunk up each side of domain 
c     routine fcurve1 defines a side based on 
c     x in [pars(1),pars(2), y in [pars(3)mparrs(4)] 
c  
c     side 1
      nch = 0
      pars(1) = 1.0d0
      pars(2) = 0.0d0
      pars(3) = 0.0d0
      pars(4) = 1.0d0
      call chunkfunczk(zks(0),npw,eps,ifclosed,irefl,irefr,chsmall,
     1       ta,tb,fcurve1,pars,nover,k,nch1,chunks,adjs,ders,ders2,hs)
      do i = 1,nch1
         iregionl(nch+i) = 1
         iregionr(nch+i) = 0
      enddo
      nch = nch+nch1
       call prinf('nch1=*',nch1,1)
ccc      if (2.ne.3) goto 111
c
c     side 2
      pars(1) = 0.0d0
      pars(2) = -1.0d0
      pars(3) = 1.0d0
      pars(4) = 0.0d0
      call chunkfunczk(zks(0),npw,eps,ifclosed,irefl,irefr,chsmall,
     1       ta,tb,fcurve1,pars,nover,k,
     1       nch1,chunks(1,1,nch+1),adjs(1,nch+1),
     2       ders(1,1,nch+1),ders2(1,1,nch+1),hs(nch+1))
       call prinf('nch1=*',nch1,1)
      do i = 1,nch1
         iregionl(nch+i) = 1
         iregionr(nch+i) = 0
      enddo
      nch = nch+nch1
c
c     side 3
      pars(1) = -1.0d0
      pars(2) = 0.0d0
      pars(3) = 0.0d0
      pars(4) = -1.0d0
      call chunkfunczk(zks(0),npw,eps,ifclosed,irefl,irefr,chsmall,
     1       ta,tb,fcurve1,pars,nover,k,
     1       nch1,chunks(1,1,nch+1),adjs(1,nch+1),
     2       ders(1,1,nch+1),ders2(1,1,nch+1),hs(nch+1))
       call prinf('nch1=*',nch1,1)
      do i = 1,nch1
         iregionl(nch+i) = 2
         iregionr(nch+i) = 0
      enddo
      nch = nch+nch1
c
c     side 4
      pars(1) = 0.0d0
      pars(2) = 1.0d0
      pars(3) = -1.0d0
      pars(4) = 0.0d0
      call chunkfunczk(zks(0),npw,eps,ifclosed,irefl,irefr,chsmall,
     1       ta,tb,fcurve1,pars,nover,k,
     1       nch1,chunks(1,1,nch+1),adjs(1,nch+1),
     2       ders(1,1,nch+1),ders2(1,1,nch+1),hs(nch+1))
       call prinf('nch1=*',nch1,1)
      do i = 1,nch1
         iregionl(nch+i) = 2
         iregionr(nch+i) = 0
      enddo
      nch = nch+nch1
c
ccc      if (2.ne.3) goto 111
c     side 5
      pars(1) = -1.0d0
      pars(2) = 1.0d0
      pars(3) = 0.0d0
      pars(4) = 0.0d0
      call chunkfunczk(zks(0),npw,eps,ifclosed,irefl,irefr,chsmall,
     1       ta,tb,fcurve1,pars,nover,k,
     1       nch1,chunks(1,1,nch+1),adjs(1,nch+1),
     2       ders(1,1,nch+1),ders2(1,1,nch+1),hs(nch+1))
       call prinf('nch1=*',nch1,1)
      do i = 1,nch1
         iregionl(nch+i) = 1
         iregionr(nch+i) = 2
      enddo
      nch = nch+nch1
c
111   continue
       call prinf('nch=*',nch,1)
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
      allocate(xmat(2*ntot,2*ntot))
      
      t1 = second()


      allocate(srcvals(8,k,nch),srccoefs(6,k,nch),whts(k,nch))

      call chunks_to_srcinfo(k,nch,chunks,ders,ders2,hs,srcvals,
     1  srccoefs,whts)

      n = k*nch
c
      a0 = dcmplx(1.0d0,0.1d0)
      b0 = dcmplx(1.2d0,0.04d0)
      a1 = dcmplx(1.3d0,0.07d0)
      b1 = dcmplx(1.4d0,0.03d0)
      a2 = dcmplx(1.1d0,0.06d0)
      b2 = dcmplx(1.8d0,0.02d0)
ccc      a2 = a1
ccc      b2 = b1
      alphas(0) = a0
      alphas(1) = a1
      alphas(2) = a2
      betas(0) = b0
      betas(1) = b1
      betas(2) = b2

c
c  get the matrix
c


      call cpu_time(t1)
C$       t1 = omp_get_wtime()

      nreg = 2
      call helm2d_dielec_mat_nreg(k,nch,n,srcvals,srccoefs,zks,
     1     alphas,betas,nreg,iregionl,iregionr,xmat)

      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      
      call prin2('matrix generation time=*',t2-t1,1)

      allocate(zrhs(k,nch,2),zsoln(k,nch,2))

      do ich=1,nch
        zkl = zks(iregionl(ich))
        zkr = zks(iregionr(ich))
        al = alphas(iregionl(ich))
        ar = alphas(iregionr(ich))
        bl = betas(iregionl(ich))
        br = betas(iregionr(ich))
        zq = 0.5*(ar/br+al/bl)
        call prin2('zq is *',zq,2)
        do j=1,k
           if (iregionr(ich).eq.0) then
              xsr(1) = xyin1(1)
              xsr(2) = xyin1(2)
           else
              xsr(1) = xyout(1)
              xsr(2) = xyout(2)
           endif
           if (iregionl(ich).eq.0) then
              xsl(1) = xyin1(1)
              xsl(2) = xyin1(2)
           else
              xsl(1) = xyout(1)
              xsl(2) = xyout(2)
           endif
           call h2d_slp(xsr,8,srcvals(1,j,ich),0,dpars,1,zkr,0,
     1       ipars,zz)
           zrhs(j,ich,1) = zz*ar/zq
           if (ich.eq.1) write(47,*) 'regionr is ',iregionr(ich)
           if (ich.eq.1) write(47,*) 'regionl is ',iregionl(ich)
           if (ich.eq.1) write(47,*) 'xsr is ',xsr
           if (ich.eq.1) write(47,*) 'xsl is ',xsl
           if (ich.eq.1) write(47,*) 'j,zrhs 1st ',j,zrhs(j,ich,1)
          
           call h2d_slp(xsl,8,srcvals(1,j,ich),0,dpars,1,zkl,0,
     1       ipars,zz)
          zrhs(j,ich,1) = zrhs(j,ich,1) - zz*al/zq
           if (ich.eq.1) write(47,*) 'j,zrhs 2nd ',j,zrhs(j,ich,1)
        enddo
      enddo

      do ich=1,nch
        zkl = zks(iregionl(ich))
        zkr = zks(iregionr(ich))
        al = alphas(iregionl(ich))
        ar = alphas(iregionr(ich))
        bl = betas(iregionl(ich))
        br = betas(iregionr(ich))
        zq = 0.5*(ar/br+al/bl)
        do j=1,k
           if (iregionr(ich).eq.0) then
              xsr(1) = xyin1(1)
              xsr(2) = xyin1(2)
           else
              xsr(1) = xyout(1)
              xsr(2) = xyout(2)
           endif
           if (iregionl(ich).eq.0) then
              xsl(1) = xyin1(1)
              xsl(2) = xyin1(2)
           else
              xsl(1) = xyout(1)
              xsl(2) = xyout(2)
           endif
           call h2d_sprime(xsr,8,srcvals(1,j,ich),0,dpars,1,zkr,0,
     1       ipars,zz)
           zrhs(j,ich,2) = zz*br
           call h2d_sprime(xsl,8,srcvals(1,j,ich),0,dpars,1,zkl,0,
     1       ipars,zz)
          zrhs(j,ich,2) = zrhs(j,ich,2) - zz*bl
        enddo
      enddo

ccc      zid = b2/2
      zid = 1.0d0
      numit = 200
      niter = 0

      eps = 1.0d-15
c
      call prin2(' zrhs is *',zrhs,2*n*2)
      call zgmres_solver(2*n,xmat,zid,zrhs,numit,eps,niter,errs,rres,
     1   zsoln)
      call prin2(' sol is *',zsoln,2*n*2)
c
c
c    test soln at an exterior point
c
      targ(1) = 5.31d0
      targ(2) = 3.33d0
c
      write(6,*) 'targ,source,zk ',targ,xyin,zks(0)

      call h2d_slp(xyin1,2,targ,0,dpars,1,zks(0),0,ipars,potex)

      
      pot = 0
      do ich=1,nch
        do j=1,k
          zz = 0
          call h2d_dlp(srcvals(1,j,ich),2,targ,0,dpars,1,zks(0),0,
     1      ipars,zz)
          pot = pot + zsoln(j,ich,1)*zz*whts(j,ich)/betas(0)
          call h2d_slp(srcvals(1,j,ich),2,targ,0,dpars,1,zks(0),0,
     1      ipars,zz)
          pot = pot - zsoln(j,ich,2)*zz*whts(j,ich)/betas(0)
cc          call prin2('zz=*',zz,2)
        enddo
      enddo

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

        subroutine fcurve1(t,pars,x,y,dxdt,dydt,dxdt2,dydt2)
        implicit real *8 (a-h,o-z)
        real *8 pars(1)
c
c       polygonal segment defined by 
c       x in [pars(1),pars(2), y in [pars(3)mparrs(4)] 
c       assuming t in [0,1].
c
        x = pars(1) + t*(pars(2)-pars(1))
        y = pars(3) + t*(pars(4)-pars(3))

        dxdt=pars(2)-pars(1)
        dydt=pars(4)-pars(3)
c
        dxdt2=0.0d0
        dydt2=0.0d0
c
        return
        end
c
