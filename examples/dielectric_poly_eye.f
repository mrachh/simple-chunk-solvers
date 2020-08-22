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
      real *8 parsall(4,1000)
      real *8, allocatable :: chunks(:,:,:),ders(:,:,:),ders2(:,:,:)
      real *8, allocatable :: hs(:)
      real *8, allocatable :: whts(:,:)
      real *8, allocatable :: xtarg(:)
      real *8, allocatable :: ytarg(:)
      real *8  verts(2,1000)
      integer iverts(2,100)
      integer irbdry(10,4)
      integer nsegs(10)
      integer irside(100)
      integer ilside(100)
      integer, allocatable :: idt(:)
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
      zks(3) = 1.7d0 + 0.0d0*eye
      zks(4) = 1.6d0 + 0.0d0*eye

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
      chsmall = 1.0d-2

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
c     read in vertices from file
c     line 1   nverts
c     line 2   x,y  of first vertex
c     ...
c     line nverts+1  x,y  of last vertex
c     line nverts+2  nreg
c     line nverts+3  nseg
c     line nverts+4  i1  i2       irl irr (verts for first seg, left/right region)
c     line nverts+5  i1  i2       irl irr (verts for second seg, left/right region)
c     ...
c     line nverts+3+nseg  i1  i2  irl irr (vertices for last seg, left/right region)
c     line nverts+3+nseg+1    nsegs
c     
c
c     chunk up each side of domain 
c     routine fcurve1 defines a side based on 
c     x in [pars(1),pars(2), y in [opars(3)mparrs(4)] 
c  
      open (unit=19,file='input_eye.dat',status='unknown')
c
      read(19,*) nverts
      do i = 1,nverts
         read(19,*) verts(1,i), verts(2,i)
      enddo
      read(19,*) nreg
      read(19,*) nseg
      do ii = 1,nseg
         read(19,*) iverts(1,ii),iverts(2,ii),ilside(ii),irside(ii)
      enddo
      do ireg = 1,nreg
         read(19,*) nsegs(ireg)
         do ii = 1,nsegs(ireg)
            read(19,*) irbdry(ii,ireg)
         enddo
      enddo
c
      write(6,*) 'nverts = ',nverts
      write(6,*) 'xverts = ',(verts(1,i),i=1,nverts)
      write(6,*) 'yverts = ',(verts(2,i),i=1,nverts)
      write(6,*) 'nreg = ',nreg
      write(6,*) 'nseg = ',nseg
      write(6,*) 'iverts = ',(iverts(1,i),i=1,nseg)
      write(6,*) 'iverts = ',(iverts(2,i),i=1,nseg)
      write(6,*) 'ilside = ',(ilside(i),i=1,nseg)
      write(6,*) 'irside = ',(irside(i),i=1,nseg)
      write(6,*) 'nsegs = ',(nsegs(i),i=1,nreg)
      write(6,*) 'irbdry 1',(irbdry(ii,1),ii=1,nsegs(1))
      write(6,*) 'irbdry 2',(irbdry(ii,2),ii=1,nsegs(2))
      write(6,*) 'irbdry 3',(irbdry(ii,3),ii=1,nsegs(3))
      write(6,*) 'irbdry 4',(irbdry(ii,4),ii=1,nsegs(4))
c
c
      do i = 1,nverts
         write(6,*) verts(1,i), verts(2,i)
      enddo
c
c
      nch = 0
      do iii = 1,nseg
         parsall(1,iii) = verts(1,iverts(1,iii))
         parsall(2,iii) = verts(1,iverts(2,iii))
         parsall(3,iii) = verts(2,iverts(1,iii))
         parsall(4,iii) = verts(2,iverts(2,iii))
         call chunkfunczk(zks(0),npw,eps,ifclosed,irefl,irefr,chsmall,
     1       ta,tb,fcurve1,parsall(1,iii),nover,k,
     1       nch1,chunks(1,1,nch+1),adjs(1,nch+1),
     2       ders(1,1,nch+1),ders2(1,1,nch+1),hs(nch+1))
         do i = 1,nch1
            iregionl(nch+i) = ilside(iii)
            iregionr(nch+i) = irside(iii)
         enddo
         nch = nch+nch1
          call prinf('nch1=*',nch1,1)
c
      enddo
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
      a3 = dcmplx(1.6d0,0.06d0)
      b3 = dcmplx(1.08d0,0.02d0)
ccc      a2 = a1
ccc      b2 = b1
      alphas(0) = a0
      alphas(1) = a1
      alphas(2) = a2
      alphas(3) = a3
      alphas(4) = a3
      betas(0) = b0
      betas(1) = b1
      betas(2) = b2
      betas(3) = b3
      betas(4) = b3
c
c     testing just getdomains
c
ccc      if (2.ne.3) goto 222

c
c  get the matrix
c

      call cpu_time(t1)
C$       t1 = omp_get_wtime()

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
        call prinf('ich is *',ich,1)
        call prin2('zkl is *',zkl,2)
        call prin2('zkr is *',zkr,2)
        call prin2('al is *',al,2)
        call prin2('ar is *',ar,2)
        call prin2('bl is *',bl,2)
        call prin2('br is *',br,2)
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
ccc      stop
c
      call cpu_time(t1)
C$       t1 = omp_get_wtime()

      call zgmres_solver(2*n,xmat,zid,zrhs,numit,eps,niter,errs,rres,
     1   zsoln)

      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      
      call prin2('matrix solve time=*',t2-t1,1)

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
c
c
c
222   continue
c
c
      hh = 0.1d0
      next = 0
      nnx = 200
      nny = 200
c
      allocate(xtarg(nnx*nny))
      allocate(ytarg(nnx*nny))
      allocate(idt(nnx*nny))
c
      do j = 1,nny
      do i = 1,nnx
         next = next+1
         xtarg(next) = -10.0000000001d0 + i*hh
         ytarg(next) = -15.0000000001d0 + j*hh
      enddo
      enddo

      nsmax = 10
      write(6,*) 'in main nsmax ',nsmax
      write(6,*) 'in main nreg',nreg
      write(6,*) 'in main nseg',nseg
      write(6,*) 'in main nsegs',(nsegs(ii),ii=1,nreg)
      write(6,*) 'nnx ',nnx
      write(6,*) 'nny ',nny
      do iseg = 1,nseg
         write(6,*) 'in main parsall',(parsall(i,iseg),i=1,4)
         write(13,*) 'in main parsall',(parsall(i,iseg),i=1,4)
      enddo
c
      call domainflag(nsmax,nreg,nseg,nsegs,parsall,irbdry,xtarg,ytarg,
     1     nnx,nny,idt)
 
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
