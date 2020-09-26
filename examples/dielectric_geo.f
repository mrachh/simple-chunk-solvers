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
      real *8 parsall(100,1000)
      real *8, allocatable :: chunks(:,:,:),ders(:,:,:),ders2(:,:,:)
      real *8, allocatable :: hs(:)
      real *8, allocatable :: whts(:,:)
      real *8, allocatable :: xtarg(:)
      real *8, allocatable :: ytarg(:)
      real *8  verts(2,1000)
      real *8  xs(1000)
      real *8  ys(1000)
      real *8  kurv(1000)
      real *8  sval(1000)
      real *8  xcoeffs(1000)
      real *8  ycoeffs(1000)
      integer iverts(2,100)
      integer irbdry(10,4)
      integer nsegs(10)
      integer irside(100)
      integer ilside(100)
      integer irefinel(100)
      integer irefiner(100)
      integer, allocatable :: idt(:)
      integer, allocatable :: adjs(:,:)
      integer, allocatable :: iregionl(:)
      integer, allocatable :: iregionr(:)
      real *8, allocatable :: srccoefs(:,:,:),srcvals(:,:,:)
      real *8 w1(100), w2(100000)

      complex *16 pars1(10),pars2(10)

      real *8 xyin1(2),xyin2(2),xyout(2),xsr(2),xsl(2)
      real *8 srcinfoin(8)
      real *8 srcinfoout(8)
      complex *16 grad(100)

      complex *16, allocatable :: zrhs(:,:,:)
      complex *16, allocatable :: zsoln(:,:,:)
      complex *16, allocatable :: upotall(:)
      complex *16, allocatable :: potexall(:)
      real *8 errs(1000),targ(2),src_normal(2),targ_normal(2)
      complex *16 zz,pot,potex,zpars(6),zid
      complex *16 pots,potd
      complex *16 a0,a1,a2,b0,b1,b2
      complex *16 zpar3(3),zq
      complex *16 uinfun,uinfun_p
      external fcurve1,h2d_comb
      data eye/(0.0d0,1.0d0)/

      call prini(6,13)

      done = 1
      pi = 4*atan(done)

      eps = 1.0d-3
      ifclosed = 1

      lused = 0
      k = 16
      zks(0) = 1.0d0

      nover = 1
      npw = 10
      ifclosed = 0
      chsmall = 5.0d-3
c
      allocate(chunks(2,k,maxc),ders(2,k,maxc),ders2(2,k,maxc))
      allocate(adjs(2,maxc),hs(maxc))
      allocate(iregionl(maxc),iregionr(maxc))
c
c     read in vertices from file
c     line 1   nregions
c     line 2   nseg
c     then...
c     npts  (1st seg)
c       x,y  of vertices
c       iregionl, iregionr
c       irefinel, irefiner
c     npts  (2nd seg)
c       x,y  of vertices
c       iregionl, iregionr
c       irefinel, irefiner
c     ...
c     npts  (last seg)
c       x,y  of vertices
c       iregionl, iregionr
c       irefinel, irefiner
c     nsegs region 1
c      seg #s  (pos oriented)
c     nsegs region 2
c      seg #s  (pos oriented)
c     ...
c     nsegs region (last_ 
c      seg #s  (pos oriented)
c
ccc      open (unit=5,file='semicircs.dat',status='unknown')
c
      read(5,*) nreg
      read(5,*) nseg
      write(6,*) 'nreg is ',nreg
      write(6,*) 'nseg is ',nseg
      do ii = 1,nseg
         read(5,*) npts
         do j = 1,npts
            read(5,*) xs(j), ys(j)
         enddo
         call projectonpoly(xs,ys,npts,sval,xcoeffs,ycoeffs,w1,w2)
         parsall(1,ii) = npts
         parsall(2,ii) = sval(1)
         parsall(3,ii) = sval(npts)
c
         do j = 1,npts
            parsall(2*j+2,ii) = xcoeffs(j)
            parsall(2*j+3,ii) = ycoeffs(j)
         enddo
c
         npts2 = 2*npts
         hh = (sval(npts)-sval(1))/(npts2-1)
         do i = 1,npts2
            sv = hh*(i-1)
            call fcurve1(sv,parsall(1,ii),xx,yy,dxdt,dydt,
     1        dxdt2,dydt2)
            rn = (dxdt**2 + dydt**2)**(1.5d0)
            kurv(i) = (dxdt*dydt2 - dydt*dxdt2)/rn
         enddo
         call prin2(' kurv is *',kurv,npts2)
c
         read(5,*) ilside(ii),irside(ii)
         read(5,*) irefinel(ii),irefiner(ii)
         write(6,*) 'segment is ',ii
         write(6,*) 'ilside, irside is ',ilside(ii), irside(ii)
         write(6,*) 'irefl, irefr is ',irefinel(ii), irefiner(ii)
      enddo
      do ireg = 1,nreg
         read(5,*) nsegs(ireg)
         do ii = 1,nsegs(ireg)
            read(5,*) irbdry(ii,ireg)
         enddo
      enddo
c
c
      nch = 0
      do iii = 1,nseg
         call chunkfunczk(zks(0),npw,eps,ifclosed,
     1       irefinel(iii),irefiner(iii),chsmall,
     1       parsall(2,iii),parsall(3,iii),fcurve1,parsall(1,iii),
     1       nover,k,nch1,chunks(1,1,nch+1),adjs(1,nch+1),
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
      stop
c
      return
      end
c-----------------------------------------------------------
c
c
c
      subroutine fcurve1(t,pars,x,y,dxdt,dydt,dxdt2,dydt2)
      implicit real *8 (a-h,o-z)
      real *8 pars(*)
      real *8, allocatable :: xcoeffs(:)
      real *8, allocatable :: ycoeffs(:)
      real *8, allocatable :: xdiff(:)
      real *8, allocatable :: ydiff(:)
      real *8, allocatable :: xdiff2(:)
      real *8, allocatable :: ydiff2(:)
c
c       polygonal segment defined by 
c       x in [pars(1),pars(2), y in [pars(3)mparrs(4)] 
c       assuming t in [0,1].
c
      npts = nint(pars(1))
      a = pars(2)
      b = pars(3)
c
      allocate(xcoeffs(npts))
      allocate(ycoeffs(npts))
      allocate(xdiff(npts))
      allocate(ydiff(npts))
      allocate(xdiff2(npts))
      allocate(ydiff2(npts))
c
      do i = 1,npts
         xcoeffs(i) = pars(2*i+2)
         ycoeffs(i) = pars(2*i+3)
      enddo
c
ccc      call prin2(' xcoeffs *',xcoeffs,npts)
ccc      call prin2(' ycoeffs *',ycoeffs,npts)
      call chbdif(xcoeffs,xdiff,npts,a,b)
      call chbdif(ycoeffs,ydiff,npts,a,b)
      call chbdif(xdiff,xdiff2,npts,a,b)
      call chbdif(ydiff,ydiff2,npts,a,b)
ccc      call prin2(' xdiff *',xdiff,npts)
ccc      call prin2(' xdiff2 *',xdiff2,npts)
ccc      call prin2(' ydiff *',ydiff,npts)
ccc      call prin2(' ydiff2 *',ydiff2,npts)
c
      call cheval(t,x,xcoeffs,npts,a,b)
      call cheval(t,y,ycoeffs,npts,a,b)
c
      call cheval(t,dxdt,xdiff,npts,a,b)
      call cheval(t,dydt,ydiff,npts,a,b)
c
      call cheval(t,dxdt2,xdiff2,npts,a,b)
      call cheval(t,dydt2,ydiff2,npts,a,b)
c
      return
      end
c

