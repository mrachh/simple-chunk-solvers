c
c     driver solves dielectric interface problem
c     on multicomponent domain specified by input file 
c     in standard format (see below).
c
c     A simple example is:
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
c     The input file states: 
c     1) number of regions
c     2) number of segments (each defined by a seq. of points)
c     3) number of points on each segment and x,y coordinates of those points.
c     4) integer parameters determining whether to refine at "left" (initial) 
c        endpoint or "right" (final) endpoint
c     5) integer parameters stating subdomain/subregion identitiy on left/right
c        of each segment
c     6) number of segments defining boundary of each subregion
c     7) segment numbers defining boundary of each subregion along
c        positively oriented path
c     8) zk0 -> Helmholtz parameter for exterior (vacuum) region
c     9) indices of refraction for each subregion
c    10) angle of orientation of incoming plane wave
c    11) polarization (ipol = 0 -> TE, ipol=1 -> TM)
c
      implicit real *8 (a-h,o-z)

      parameter (maxc = 10000)
 
      complex *16 zks(0:10),alphas(0:10),betas(0:10),eye
      complex *16 rns(0:10)
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
ccc      zks(0) = 1.0d1
ccc      zks(1) = 1.4d1 + 0.0d0*eye
cccccc      zks(1) = 1.0d0
ccc      zks(2) = 1.3d1 + 0.0d0*eye
ccc      zks(3) = 1.4d1 + 0.0d0*eye
ccc      zks(4) = 1.4d1 + 0.1d0*eye
c
      nover = 1
      npw = 20
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
ccc         call prin2(' kurv is *',kurv,npts2)
c
         read(5,*) ilside(ii),irside(ii)
         read(5,*) irefinel(ii),irefiner(ii)
         write(6,*) 'ii is ',ii
         write(6,*) 'ilside, irside is ',ilside(ii), irside(ii)
         write(6,*) 'irefl, irefr is ',irefinel(ii), irefiner(ii)
      enddo
      do ireg = 1,nreg
         read(5,*) nsegs(ireg)
         do ii = 1,nsegs(ireg)
            read(5,*) irbdry(ii,ireg)
         enddo
         write(6,*) 'irbdry is ',(irbdry(ii,ireg),ii=1,nsegs(ireg))
      enddo
      read(5,*) xxx,yyy
      write(6,*) ' zk0 is ',zks(0)
      zks(0) = xxx+ eye*yyy
      do ireg = 1,nreg
         read(5,*) rn1,rn2
         rns(ireg) = rn1+ eye*rn2
         zks(ireg) = rns(ireg)*zks(0)
      enddo
      read(5,*) thetain
      read(5,*) ipol
      write(6,*) ' zk1 is ',zks(1)
      write(6,*) ' zk2 is ',zks(2)
      write(6,*) ' thetain is ',thetain
      write(6,*) ' ipol is ',ipol
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
c
      xyin1(1) = 0.3d0
      xyin1(2) = 0.2d0
      xyin2(1) = 0.2d0
      xyin2(2) = -0.3d0
      xyout(1) = 6.3d0
      xyout(2) = 6.2d0
c
c     true scattering problem: ifin = 1, iscat = 1
c
      ifin = 0
      iscat = 0


      ntot = k*nch
      allocate(xmat(2*ntot,2*ntot))
      
      t1 = second()


      allocate(srcvals(8,k,nch),srccoefs(6,k,nch),whts(k,nch))

      call chunks_to_srcinfo(k,nch,chunks,ders,ders2,hs,srcvals,
     1  srccoefs,whts)

      n = k*nch
c
ccc      a0 = dcmplx(1.0d0,0.1d0)
ccc      b0 = dcmplx(1.2d0,0.04d0)
ccc      a1 = dcmplx(1.3d0,0.07d0)
ccc      b1 = dcmplx(1.4d0,0.03d0)
ccc      a1 = a0
ccc      a2 = dcmplx(1.1d0,0.06d0)
ccc      b2 = dcmplx(1.8d0,0.02d0)
ccc      a2 = a0
ccc      b2 = b0
ccc      a3 = dcmplx(1.6d0,0.06d0)
ccc      b3 = dcmplx(1.08d0,0.02d0)
ccc      a2 = a1
ccc      b2 = b1
      alphas(0) = 1.0d0
      betas(0) = 1.0d0
      do ireg = 1,nreg
         alphas(ireg) = 1.0d0
         betas(ireg) = 1.0d0
         if (ipol.eq.1) betas(ireg) = 1.0d0/(rns(ireg)**2)
      enddo
c
c     if treu scattering problem, angle of incidence of plane wave
c     theta = 0  -> normal,   0 < theta < pi/2 is rightgoing
c                         -pi/2 < theta < 0 is leftgoing
c
ccc      thetain = 0.0d0
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
ccc        write(6,*) 'ich,zkl,zkr',ich,zkl,zkr
ccc        write(13,*) 'ich,zkl,zkr',ich,zkl,zkr
ccc        write(6,*) 'irr,irl',iregionr(ich),iregionl(ich)
ccc        write(13,*) 'irr,irl',iregionr(ich),iregionl(ich)
ccc        call prinf('ich is *',ich,1)
ccc        call prin2('zkl is *',zkl,2)
ccc        call prin2('zkr is *',zkr,2)
ccc        call prin2('al is *',al,2)
ccc        call prin2('ar is *',ar,2)
ccc        call prin2('bl is *',bl,2)
ccc        call prin2('br is *',br,2)
ccc        call prin2('zq is *',zq,2)
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
ccc           call prinf('iscat is *',iscat,1)
           if (iscat.eq.0) then 
              call h2d_slp(xsr,8,srcvals(1,j,ich),0,dpars,1,zkr,0,
     1          ipars,zz)
              zrhs(j,ich,1) = zz*ar/zq
ccc           if (ich.eq.1) write(47,*) 'regionr is ',iregionr(ich)
ccc           if (ich.eq.1) write(47,*) 'regionl is ',iregionl(ich)
ccc           if (ich.eq.1) write(47,*) 'xsr is ',xsr
ccc           if (ich.eq.1) write(47,*) 'xsl is ',xsl
ccc           if (ich.eq.1) write(47,*) 'j,zrhs 1st ',j,zrhs(j,ich,1)
          
              call h2d_slp(xsl,8,srcvals(1,j,ich),0,dpars,1,zkl,0,
     1          ipars,zz)
             zrhs(j,ich,1) = zrhs(j,ich,1) - zz*al/zq
ccc           if (ich.eq.1) write(47,*) 'j,zrhs 2nd ',j,zrhs(j,ich,1)
           else
              zrhs(j,ich,1) = 0.0d0
              if (iregionr(ich).eq.0) then
                 zrhs(j,ich,1) = zrhs(j,ich,1) 
     1                   -ar*uinfun(srcvals(1,j,ich),zkr,thetain)/zq
ccc              call prin2('rhs =*',zrhs(j,ich,1),2)
              endif
              if (iregionl(ich).eq.0) then
                 zrhs(j,ich,1) = zrhs(j,ich,1) 
     1                   +al*uinfun(srcvals(1,j,ich),zkl,thetain)/zq
              endif
           endif
        enddo
      enddo
c
      do ich=1,nch
ccc        call prinf(' iregionr = *',iregionr(ich),1)
ccc        call prinf(' iregionl = *',iregionl(ich),1)
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
           if (iscat.eq.0) then 
              call h2d_sprime(xsr,8,srcvals(1,j,ich),0,dpars,1,zkr,0,
     1          ipars,zz)
              zrhs(j,ich,2) = zz*br
              call h2d_sprime(xsl,8,srcvals(1,j,ich),0,dpars,1,zkl,0,
     1          ipars,zz)
              zrhs(j,ich,2) = zrhs(j,ich,2) - zz*bl
           else
              zrhs(j,ich,2) = 0.0d0
              if (iregionr(ich).eq.0) then
                 zrhs(j,ich,2) = zrhs(j,ich,2) 
     1                   -br*uinfun_p(srcvals(1,j,ich),zkr,thetain)
              endif
              if (iregionl(ich).eq.0) then
                 zrhs(j,ich,2) = zrhs(j,ich,2) 
     1                   +bl*uinfun_p(srcvals(1,j,ich),zkl,thetain)
              endif
ccc              zrhs(j,ich,2) = 0.0d0
           endif
        enddo
      enddo

ccc      zid = b2/2
      zid = 1.0d0
      numit = 200
      niter = 0

      eps = 1.0d-9
c
ccc      call prin2(' zrhs is *',zrhs,2*n*2)
ccc      stop
c
      call cpu_time(t1)
C$       t1 = omp_get_wtime()

      call zgmres_solver(2*n,xmat,zid,zrhs,numit,eps,niter,errs,rres,
     1   zsoln)

      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      
      call prin2('matrix solve time=*',t2-t1,1)

ccc      call prin2(' sol is *',zsoln,2*n*2)
      do ich = 1,nch
         do j = 1,k
            write(33,*) abs(zsoln(j,ich,1) - zrhs(j,ich,1))
            write(33,*) abs(zsoln(j,ich,2) - zrhs(j,ich,2))
         enddo
      enddo
c
c
c    test soln at an exterior point
c
      targ(1) = 5.31d0
      targ(2) = -4.01d0
c
ccc      write(6,*) 'targ,source,zk ',targ,xyin,zks(0)

      call  prin2(' xyin1 *',xyin1,2)
      call  prin2(' targ *',targ,2)
      call h2d_slp(xyin1,2,targ,0,dpars,1,zks(0),0,ipars,potex)
ccc      potex = uinfun(targ,zks(0),theta)
      call prin2('potex=*',potex,2)

      
      pot = 0
      pots = 0
      potd = 0
      do ich=1,nch
        do j=1,k
          zz = 0
          call h2d_dlp(srcvals(1,j,ich),2,targ,0,dpars,1,zks(0),0,
     1      ipars,zz)
          pot = pot + zsoln(j,ich,1)*zz*whts(j,ich)/betas(0)
          potd = potd + zsoln(j,ich,1)*zz*whts(j,ich)/betas(0)
ccc          call prin2('zsoln=*',zsoln(j,ich,1),2)
ccc          call prin2('potd=*',potd,2)
          call h2d_slp(srcvals(1,j,ich),2,targ,0,dpars,1,zks(0),0,
     1      ipars,zz)
          pot = pot - zsoln(j,ich,2)*zz*whts(j,ich)/betas(0)
          pots = pots - zsoln(j,ich,2)*zz*whts(j,ich)/betas(0)
cc          call prin2('zz=*',zz,2)
        enddo
      enddo

      call prin2('pots=*',pots,2)
      call prin2('potd=*',potd,2)
      call prin2('pot=*',pot,2)
      call prin2('abs error=*',abs(potex-pot),1)
      call prin2('rel error=*',abs(potex-pot)/abs(potex),1)
c
c
c    test soln at an interior point
c
      targ(1) = 0.0d0
      targ(2) = 0.2d0
c
      call h2d_slp(xyout,2,targ,0,dpars,1,zks(1),0,ipars,potex)

      
      pot = 0
      do ich=1,nch
        do j=1,k
          zz = 0
          call h2d_dlp(srcvals(1,j,ich),2,targ,0,dpars,1,zks(1),0,
     1      ipars,zz)
          pot = pot + zsoln(j,ich,1)*zz*whts(j,ich)/betas(1)
          call h2d_slp(srcvals(1,j,ich),2,targ,0,dpars,1,zks(1),0,
     1      ipars,zz)
          pot = pot - zsoln(j,ich,2)*zz*whts(j,ich)/betas(1)
cc          call prin2('zz=*',zz,2)
        enddo
      enddo

      call prin2(' region 1 targ pot=*',pot,2)
      call prin2(' region 1 targ potex=*',potex,2)
      call prin2('abs error=*',abs(potex-pot),1)
      call prin2('rel error=*',abs(potex-pot)/abs(potex),1)
c
c
c
c
222   continue
c
c
      hhx = 0.1d0
      hhy = 0.15d0
      next = 0
      nnx = 200
      nny = 200
ccc      nnx = 2
ccc      nny = 2
c
      allocate(xtarg(nnx*nny))
      allocate(ytarg(nnx*nny))
      allocate(idt(nnx*nny))
      allocate(upotall(nnx*nny))
      allocate(potexall(nnx*nny))
c
      do j = 1,nny
      do i = 1,nnx
         next = next+1
         xtarg(next) = -10.02d0 + i*hhx
         ytarg(next) = -15.02d0 + j*hhy
      enddo
      enddo
ccc      xtarg(1) = 0.0d0
ccc      ytarg(1) = 0.2d0
ccc      xtarg(2) = 0.0d0
ccc      ytarg(2) = 0.2d0
ccc      xtarg(3) = 0.0d0
ccc      ytarg(3) = 0.2d0
ccc      xtarg(4) = 0.0d0
ccc      ytarg(4) = 0.2d0
ccc      xtarg(1) = 0.31d0
ccc      ytarg(1) = -0.00001d0

      nsmax = 10
ccc      write(6,*) 'in main nsmax ',nsmax
ccc      write(6,*) 'in main nreg',nreg
ccc      write(6,*) 'in main nseg',nseg
ccc      write(6,*) 'in main nsegs',(nsegs(ii),ii=1,nreg)
ccc      write(6,*) 'nnx ',nnx
ccc      write(6,*) 'nny ',nny
ccc      do iseg = 1,nseg
ccc         write(6,*) 'in main parsall',(parsall(i,iseg),i=1,4)
ccc         write(13,*) 'in main parsall',(parsall(i,iseg),i=1,4)
ccc      enddo
c
ccc      call domainflag(nsmax,nreg,nseg,nsegs,parsall,irbdry,xtarg,ytarg,
ccc     1     nnx,nny,idt)
      nnn = nch*k
      ntarget = nnx*nny
      call prinf(' calling domainflagsm nch is *',nch,1)
      t0 = second()
      call domainflagsm(k,nch,nnn,srccoefs,srcvals,
     1      iregionl,iregionr,nreg,ntarget,xtarg,ytarg,idt) 
      t1 = second()
      call prin2(' time for domainflag is *',t1-t0,1)
ccc      do i = 1,ntarget
ccc         idt(i) = 0
ccc         rr = dsqrt(xtarg(i)**2 + ytarg(i)**2)
ccc         if (rr .lt.2) idt(i) = 1
ccc         if ((rr .lt.2) .and. (ytarg(i).lt.0.0d0)) idt(i) = 2
ccc      enddo
      open(unit = 20,file = 'energy.m',status='unknown');
      write(20,*) 'edens = ['
      do i = 1,ntarget
         write(20,*) xtarg(i), ytarg(i), idt(i)
      enddo
      write(20,*) '];'
      write(20,*) 'X = reshape(edens(:,1),',nnx,',',nny,')'
      write(20,*) 'Y = reshape(edens(:,2),',nnx,',',nny,')'
      write(20,*) 'Z = reshape(edens(:,3),',nnx,',',nny,')'
      write(20,*) 'mesh(X,Y,Z)'
ccc      call prinf( 'idt is *',idt,4)
c
ccc      call prinf(' ntarget = *',ntarget,1)
ccc      call prinf(' nreg = *',nreg,1)
      call prinf(' calling evalpot = *',n,1)
      t1 = second()
      call evalpotvol(k,nch,n,srccoefs,srcvals,zks,alphas,betas,
     1  zsoln(1,1,2),zsoln(1,1,1),nreg,idt,ntarget,xtarg,ytarg,upotall,
     2  ifin,uinfun,thetain)
ccc      call prinf(' after eval vol = *',ntarget,1)
ccc      call prin2( 'upotall is *',upotall,2*4)
      call prinf( 'idt(1) is *',idt(1),1)
      call prin2( 'upotall(1) is *',upotall(1),2)
c
      do i=1,ntarget
         ireg = idt(i)
         zkl = zks(ireg)
         if (ireg .eq. 0) then
            xsl(1) = xyin1(1)
            xsl(2) = xyin1(2)
         else
            xsl(1) = xyout(1)
            xsl(2) = xyout(2)
         endif
         targ(1) = xtarg(i)
         targ(2) = ytarg(i)
         if (iscat.eq.0) then
            call h2d_slp(xsl,2,targ,0,dpars,1,zkl,0,ipars,potexall(i))
         else
            potexall(i) = uinfun(targ,zkl,theta)
         endif
      enddo
ccc      call prin2( 'potexall is *',potexall,2*4)
c
c 
      open(unit = 19,file = 'pot.m',status='unknown');
      open(unit = 27,file = 'potex.m',status='unknown');
ccc      open(unit = 20,file = 'energy.m',status='unknown');
      write(19,*) 'rpot = ['
ccc      write(20,*) 'edens = ['
      write(27,*) 'epot = ['
      do i = 1,ntarget
         write(19,*) xtarg(i), ytarg(i), dreal(upotall(i))
         write(27,*) xtarg(i), ytarg(i), dreal(potexall(i))
ccc         write(20,*) xtarg(i), ytarg(i), cdabs(upotall(i))**2
ccc         write(20,*) xtarg(i), ytarg(i), idt(i)
      enddo
      t2 = second()
      call prin2(' time to eval at targs *',t2-t1,1)
      write(19,*) '];'
      write(19,*) 'X = reshape(rpot(:,1),',nnx,',',nny,')'
      write(19,*) 'Y = reshape(rpot(:,2),',nnx,',',nny,')'
      write(19,*) 'Z = reshape(rpot(:,3),',nnx,',',nny,')'
      write(19,*) 'mesh(X,Y,Z)'
ccc      write(20,*) '];'
ccc      write(20,*) 'X = reshape(edens(:,1),',nnx,',',nny,')'
ccc      write(20,*) 'Y = reshape(edens(:,2),',nnx,',',nny,')'
ccc      write(20,*) 'Z = reshape(edens(:,3),',nnx,',',nny,')'
ccc      write(20,*) 'mesh(X,Y,Z)'
      write(27,*) '];'
      write(27,*) 'X = reshape(epot(:,1),',nnx,',',nny,')'
      write(27,*) 'Y = reshape(epot(:,2),',nnx,',',nny,')'
      write(27,*) 'Z = reshape(epot(:,3),',nnx,',',nny,')'
      write(27,*) 'mesh(X,Y,Z)'
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
c
        function uinfun(targ,zk,theta)
        implicit none
        real *8 targ(2),theta,dot
        complex *16 eye,zk, uinfun
        data eye/(0.0d0,1.0d0)/
c
        dot = dsin(theta)*targ(1) - dcos(theta)*targ(2)
        uinfun = cdexp(eye*zk*dot)
        return
        end
c
        function uinfun_p(srcval,zk,theta)
        implicit none
        real *8 srcval(8),theta,dot
        complex *16 eye,zk,uinfun_p,gx,gy
        data eye/(0.0d0,1.0d0)/
c
        dot = dsin(theta)*srcval(1) - dcos(theta)*srcval(2)
        gx = eye*zk*dsin(theta)*cdexp(eye*zk*dot)
        gy = -eye*zk*dcos(theta)*cdexp(eye*zk*dot)
        uinfun_p = gx*srcval(7) + gy*srcval(8)
        return
        end


