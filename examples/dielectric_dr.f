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
      parameter (maxseg = 10)
      parameter (maxreg = 4)
      parameter (maxp = 10)
 
      complex *16 zks(0:10),alphas(0:10),betas(0:10),eye
      complex *16 rns(0:10)
      complex *16 zkl,zkr,al,ar,bl,br
      real *8 parsall(2*maxp+3,maxseg)
      real *8, allocatable :: chunks(:,:,:),ders(:,:,:),ders2(:,:,:)
      real *8, allocatable :: hs(:)
      real *8, allocatable :: whts(:,:)
      real *8, allocatable :: xtarg(:)
      real *8, allocatable :: ytarg(:)
      real *8  xsall(maxp,maxseg)
      real *8  ysall(maxp,maxseg)
      integer irbdry(maxseg,4)
      integer nsegs(maxseg)
      integer irside(maxseg)
      integer ilside(maxseg)
      integer irefinel(maxseg)
      integer irefiner(maxseg)
      integer, allocatable :: idt(:)
      integer, allocatable :: adjs(:,:)
      integer, allocatable :: iregionl(:)
      integer, allocatable :: iregionr(:)
      real *8, allocatable :: srccoefs(:,:,:),srcvals(:,:,:)
c
      real *8 xyin1(2),xyin2(2),xyout(2),xsr(2),xsl(2)

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
c
c     true scattering problem: iscat = 1
c     total field : ifin = 1 otherwise, scattered field
c
      open(unit = 55,file = 'controls.dat',status='unknown')
      read(55,*) ifin
      read(55,*) iscat
      read(55,*) xyin1(1)
      read(55,*) xyin1(2)
      read(55,*) xyout(1)
      read(55,*) xyout(2)
      read(55,*) eps_discrete
      read(55,*) k
      read(55,*) npw
      read(55,*) chsmall
      read(55,*) ifpr_geom
      read(55,*) numit
      read(55,*) eps
      read(55,*) xlow
      read(55,*) ylow
      read(55,*) hhx
      read(55,*) hhy
      read(55,*) nnx
      read(55,*) nny
c
      call getspecs(nreg,nseg,zks,rns,ilside,irside,
     1               irefinel,irefiner,parsall,xsall,ysall,
     1               maxp,maxseg,maxreg,irbdry,nsegs,thetain,ipol)
c
      nover = 1
      allocate(chunks(2,k,maxc),ders(2,k,maxc),ders2(2,k,maxc))
      allocate(adjs(2,maxc),hs(maxc))
      allocate(iregionl(maxc),iregionr(maxc))
      allocate(srcvals(8,k,maxc),srccoefs(6,k,maxc),whts(k,maxc))
c
      call getdiscrete(nreg,nseg,zks,eps_discrete,k,nover,npw,rns,
     1            ilside,irside,iregionl,iregionr,irefinel,irefiner,
     1            alphas,betas,ifpr_geom,maxreg,maxc,
     1            chsmall,parsall,maxp,chunks,adjs,ders,ders2,hs,
     1            srcvals,srccoefs,whts,nch,maxseg,ipol)
c
ccc      if (2.ne.3) goto 222
      allocate(zrhs(k,nch,2),zsoln(k,nch,2))

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
           if (iscat.eq.0) then 
              call h2d_slp(xsr,8,srcvals(1,j,ich),0,dpars,1,zkr,0,
     1          ipars,zz)
              zrhs(j,ich,1) = zz*ar/zq
          
              call h2d_slp(xsl,8,srcvals(1,j,ich),0,dpars,1,zkl,0,
     1          ipars,zz)
             zrhs(j,ich,1) = zrhs(j,ich,1) - zz*al/zq
           else
              zrhs(j,ich,1) = 0.0d0
              if (iregionr(ich).eq.0) then
                 zrhs(j,ich,1) = zrhs(j,ich,1) 
     1                   -ar*uinfun(srcvals(1,j,ich),zkr,thetain)/zq
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
           endif
        enddo
      enddo
c
      call dielectric_solver(nreg,zks,alphas,betas,
     1               maxc,k,nch,srccoefs,srcvals,errs,numit,eps,
     1               iregionl,iregionr,zrhs,zsoln)
c
      if (iscat.eq.0) then 
c
c    test soln at an exterior point
c
         targ(1) = xyout(1)
         targ(2) = xyout(2)
c
         call  prin2(' xyin1 *',xyin1,2)
         call  prin2(' targ *',targ,2)
         call h2d_slp(xyin1,2,targ,0,dpars,1,zks(0),0,ipars,potex)
         call prin2('potex=*',potex,2)
c
         pot = 0
         pots = 0
         potd = 0
         do ich=1,nch
           do j=1,k
             zz = 0
             call h2d_dlp(srcvals(1,j,ich),2,targ,0,dpars,1,zks(0),0,
     1            ipars,zz)
             pot = pot + zsoln(j,ich,1)*zz*whts(j,ich)/betas(0)
             potd = potd + zsoln(j,ich,1)*zz*whts(j,ich)/betas(0)
             call h2d_slp(srcvals(1,j,ich),2,targ,0,dpars,1,zks(0),0,
     1            ipars,zz)
             pot = pot - zsoln(j,ich,2)*zz*whts(j,ich)/betas(0)
             pots = pots - zsoln(j,ich,2)*zz*whts(j,ich)/betas(0)
           enddo
         enddo
         call prin2('pots=*',pots,2)
         call prin2('potd=*',potd,2)
         call prin2('pot=*',pot,2)
         call prin2('abs error=*',abs(potex-pot),1)
         call prin2('rel error=*',abs(potex-pot)/abs(potex),1)
c
c    test soln at an interior point
c
         targ(1) = xyin1(1)
         targ(2) = xyin1(2)
c
         call h2d_slp(xyout,2,targ,0,dpars,1,zks(1),0,ipars,potex)
      
         pot = 0
         do ich=1,nch
           do j=1,k
             zz = 0
             call h2d_dlp(srcvals(1,j,ich),2,targ,0,dpars,1,zks(1),0,
     1            ipars,zz)
             pot = pot + zsoln(j,ich,1)*zz*whts(j,ich)/betas(1)
             call h2d_slp(srcvals(1,j,ich),2,targ,0,dpars,1,zks(1),0,
     1            ipars,zz)
             pot = pot - zsoln(j,ich,2)*zz*whts(j,ich)/betas(1)
           enddo
         enddo
         call prin2(' region 1 targ pot=*',pot,2)
         call prin2(' region 1 targ potex=*',potex,2)
         call prin2('abs error=*',abs(potex-pot),1)
         call prin2('rel error=*',abs(potex-pot)/abs(potex),1)
      endif
222   continue
c
      allocate(xtarg(nnx*nny))
      allocate(ytarg(nnx*nny))
      allocate(idt(nnx*nny))
      allocate(upotall(nnx*nny))
      allocate(potexall(nnx*nny))
c
      next = 0
      do j = 1,nny
      do i = 1,nnx
         next = next+1
         xtarg(next) = xlow + i*hhx
         ytarg(next) = ylow + j*hhy
      enddo
      enddo

      nnn = nch*k
      ntarget = nnx*nny
      call prinf(' calling domainflagsm nch is *',nch,1)
      t0 = second()
      call domainflagsm(k,nch,nnn,srccoefs,srcvals,
     1      iregionl,iregionr,nreg,ntarget,xtarg,ytarg,idt) 
      t1 = second()
      call prin2(' time for domainflagsm is *',t1-t0,1)
      open(unit = 61,file = 'domains.m',status='unknown')
      write(61,*) 'doms = ['
      do i = 1,ntarget
         write(61,*) xtarg(i), ytarg(i), idt(i)
      enddo
      write(61,*) '];'
      write(61,*) 'X = reshape(doms(:,1),',nnx,',',nny,');'
      write(61,*) 'Y = reshape(doms(:,2),',nnx,',',nny,');'
      write(61,*) 'Z = reshape(doms(:,3),',nnx,',',nny,');'
      write(61,*) 'mesh(X,Y,Z)'
      call prinf(' calling evalpot = *',n,1)
      t1 = second()
      call evalpotvol(k,nch,n,srccoefs,srcvals,zks,alphas,betas,
     1  zsoln(1,1,2),zsoln(1,1,1),nreg,idt,ntarget,xtarg,ytarg,upotall,
     2  ifin,uinfun,thetain)
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
            potexall(i) = 0.0d0
         endif
      enddo
ccc      call prin2( 'potexall is *',potexall,2*4)
c
c 
      open(unit = 19,file = 'pot.m',status='unknown')
      open(unit = 20,file = 'energy.m',status='unknown')
      if (iscat.eq.0) open(unit = 27,file = 'potex.m',status='unknown')
      write(19,*) 'rpot = ['
      write(20,*) 'edens = ['
      if (iscat.eq.0) write(27,*) 'epot = ['
      do i = 1,ntarget
         write(19,*) xtarg(i), ytarg(i), dreal(upotall(i))
         write(20,*) xtarg(i), ytarg(i), dreal(upotall(i))**2
         if (iscat.eq.0) 
     1              write(27,*) xtarg(i), ytarg(i), dreal(potexall(i))
      enddo
      t2 = second()
      call prin2(' time to eval at targs *',t2-t1,1)
      write(19,*) '];'
      write(19,*) 'X = reshape(rpot(:,1),',nnx,',',nny,');'
      write(19,*) 'Y = reshape(rpot(:,2),',nnx,',',nny,');'
      write(19,*) 'Z = reshape(rpot(:,3),',nnx,',',nny,');'
      write(19,*) 'mesh(X,Y,Z)'
      write(20,*) '];'
      write(20,*) 'X = reshape(edens(:,1),',nnx,',',nny,');'
      write(20,*) 'Y = reshape(edens(:,2),',nnx,',',nny,');'
      write(20,*) 'Z = reshape(edens(:,3),',nnx,',',nny,');'
      write(20,*) 'mesh(X,Y,Z)'
      if (iscat.eq.0) then
         write(27,*) '];'
         write(27,*) 'X = reshape(epot(:,1),',nnx,',',nny,');'
         write(27,*) 'Y = reshape(epot(:,2),',nnx,',',nny,');'
         write(27,*) 'Z = reshape(epot(:,3),',nnx,',',nny,');'
         write(27,*) 'mesh(X,Y,Z)'
      endif
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
c     incoming plane wave
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
c     normal derivative of incoming plane wave
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
