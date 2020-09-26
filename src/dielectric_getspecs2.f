c
c     getgeometry for dielectric interface problem
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
      subroutine getspecs(nreg,nseg,zks,rns,ilside,irside,
     1               irefinel,irefiner,parsall,xsall,ysall,
     1               maxp,maxseg,maxreg,thetain,ipol)
ccc     1               maxp,maxseg,maxreg,irbdry,nsegs,thetain,ipol)

      implicit real *8 (a-h,o-z)

      parameter (maxc = 10000)
 
ccc      integer irbdry(maxseg,maxreg)
      integer ilside(maxseg)
      integer irside(maxseg)
      integer irefinel(maxseg)
      integer irefiner(maxseg)
ccc      integer nsegs(maxreg)
      complex *16 zks(0:maxreg),eye
      complex *16 rns(0:maxreg)
      real *8  xsall(maxp,maxseg)
      real *8  ysall(maxp,maxseg)
      real *8  parsall(2*maxp+3,*)
      real *8  xcoeffs(1000)
      real *8  ycoeffs(1000)
      real *8  kurv(1000)
      real *8  sval(1000)
      real *8 xyin1(2),xyin2(2),xyout(2),xsr(2),xsl(2)
      real *8 w1(100), w2(100000)
      data eye/(0.0d0,1.0d0)/
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
            read(5,*) xsall(j,ii), ysall(j,ii)
         enddo
         call projectonpoly(xsall(1,ii),ysall(1,ii),npts,sval,
     1         xcoeffs,ycoeffs,w1,w2)
         parsall(1,ii) = npts
         parsall(2,ii) = sval(1)
         parsall(3,ii) = sval(npts)
c
         do j = 1,npts
            parsall(2*j+2,ii) = xcoeffs(j)
            parsall(2*j+3,ii) = ycoeffs(j)
         enddo
c
         ikurv = 0
         if (ikurv.eq.1) then
            npts2 = 2*npts
            hh = (sval(npts)-sval(1))/(npts2-1)
            do i = 1,npts2
               sv = hh*(i-1)
               call fcurve1(sv,parsall(1,ii),xx,yy,dxdt,dydt,
     1           dxdt2,dydt2)
               rn = (dxdt**2 + dydt**2)**(1.5d0)
               kurv(i) = (dxdt*dydt2 - dydt*dxdt2)/rn
            enddo
         endif
ccc         call prin2(' kurv is *',kurv,npts2)
c
         read(5,*) ilside(ii),irside(ii)
         read(5,*) irefinel(ii),irefiner(ii)
         write(6,*) 'ii is ',ii
         write(6,*) 'ilside, irside is ',ilside(ii), irside(ii)
         write(6,*) 'irefl, irefr is ',irefinel(ii), irefiner(ii)
      enddo
ccc      do ireg = 1,nreg
ccc         read(5,*) nsegs(ireg)
ccc         do ii = 1,nsegs(ireg)
ccc            read(5,*) irbdry(ii,ireg)
ccc         enddo
ccc         write(6,*) 'irbdry is ',(irbdry(ii,ireg),ii=1,nsegs(ireg))
ccc      enddo
      read(5,*) xxx,yyy
      zks(0) = xxx+ eye*yyy
      write(6,*) ' zk0 is ',zks(0)
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
c     theta = 0  -> normal,   0 < theta < pi/2 is rightgoing
c                         -pi/2 < theta < 0 is leftgoing
c
c
      return
      end
c-----------------------------------------------------------
c

