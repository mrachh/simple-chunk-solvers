c
c     getgeometry for dirichlet scattering
c
c     A simple example is:
c
c     (-1,0)_/\__/\__/\______ (1,0)    sin series
c          |                 |
c          |                 |
c          |                 |
c   (-1,-D)|_________________| (1.-D)
c
c
c    The input file states: 
c
c    angle of orientation of incoming plane wave
c
      subroutine getspecs(iscat,itot,depth,zk,parsall,maxp,
     1               norder,thetain)

      implicit real *8 (a-h,o-z)
      complex *16 zk,eye
      real *8  parsall(maxp+3,4)
c
c     read in vertices from file
c     line 1   nseg
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
      open (unit=5,file='daria_geo.dat',status='unknown')
c
      pi = 4*datan(1.0d0)
      read(5,*) iscat
      read(5,*) itot
      read(5,*) depth
      read(5,*) norder
      read(5,*) xx
      zk = dcmplx(xx,0.0d0)
      write(6,*) 'D is ',depth
      parsall(1,1) =  0.0d0
      parsall(3,1) =  0.0d0
      parsall(2,1) =  0.0d0
      parsall(4,1) = -depth
c
      parsall(1,2) = 0.0d0
      parsall(3,2) = -depth
      parsall(2,2) = pi
      parsall(4,2) = -depth
c
      parsall(1,3) = pi
      parsall(3,3) = -depth
      parsall(2,3) = pi
      parsall(4,3) = 0.0d0
c
      parsall(1,4) = pi
      parsall(2,4) = 0.0d0
      parsall(3,4) = norder
      do i = 1,norder
         read(5,*) parsall(3+i,4) 
      enddo
      read(5,*) thetain
      return
      end
c-----------------------------------------------------------
c    t:   1 -> -1    => x, y is sum_{l=1}^porder a_l sin( pi l x)

