      implicit real *8 (a-h,o-z)
      real *8 xs(1000), ys(1000)
      real *8 xs2(10000), ys2(1000)
      real *8 xs3(1000), ys3(1000)
      real *8 dxdtp(1000), dydtp(1000)
      real *8 dxdt2p(1000), dydt2p(1000)
      real *8 kurvp(1000), kurv(1000)
      real *8 dxdt(1000), dydt(1000)
      real *8 dxdt2(1000), dydt2(1000)
      real *8 sval(1000)
      real *8 errx(1000)
      real *8 xcoeffs(100), ycoeffs(100)
      real *8 w1(100), w2(100000)
      real *8 pars(100)
      external fcurve1
c
c
      pi = 4*datan(1.0d0)
c
      rad = 2.0d0
      npts = 6
      frac = 0.125
      arc = frac*2*pi
      h = arc/(npts-1)
      write(6,*) ' frac is ',frac
      write(17,*) 2
      write(17,*) 9
      do jj = 1,8
         write(17,*) npts
         do i = 1,npts
            arg = (i-1)*h + (jj-1)*pi/4.0d0
            xs(i) = rad*dcos(arg) 
            ys(i) = rad*dsin(arg) 
            write(17,*) xs(i), ys(i)
         enddo
         if (jj.le.4) then 
            write(17,*) 1, 0
         else
            write(17,*) 2, 0
         endif
         il = 0
         ir = 0
         if (jj.eq.1) il = 1
         if (jj.eq.4) ir = 1
         if (jj.eq.5) il = 1
         if (jj.eq.8) ir = 1
         write(17,*) il, ir
      enddo
c
      jjj = 2
      write(17,*) jjj
      write(17,*) 2.0d0, 0.0d0
      write(17,*) -2.0d0, 0.0d0
      write(17,*) 2, 1
      write(17,*) 1, 1
      write(17,*) 5
      write(17,*) 1
      write(17,*) 2
      write(17,*) 3
      write(17,*) 4
      write(17,*) -9
      write(17,*) 5
      write(17,*) 5
      write(17,*) 6
      write(17,*) 7
      write(17,*) 8
      write(17,*) 9
     
      stop
      end
c
c
