c
c-----------------------------------------------------------
c
cc Copyright (C) 2020: Leslie Greengard 
c
c     parsall are x1,x2,y1,y2 corrds of segmnt endpoints.
c 
c     idomain(j,k) = ordered segments for domain k
c     if j is pos, it means segment is positively 
c     oriented and points added to boundary discretization in that order.
c     if j is neg, it means segment is negatively 
c     oriented and points need to be added to boundary discretization in 
c     reverse order.
c
c     1) For region = 1,nreg: flag points i in region j. Set idt(i) = j.
c     2) All other targets are in exterior -> idt(i) = 0.
c
c   This file takes a set of curves Gamma1,..,GammaN 
c     defining domains  D1...DN and their exterior D0
c
c     and a collection of targets t1,...tM  and assigns an
c     ingteger flag ID(i) = j if the target is in domain j
c     (defined by polygonal curves (xs,ys, i = 1,ns)
c
c     Uses the fact that \int_Gamma_i dG/dn = 1 if target is in
c     D_i and zero otherwise.
c
c      |-------|        |-------|
c      | D1    |        | D2    |  D0
c      |       |        |       |
c      |-------|        |-------|
c
c
      subroutine domainflag(nsmax,nreg,nseg,nsegs,parsall,
     1             irbdry,xt,yt,nnx,nny,idt)
c
      implicit real *8 (a-h,o-z)
      integer irbdry(nsmax,nreg)
      integer nsegs(nreg)
      integer  nb,ntot,nstart,nseg
      integer  idt(nnx*nny)
      real *8 parsall(4,*)
      real *8 xs(1000),ys(1000)
      real *8 xt(nnx*nny),yt(nnx*nny)
      real *8 val
c
ccc      write(6,*) 'nsmax',nsmax
ccc      write(6,*) 'nnx',nnx
ccc      write(6,*) 'nny',nny
      ntarg = nnx*nny
ccc      write(6,*) 'ntarg',ntarg
      do ii = 1,ntarg
        idt(ii) = 0
      enddo
cc
ccc      call prinf(' ns is *',ns,nb)
c
ccc      write(6,*) 'nsmax',nsmax
ccc      write(6,*) 'nreg',nreg
ccc      write(6,*) 'nseg',nseg
ccc      write(6,*) 'nsegs',(nsegs(i),i=1,nreg)
ccc      do iseg = 1,nseg
ccc         write(6,*) 'parsall',(parsall(i,iseg),i=1,4)
ccc         write(13,*) 'parsall',(parsall(i,iseg),i=1,4)
ccc      enddo
ccc      do ireg = 1,nreg
ccc         write(6,*) 'irbdry',(irbdry(i,ireg),i=1,nsegs(ireg))
ccc         write(13,*) 'irbdry',(irbdry(i,ireg),i=1,nsegs(ireg))
ccc      enddo
c
      do ireg = 1,nreg
         inext = 0
         do ii = 1,nsegs(ireg)
            iseg = irbdry(ii,ireg)
            if (iseg.gt.0) then
               inext = inext+1
               xs(inext) = parsall(1,iseg)
               ys(inext) = parsall(3,iseg)
               inext = inext+1
               xs(inext) = 0.5*(parsall(1,iseg)+parsall(2,iseg))
               ys(inext) = 0.5*(parsall(3,iseg)+parsall(4,iseg))
            else
               inext = inext+1
               xs(inext) = parsall(2,-iseg)
               ys(inext) = parsall(4,-iseg)
               inext = inext+1
               xs(inext) = 0.5*(parsall(1,-iseg)+parsall(2,-iseg))
               ys(inext) = 0.5*(parsall(3,-iseg)+parsall(4,-iseg))
            endif
         enddo
ccc         write(6,*) 'xs',(xs(ii),ii=1,inext)
ccc         write(6,*) 'ys',(ys(ii),ii=1,inext)
ccc         write(6,*) 'ireg,inext',ireg,inext
ccc         write(6,*) 'xs',(xs(ii),ii=1,inext)
ccc         write(6,*) 'ys',(ys(ii),ii=1,inext)
ccc         xx = 0.9d0
ccc         yy = -0.0005d0
ccc         call dlpcomp(xs,ys,inext,
ccc     1           xx,yy,ival)
ccc         write(6,*) 'at 0,.5 ival',ival
         do j=1,ntarg
            call dlpcomp(xs,ys,inext,
     1           xt(j),yt(j),ival)
            idt(j) = idt(j) + ival*ireg
         enddo
      enddo
      iwrite = 0
      if (iwrite.eq.1) then
         open(unit = 19,file = 'domainpot.m',status='unknown');
         write(19,*) 'rpot = ['
         do i = 1,ntarg
            write(19,*) xt(i), yt(i), idt(i)
         enddo
         write(19,*) '];'
         write(19,*) 'X = reshape(rpot(:,1),',nnx,',',nny,')'
         write(19,*) 'Y = reshape(rpot(:,2),',nnx,',',nny,')'
         write(19,*) 'Z = reshape(rpot(:,3),',nnx,',',nny,')'
         write(19,*) 'mesh(X,Y,Z)'
      endif
c
      return
      end
c
c
c
      subroutine dlpcomp(xs,ys,ns,xt,yt,ival)
c
c     assumes curve is pos oriented:
c
c
c           p -------- Q_{i+1}
c             -  
c               -
c                 -
c                   - Q_i
c
c     v = (Q_i - P)/||Q_i-P||  = (v1,v2)
c     w = (Q_{i+1} - P)/||Q_{i+1}-P|| = (w1,w2)
c     n = v \times w = (0,0,v1w2-v2w1)
c     dtheta = asin(v1w2-v2w1)
c
c
c
      implicit real *8 (a-h,o-z)
      integer  ns,ival
      real *8 xs(ns),ys(ns)
      real *8 xt,yt
      real *8 val,pi
      complex *16 z1,z2,z3
c
      pi = 4.0d0*datan(1.0d0)
      val = 0.0d0
      do ii=1,ns-1
         v1 = xs(ii)-xt
         v2 = ys(ii)-yt
         rr = sqrt(v1*v1+v2*v2)
         z1 = dcmplx(v1,v2)
         v1 = v1/rr
         v2 = v2/rr
         w1 = xs(ii+1)-xt
         w2 = ys(ii+1)-yt
         rr = sqrt(w1*w1+w2*w2)
         z2 = dcmplx(w1,w2)
         w1 = w1/rr
         w2 = w2/rr
ccc         write(6,*) 'v1,v2 ',v1,v2
ccc         write(6,*) 'w1,w2 ',w1,w2
ccc         write(6,*) 'cross ',v1*w2-v2*w1
         z3 = z2/z1
ccc         val = val + asin(v1*w2-v2*w1)
         x = dreal(z3)
         y = dimag(z3)
         val = val + atan2(y,x)
ccc         write(6,*) 'val ',val
      enddo
      ival = nint(val/(2.0d0*pi))
c
      return
      end
c
c
