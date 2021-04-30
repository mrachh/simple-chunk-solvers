      subroutine get_fmm_thresh(nds,ns,src,ndt,nt,trg,thresh)
c---------------------------------------
c  This subroutine returns the threshold used by an fmm
c  for ignoring self interactions for a given collection of sources
c  and targets.
c
c  The rotuine identifies the size of the smallest bounding cube
c  and sets the threshold to be 2**(-51)*(length of bounding cube edge)
c
c
c  Input arguments:
c  
c    - nds: integer
c        leading dimension for sources array (must be at least 3)
c    - ns: integer
c        number of sources
c    - src: real *8 (nds,ns)
c        source info array (src(1:3,i) should contain the x,y,z components)
c    - ndt: integer
c        leading dimension for targets array (must be at least 3)
c    - nt: integer
c        number of targets
c    - trg: real *8 (ndt,nt)c
c        targ info array (targ(1:3,i) should contain the x,y,z components)
c
c  Output arguments:
c  
c    - thresh: real *8
c        fmm threshold
c
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: nds,ns,ndt,nt
      real *8, intent(in) :: src(nds,ns),trg(ndt,nt)
      real *8, intent(out) :: thresh

      integer i
      real *8 xmin,xmax,ymin,ymax,bsizex,bsizey,bsize

      xmin = src(1,1)
      xmax = src(1,1)
      ymin = src(2,1)
      ymax = src(2,1)
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$REDUCTION(min:xmin,ymin)
C$OMP$REDUCTION(max:xmax,ymax)
C$OMP$PRIVATE(i)
      do i=1,ns
        if(src(1,i) .lt. xmin) xmin=src(1,i)
        if(src(1,i) .gt. xmax) xmax=src(1,i)
        if(src(2,i) .lt. ymin) ymin=src(2,i)
        if(src(2,i) .gt. ymax) ymax=src(2,i)
      enddo
C$OMP END PARALLEL DO    

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$REDUCTION(min:xmin,ymin)
C$OMP$REDUCTION(max:xmax,ymax)
C$OMP$PRIVATE(i)
      do i=1,nt
        if(trg(1,i) .lt. xmin) xmin=trg(1,i)
        if(trg(1,i) .gt. xmax) xmax=trg(1,i)
        if(trg(2,i) .lt. ymin) ymin=trg(2,i)
        if(trg(2,i) .gt. ymax) ymax=trg(2,i)
      enddo
C$OMP END PARALLEL DO    

      bsize = xmax-xmin
      bsizey = ymax-ymin
      if(bsizey.gt.bsize) bsize = bsizey

      thresh = bsize*1.0d-16

      return
      end
