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
c    11) polarization (ipol = 0 -> betas=1 , ipol=1 -> beta = 1/rn^2)
c
      subroutine dielectric_solver(nreg,zks,alphas,betas,
     1               maxc,k,nch,srccoefs,srcvals,errs,numit,eps,
     1               iregionl,iregionr,zrhs,zsoln)

      implicit none
c
      integer nreg,maxc,k,nch,numit,niter,nsys,ntot
      integer iregionl(maxc)
      integer iregionr(maxc)
      real *8 srccoefs(6,k,nch),srcvals(8,k,nch)
      real *8 errs(*),rres,eps,t1,t2
      complex *16 alphas(0:nreg),betas(0:nreg)
      complex *16 zks(0:nreg)
      complex *16 zrhs(k,nch,2)
      complex *16 zsoln(k,nch,2)
      complex *16, allocatable :: xmat(:,:)
c
      complex *16 zz,pot,potex,zpars(6),zid
      complex *16 pots,potd
      complex *16 zq
c
      ntot = k*nch
      nsys = 2*ntot
      call prinf('nsys =*',nsys,1)
      allocate(xmat(nsys,nsys))
c      
      call cpu_time(t1)
C$       t1 = omp_get_wtime()
      call helm2d_dielec_mat_nreg(k,nch,ntot,srcvals,srccoefs,zks,
     1     alphas,betas,nreg,iregionl,iregionr,xmat)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      call prin2('matrix generation time=*',t2-t1,1)
      zid = 1.0d0
      niter = 0
c
ccc      call prin2(' zrhs is *',zrhs,2*n*2)
c
      call cpu_time(t1)
C$       t1 = omp_get_wtime()
      call zgmres_solver(nsys,xmat,zid,zrhs,numit,eps,niter,
     1   errs,rres,zsoln)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      call prin2('matrix solve time=*',t2-t1,1)
ccc      call prin2(' sol is *',zsoln,2*n*2)
c
      return
      end
c-----------------------------------------------------------
c
