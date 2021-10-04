



      subroutine zgetmat_bdryc_nreg(k,nch,n,srcinfo,srccoefs,
     1    fker,ndt,ndd,dpars,zks,alphas,betas,nregion,
     2    iregionl,iregionr,ibl,ndi,ipars,xmat)
c
c   this subroutine generates the matrix corresponding to the
c   discretization of a weakly singular integral operator
c   whose kernel is given by fker
c
c   The expected calling sequence of fker is
c      fker(srcinfo,ndt,targinfo,ndd,dpars,ndz,zpars,ndi,ipars,val)
c
c   where 
c     - srcinfo: complex double precision
c         is the source info
c     - ndt: integer
c         <=8 is the expected
c       number of parameters to be used in the target info,
c     - targinfo: double precision(ndt) 
c         is the target info
c     - ndd: integer
c         number of double precision parameters
c     - dpars: double precision(ndd) 
c         the double precision parameters 
c     - zks: complex *16(3)
c        zk(i): Helmholtz parameter in region i
c     -  alphas: complex *16(3)
c       alphas(i): alpha parameter in region i
c     - betas: complex *16(3)
c        betas(i): beta parameter in region i
c     - nregion: integer        - number of regions
c     - iregionl: integer (nch)  - left domain for chunk
c     - iregionr: integer (nch)  - right domain for chunk
c     - ibl: integer   
c             ibl = 1 => 11 block  of system matrix
c             ibl = 2 => 12 block  of system matrix
c             ibl = 3 => 22 block  of system matrix
c             ibl = 4 => 21 block  of system matrix
c     - ndi: integer
c         number of integer parameters
c     - ipars: integer(ndi) 
c         integer parameters 
c
c   other input arguments:
c     - k: integer
c         order of discretization
c     - nch: integer
c         number of chunks
c     - n: integer
c        total number of points on the boundary = k*nch,
c        also leading dimension of output matrix
c     - srcinfo: complex *16 (8,k,nch)
c         source information
c         * srcinfo(1:2,j,ich) - x,y info at jth node on ichth chunk
c         * srcinfo(3:4,j,ich) - dxdt,dydt info at jth node on ichth
c           chunk
c         * srcinfo(5:6,j,ich) - d2xdt2,d2ydt2 info at jth node on ichth
c           chunk
c         * srcinfo(7:8,j,ich) - normal information at jth node on ichth
c           chunk
c     - srccoefs: complex *16 (6,k,nch)
c         Legendre expansion coeffs of the first 6 arrays in srcinfo
c         array
c   output parameters:
c     - xmat: complex *16(n,n)
c         discretized matrix corresponding to fker
c
c  This subroutine uses generalized gaussian quadratures
c  for self, and adaptive integration for the rest of the chunks
c
c-----

      implicit real *8 (a-h,o-z)
      integer, intent(in) :: k,nch,n,nregion
      complex *16, intent(in) :: srcinfo(8,k,nch),srccoefs(6,k,nch)
      integer, intent(in) :: ndt,ndd,ndi
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zks(0:nregion)
      complex *16, intent(in) :: alphas(0:nregion)
      complex *16, intent(in) :: betas(0:nregion)
      integer, intent(in) :: iregionl(nch)
      integer, intent(in) :: iregionr(nch)
      integer, intent(in) :: ipars(ndi)

      complex *16, intent(out) :: xmat(n,n)
      external fker

c
c    breemer \self quad variables
c
      integer nquad0,nquad1
      real *8, allocatable :: xs1(:),whts1(:),xs0(:,:),whts0(:,:)
      real *8, allocatable :: ainterp(:,:,:)
      real *8, allocatable :: ts0(:),w0(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: tsquad(:),wquad(:)
      real *8, allocatable :: tstmp(:)
      real *8, allocatable :: work(:)
      complex *16, allocatable :: xtmp(:,:),xtmp2(:)

c
c       temporary variables 
c
      integer i,j,l,ndz
      integer ich,jch,isch,itch
      integer ismat,itmat
      integer lwork
      complex *16 zpars(6),ar,al,br,bl,zq

      call get_quads_info(k, nquad1, nquad0)
      allocate(xs1(nquad1), whts1(nquad1))
      allocate(xs0(k,k), whts0(k,k))
      call get_quads(k, nquad1, xs1, whts1, nquad0, xs0, whts0)

      itype = 2
      allocate(ts0(k),w0(k),umat(k,k),vmat(k,k))
      call legeexps(itype,k,ts0,umat,vmat,w0)

      allocate(tstmp(k+10))
      
      allocate(xtmp(k,k),xtmp2(k))
      allocate(ainterp(k,k,k))

      lwork = k*k*10+10
      allocate(work(lwork))
      do i=1,k
        call lematrin(k,k,xs0(1,i),ainterp(1,1,i),tstmp,work)
      enddo

      m = 20
      itype = 1
      allocate(tsquad(m),wquad(m))
      call legeexps(itype,m,tsquad,utmp,vtmp,wquad) 
c
      ndz = 6

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(isch,ismat,itch,itmat)
C$OMP$PRIVATE(i,j,xtmp,xtmp2)
      do isch=1,nch
        ismat = (isch-1)*k + 1
        do itch=1,nch
          itmat = (itch-1)*k + 1
          do i=1,k
            do j=1,k
              xtmp(j,i) = 0
            enddo
          enddo
          zpars(1) = zks(iregionr(itch))
          zpars(2) = zks(iregionl(itch))
          ar = alphas(iregionr(itch))
          al = alphas(iregionl(itch))
          br = betas(iregionr(itch))
          bl = betas(iregionl(itch))
          zq = 0.5*(ar/br+al/bl)
          if (ibl.eq.1) then
             zpars(3) = 0.0d0
             zpars(4) = 0.0d0
             zpars(5) = (ar/br)/zq
             zpars(6) = -(al/bl)/zq
          else if (ibl.eq.2) then
             zpars(3) = -(ar/br)/zq
             zpars(4) = (al/bl)/zq
             zpars(5) = 0.0d0
             zpars(6) = 0.0d0
          else if (ibl.eq.3) then
             zpars(3) = -1.0d0
             zpars(4) = 1.0d0
             zpars(5) = 0.0d0
             zpars(6) = 0.0d0
          else
             zpars(3) = 0.0d0
             zpars(4) = 0.0d0
             zpars(5) = 1.0d0
             zpars(6) = -1.0d0
          endif
c
          if(isch.eq.itch) then
            call zself_buildmatc(k,srccoefs(1,1,isch),
     1         srcinfo(1,1,itch),fker,ndt,ndd,dpars,ndz,zpars,
     2         ndi,ipars,xs0,whts0,ainterp,xtmp)
          else
            do j=1,k
              xtmp2 = 0
              call zadapquadc(k,srccoefs(1,1,isch),srcinfo(1,j,itch),
     1          fker,8,ndd,dpars,ndz,zpars,ndi,ipars,umat,m,tsquad,
     2          wquad,xtmp2)
              do i=1,k
                xtmp(j,i) = xtmp2(i)
              enddo
            enddo
          endif
          call zinsertmat(k,k,xtmp,itmat,ismat,n,n,xmat)
        enddo
      enddo
C$OMP END PARALLEL DO      

      return
      end
c
c
c
