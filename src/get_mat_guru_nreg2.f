c
      subroutine zgetmat_bdry_nreg(k,nch,n,srcinfo,srccoefs,
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
c     - srcinfo: double precision(8) 
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
c     - srcinfo: real *8 (8,k,nch)
c         source information
c         * srcinfo(1:2,j,ich) - x,y info at jth node on ichth chunk
c         * srcinfo(3:4,j,ich) - dxdt,dydt info at jth node on ichth
c           chunk
c         * srcinfo(5:6,j,ich) - d2xdt2,d2ydt2 info at jth node on ichth
c           chunk
c         * srcinfo(7:8,j,ich) - normal information at jth node on ichth
c           chunk
c     - srccoefs: real *8 (6,k,nch)
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
      real *8, intent(in) :: srcinfo(8,k,nch),srccoefs(6,k,nch)
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
      complex *16, allocatable :: xtmp(:,:),xtmp2(:),xtmp3(:,:)
      
      real *8, allocatable :: cm(:,:),rads(:),radtmp(:)
      real *8, allocatable :: srcover(:,:,:)
      integer, allocatable :: row_ptr(:),col_ind(:),col_ptr(:)
      integer, allocatable :: row_ind(:),iarr(:)
      real *8 rfac
      real *8, allocatable :: tsover(:),wover(:),uover(:,:),vover(:,:)
      real *8, allocatable :: ximat(:,:)
      complex *16, allocatable :: zximat(:,:)
      complex *16 zalpha,zbeta
      

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

      lwork = k*k*10+10*k + 100
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
c
c  find all close interactions
c
      allocate(cm(2,nch),rads(nch),radtmp(nch))
      
      call get_cm_rad(k,nch,srcinfo,cm,rads)
      rfac = 1.5d0
      do i=1,nch
        radtmp(i) = rads(i)*rfac
      enddo

cc      call prin2('cm=*',cm,24)
cc      call prin2('rads=*',rads,24)
cc      call prin2('radtmp=*',radtmp,24)
cc      call prinf('n=*',n,1)
cc      call prinf('nch=*',nch,1)

      
      call findnearmem(cm,nch,radtmp,8,srcinfo,n,nnz)
      call prinf('nnz=*',nnz,1)

      allocate(row_ptr(n+1),col_ptr(nch+1))
      allocate(row_ind(nnz),col_ind(nnz),iarr(nnz))

      call findnear(cm,nch,radtmp,8,srcinfo,n,row_ptr,col_ind)
      
      call rsc_to_csc(nch,n,nnz,row_ptr,col_ind,col_ptr,row_ind,
     1   iarr)

      kover = 20
      allocate(tsover(kover),wover(kover),uover(kover,kover))
      allocate(vover(kover,kover))

      itype = 2
      call legeexps(itype,kover,tsover,uover,vover,wover)

      allocate(ximat(kover,k),zximat(kover,k))
      call lematrin(k,kover,tsover,ximat,tstmp,work)
      do i=1,k
        do j=1,kover
          zximat(j,i) = ximat(j,i)
        enddo
      enddo

      call prin2('finished getting ximat*',i,0)

      allocate(srcover(8,kover,nch))
c
c  oversample the geometry information
c
      alpha = 1.0d0
      beta = 0.0d0
      do ich=1,nch
        call dgemm('n','t',8,kover,k,alpha,srcinfo(1,1,ich),
     1   8,ximat,kover,beta,srcover(1,1,ich),8)
      enddo

      call prinf('finished oversamling*',i,0)

      allocate(xtmp3(k,kover))
c
c  end of oversampling geometry
c

c
c  generate matrix entries corresponding to smooth quadrature
c  where self is corrected via ggq
c 
      zalpha = 1.0d0
      zbeta = 0.0d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(isch,ismat,itch,itmat)
C$OMP$PRIVATE(i,j,xtmp,xtmp2,xtmp3)
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

          if(isch.eq.itch) then
            call zself_buildmat(k,srccoefs(1,1,isch),
     1         srcinfo(1,1,itch),fker,ndt,ndd,dpars,ndz,zpars,
     2         ndi,ipars,xs0,whts0,ainterp,xtmp)
          else
            do j=1,kover
              dst = sqrt(srcover(3,j,isch)**2 + srcover(4,j,isch)**2)*
     1           wover(j)
              do i=1,k
                call fker(srcover(1,j,isch),8,srcinfo(1,i,itch),
     1            ndd,dpars,ndz,zpars,ndi,ipars,xtmp3(i,j))
                xtmp3(i,j) = xtmp3(i,j)*dst
              enddo
            enddo
            call zgemm('n','n',k,k,kover,zalpha,xtmp3,k,zximat,
     1        kover,zbeta,xtmp,k)
          endif
          call zinsertmat(k,k,xtmp,itmat,ismat,n,n,xmat)
        enddo
      enddo
C$OMP END PARALLEL DO     
 
      call prinf('finished n2 work*',i,0)

c
c  now correct the near part of the matrices
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(isch,istart,ii,itarg,itch)
C$OMP$PRIVATE(j,xtmp2,i)
      do isch=1,nch
        istart = (isch-1)*k
        nn = col_ptr(isch+1)-col_ptr(isch)
cc        call prinf('rel targs=*',row_ind(col_ptr(isch)),nn)
        do ii=col_ptr(isch),col_ptr(isch+1)-1
c
c  identify chunk number and node number of given target
c
          itarg = row_ind(ii)
          itch = ceiling((itarg+0.0d0)/k)
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
          if(isch.eq.itch) goto 1000
          j = itarg - (itch-1)*k
          xtmp2 = 0
          call zadapquad(k,srccoefs(1,1,isch),srcinfo(1,j,itch),
     1       fker,8,ndd,dpars,ndz,zpars,ndi,ipars,umat,m,tsquad,
     2       wquad,xtmp2)
          do i=1,k
            xmat(itarg,istart+i) = xtmp2(i)
          enddo
 1000     continue          
        enddo
      enddo
C$OMP END PARALLEL DO      

      return
      end
c
c
c
