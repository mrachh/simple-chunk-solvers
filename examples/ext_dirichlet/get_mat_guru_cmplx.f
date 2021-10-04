



      subroutine zgetmat_bdryc(k,nch,n,srcinfo,srccoefs,
     1    fker,ndt,ndd,dpars,ndz,zpars,ndi,ipars,xmat)
c
c   this subroutine generates the matrix corresponding to the
c   discretization of a weakly singular integral operator
c   whose kernel is given by fker
c
c   The expected calling sequence of fker is
c      fker(srcinfo,ndt,targinfo,ndd,dpars,ndz,zpars,ndi,ipars,val)
c
c   where 
c     - srcinfo: double complex
c         is the source info
c     - ndt: integer
c         <=8 is the expected
c       number of parameters to be used in the target info,
c     - targinfo: double complex
c         is the target info
c     - ndd: integer
c         number of double precision parameters
c     - dpars: double precision(ndd) 
c         the double precision parameters 
c     - ndz: integer
c         number of complex parameters
c     - zpars: double precision(ndz)
c         complex *16 parameters
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
      integer, intent(in) :: k,nch,n
      complex *16, intent(in) :: srcinfo(8,k,nch),srccoefs(6,k,nch)
      integer, intent(in) :: ndt,ndd,ndi,ndz
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
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
      real *8, allocatable :: rsrcinfo(:,:,:)
ccc      real *8, allocatable :: rsrccoefs(:,:,:)
      
      real *8, allocatable :: cm(:,:),rads(:),radtmp(:)
      complex *16, allocatable :: srcover(:,:,:)
      integer, allocatable :: row_ptr(:),col_ind(:),col_ptr(:)
      integer, allocatable :: row_ind(:),iarr(:)
      real *8 rfac
      real *8, allocatable :: tsover(:),wover(:),uover(:,:),vover(:,:)
      real *8, allocatable :: ximat(:,:)
      complex *16, allocatable :: zximat(:,:)
      complex *16 zalpha,zbeta,dst
      

c
c       temporary variables 
c
      integer i,j,l
      integer ich,jch,isch,itch
      integer ismat,itmat
      integer lwork

      call get_quads_info(k, nquad1, nquad0)
      allocate(rsrcinfo(8,k,nch))
ccc      allocate(rsrccoeffs(6,k,nch))
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

      do i = 1,nch
         do kk = 1,k
            do j = 1,8
               rsrcinfo(j,kk,i) = real(srcinfo(j,kk,i))
            enddo
         enddo
      enddo
c
c  find all close interactions
c
      allocate(cm(2,nch),rads(nch),radtmp(nch))
      
      call get_cm_rad(k,nch,rsrcinfo,cm,rads)
      rfac = 2.0d0
      do i=1,nch
        radtmp(i) = rads(i)*rfac
      enddo

cc      call prin2('cm=*',cm,24)
cc      call prin2('rads=*',rads,24)
cc      call prin2('radtmp=*',radtmp,24)
cc      call prinf('n=*',n,1)
cc      call prinf('nch=*',nch,1)
      
      call findnearmem(cm,nch,radtmp,8,rsrcinfo,n,nnz)

      allocate(row_ptr(n+1),col_ptr(nch+1))
      allocate(row_ind(nnz),col_ind(nnz),iarr(nnz))

      call findnear(cm,nch,radtmp,8,rsrcinfo,n,row_ptr,col_ind)
      
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


      allocate(srcover(8,kover,nch))
c
c  oversample the geometry information
c
      alpha = 1.0d0
      beta = 0.0d0
      do ich=1,nch
        do j=1,kover
          do i=1,8
            srcover(i,j,ich) = 0
            do l=1,k
              srcover(i,j,ich) = srcover(i,j,ich) + 
     1           ximat(j,l)*srcinfo(i,l,ich)
            enddo
          enddo
        enddo
      enddo


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

          if(isch.eq.itch) then
            call zself_buildmatc(k,srccoefs(1,1,isch),
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
            do j=1,k
              do i=1,k
                xtmp(i,j) = 0
                do l=1,kover
                  xtmp(i,j) = xtmp(i,j) + xtmp3(i,l)*ximat(l,j) 
                enddo
              enddo
            enddo
c            call zgemm('n','n',k,k,kover,zalpha,xtmp3,k,zximat,
c     1        kover,zbeta,xtmp,k)
          endif
          call zinsertmat(k,k,xtmp,itmat,ismat,n,n,xmat)
        enddo
      enddo
C$OMP END PARALLEL DO     
 

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
          if(isch.eq.itch) goto 1000
          j = itarg - (itch-1)*k
          xtmp2 = 0
          call zadapquadc(k,srccoefs(1,1,isch),srcinfo(1,j,itch),
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
c
c
      subroutine zself_buildmatc(k,srccoefs,srcvals,fker,ndt,ndd,
     1   dpars,ndz,zpars,ndi,ipars,xs0,whts0,ainterp,xmat)
      implicit real *8 (a-h,o-z)
      integer k
      complex *16 srccoefs(6,k),srcvals(8,k)
      integer ndt,ndd,ndz,ndi
      real *8 dpars(ndd)
      complex *16 zpars(ndz),dst
      integer ipars(ndi)
      real *8 xs0(k,k),whts0(k,k)
      complex *16 xmat(k,k)
      real *8 ainterp(k,k,k)
      complex *16 wtmp2(k),wtmp3(k)
      complex *16 val
      real *8 pols(k)
      complex *16 srcinfo(8)

      real *8 tsrc
      integer i,j,l
      real *8 alpha,beta
      external fker

      alpha = 1.0d0
      beta = 0.0d0

      do i=1,k
        do j=1,k
          tsrc = xs0(j,i)
          call legepols(tsrc,k-1,pols)
          do m=1,8
            srcinfo(m) = 0
          enddo
          do l=1,k
            do m=1,6
              srcinfo(m) = srcinfo(m) + srccoefs(m,l)*pols(l)
            enddo
          enddo
          dst = sqrt(srcinfo(3)**2 + srcinfo(4)**2)
          srcinfo(7) = srcinfo(4)/dst
          srcinfo(8) = -srcinfo(3)/dst

          call fker(srcinfo,8,srcvals(1,i),ndd,dpars,ndz,zpars,ndi,
     1      ipars,wtmp2(j))
          wtmp2(j) = wtmp2(j)*whts0(j,i)*dst
        enddo

        do j=1,k
          wtmp3(j) = 0
          do l=1,k
            wtmp3(j) = wtmp3(j) + wtmp2(l)*ainterp(l,j,i)
          enddo
        enddo

        do j=1,k
          xmat(i,j) = wtmp3(j)
        enddo
      enddo
      
      return
      end
c
c
c
c
c
      subroutine zadapquadc(k,srccoefs,targ,fker,ndt,ndd,dpars,ndz,
     1   zpars,ndi,ipars,umat,m,ts0,w0,xmat)
      implicit real *8 (a-h,o-z)
      complex *16 srccoefs(6,k),targ(ndt)
      integer ndt,ndd,ndz,ndi
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)
      real *8 umat(k,k),ts0(m),w0(m)
      complex *16 xmat(k)
      complex *16 rints(k)
      real *8, allocatable :: stack(:,:)
      complex *16, allocatable :: vals(:,:)
      complex *16 value2(k),value3(k)
      external fker

      a = -1.0d0
      b = 1.0d0

      maxrec = 0
      numit = 0
      
      maxdepth = 200
      allocate(stack(2,maxdepth),vals(k,maxdepth))

      nnmax = 100000
      eps = 1.0d-12
      do j=1,maxdepth
        stack(1,j) = 0
        stack(2,j) = 0
        vals(1:k,j) = 0
      enddo



      call adinrecmc(ier,stack,a,b,k,srccoefs,targ,fker,
     1   ndt,ndd,dpars,ndz,zpars,ndi,ipars,m,ts0,w0,vals,nnmax,
     2   eps,rints,maxdepth,maxrec,numit,value2,value3)

c
c      convert rints to matrix entries
c
      call zrmatmatt(1,k,rints,k,umat,xmat)

      return
      end
c
c
c
c

      subroutine adinrecmc(ier,stack,a,b,k,srccoefs,targ,fker,
     1   ndt,ndd,dpars,ndz,zpars,ndi,ipars,m,ts0,w0,vals,nnmax,
     2   eps,rints,maxdepth,maxrec,numit,value2,value3)
      implicit real *8 (a-h,o-z)
      real *8 stack(2,maxdepth),a,b
      integer k
      complex *16 srccoefs(6,k),targ(ndt)
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)
      real *8 ts0(m),w0(m)
      complex *16 vals(k,maxdepth),rints(k),value2(k),value3(k)
      external fker
      
      


      ier = 0
      stack(1,1) = a
      stack(2,1) = b

      call oneintmuc(a,b,k,srccoefs,targ,fker,ndt,ndd,dpars,ndz,
     1   zpars,ndi,ipars,m,ts0,w0,vals(1,1))
      
c 
c       recursively integrate the thing
c 
      j=1
      do i=1,k
        rints(i)=0
      enddo
c 
      ier=0
      maxrec=0
      do 3000 i=1,nnmax
        numint=i
        if(j .gt. maxrec) maxrec=j
c 
c       subdivide the current subinterval
c 
         c=(stack(1,j)+stack(2,j))/2
         call oneintmuc(stack(1,j),c,k,srccoefs,targ,fker,ndt,ndd,
     1      dpars,ndz,zpars,ndi,ipars,m,ts0,w0,value2)
         call oneintmuc(c,stack(2,j),k,srccoefs,targ,fker,ndt,ndd,
     1      dpars,ndz,zpars,ndi,ipars,m,ts0,w0,value3)
c 
          dd=0
          do jj=1,k
            ddd=abs(value2(jj)+value3(jj)-vals(jj,j) )
            if(ddd .gt. dd) dd=ddd
          enddo
c 
          ifdone=0
          if(dd .le. eps) ifdone=1
c 
c       if the function on this subinterval has been
c       integrated with sufficient accuracy - add the
c       value to that of the global integral and move up
c       in the stack
c 
          if(ifdone  .eq. 0) goto 2000
c 
          do jj=1,k
            rints(jj)=rints(jj)+value2(jj)+value3(jj)
          enddo
          j=j-1
c 
c        if the whole thing has been integrated - return
c 
          if(j .eq. 0) return
          goto 3000
 2000 continue
c 
c       if the function on this subinterval has not been
c       integrated with sufficient accuracy - move
c       down the stack
c 
        stack(1,j+1)=stack(1,j)
        stack(2,j+1)=(stack(1,j)+stack(2,j))/2
        vals(1:k,j+1) = value2(1:k)
        vals(1:k,j) = value3(1:k)
c 
        stack(1,j)=(stack(1,j)+stack(2,j))/2
c 
        j=j+1
c 
c       if the depth of the recursion has become excessive - bomb
c 
        if(j .lt. maxdepth) goto 3000
        ier=8
        return
 3000 continue
      ier=16

      return
      end
c
c
c
c
c
      subroutine oneintmuc(a,b,k,srccoefs,targ,fker,ndt,ndd,dpars,
     1   ndz,zpars,ndi,ipars,m,t,w,vals)
      implicit real *8 (a-h,o-z)
      real *8 a,b
      complex *16 srccoefs(6,k),targ(ndt)
      complex *16 srctmp(8),dst
      real *8 dpars(ndd),pols(k)
      complex *16 zpars(ndz)
      integer ipars(ndi)
      real *8 t(m),w(m)
      complex *16 vals(k),val
      external fker

      
      do i=1,k
        vals(i) = 0
      enddo

      u = (b-a)/2
      v = (b+a)/2
      alpha = 1.0d0
      beta = 0.0d0

      do i=1,m
        tt = u*t(i)+v
        call legepols(tt,k-1,pols) 

        do j=1,8
          srctmp(j) = 0
        enddo
        do l=1,k
          do j=1,6
            srctmp(j) = srctmp(j) + srccoefs(j,l)*pols(l)
          enddo
        enddo
        dst = sqrt(srctmp(3)**2 + srctmp(4)**2)
        srctmp(7) = srctmp(4)/dst
        srctmp(8) = -srctmp(3)/dst
        val = 0
        call fker(srctmp,ndt,targ,ndd,dpars,ndz,zpars,ndi,
     1      ipars,val)
        do j=1,k
          vals(j) = vals(j) + dst*w(i)*val*pols(j) 
        enddo
      enddo

      do i=1,k
        vals(i) = vals(i)*u
      enddo

      return
      end
c
c
c
c
c
ccc      subroutine zrmatmatt(m, n, a, k, b, c)
ccc      implicit double precision (a-h,o-z)
ccc      complex *16 :: a(n,m),c(k,m)
ccc      real *8 :: b(n,k)
ccc      complex *16, allocatable :: bz(:,:)
ccc      character *1 :: transa, transb
ccc      complex *16 :: alpha, beta
ccc
ccc
ccc
ccc      do i=1,m
ccc        do j=1,k
ccc          c(j,i) = 0
ccc          do l=1,n
ccc            c(j,i) = c(j,i) + a(l,i)*b(l,j)
ccc          enddo
ccc        enddo
ccc      enddo
ccc
ccc
ccc      return
ccc      end
ccc
ccc
ccc      subroutine zinsertmat(km, kn, amat, iloc, jloc, m, n, cmat)
ccc        implicit real *8 (a-h,o-z)
ccc        complex *16 :: amat(km,kn), cmat(m,n)
ccc
ccc        do j = 1,kn
ccc          do i = 1,km
ccc            cmat(iloc-1+i,jloc-1+j) = amat(i,j)
ccc          enddo
ccc        enddo
ccc
ccc      return
ccc      end

