      subroutine helm2d_comb_dir_targ(k,nch,n,srccoefs,srcinfo,
     1   zpars,sigma,nt,targs,pot) 
      implicit real *8 (a-h,o-z)
      real *8 srcinfo(8,k,nch),srccoefs(6,k,nch)
      complex *16 sigma(k,nch),pot(nt),zpars(3)
      real *8 targs(2,nt)
      
      real *8, allocatable :: cm(:,:),rads(:),radtmp(:)
      real *8, allocatable :: srcover(:,:,:)
      integer, allocatable :: row_ptr(:),col_ind(:)
      real *8 rfac
      real *8, allocatable :: tsover(:),wover(:),uover(:,:),vover(:,:)
      real *8, allocatable :: ximat(:,:)
      complex *16, allocatable :: zximat(:,:)
      complex *16, allocatable :: sigmaover(:,:)
      complex *16 zalpha,zbeta
      complex *16, allocatable :: charges(:),dipstr(:)
      real *8, allocatable :: dipvec(:,:),sources(:,:)

      real *8, allocatable :: tsquad(:),wquad(:)
      real *8, allocatable :: tstmp(:)
      real *8, allocatable :: work(:)
      real *8, allocatable :: ts0(:),w0(:),umat(:,:),vmat(:,:)
      complex *16, allocatable :: xtmp2(:)

      complex *16 tmp(10),pottmp

      external h2d_comb
      

      itype = 2
      allocate(ts0(k),w0(k),umat(k,k),vmat(k,k))
      call legeexps(itype,k,ts0,umat,vmat,w0)

      allocate(xtmp2(k))

      m = 20
      itype = 1
      allocate(tsquad(m),wquad(m))
      call legeexps(itype,m,tsquad,utmp,vtmp,wquad)

      allocate(tstmp(k+10))
      lwork = k*k*10+10*k + 100
      allocate(work(lwork))
c
c  find all close interactions
c
      allocate(cm(2,nch),rads(nch),radtmp(nch))
      
      call get_cm_rad(k,nch,srcinfo,cm,rads)
      rfac = 1.5d0
      do i=1,nch
        radtmp(i) = rads(i)*rfac
      enddo


      
      call findnearmem(cm,nch,radtmp,2,targs,nt,nnz)
      call prinf('nnz=*',nnz,1)

      allocate(row_ptr(nt+1))
      allocate(col_ind(nnz))

      call findnear(cm,nch,radtmp,2,targs,nt,row_ptr,col_ind)
      
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

      allocate(srcover(8,kover,nch),sigmaover(kover,nch))
c
c  oversample the geometry information
c
      alpha = 1.0d0
      beta = 0.0d0

      do ich=1,nch
        call dgemm('n','t',8,kover,k,alpha,srcinfo(1,1,ich),
     1   8,ximat,kover,beta,srcover(1,1,ich),8)
        call dgemm('n','t',2,kover,k,alpha,sigma(1,ich),
     1   2,ximat,kover,beta,sigmaover(1,ich),2)
      enddo

      call prinf('finished oversamling*',i,0)
c
c   call the fmm to compute all interactions
c
      ifcharge = 1
      ifdipole = 1
      if(abs(zpars(2)).le.1.0d-16) ifcharge = 0
      if(abs(zpars(3)).le.1.0d-16) ifdipole = 0
      
      nover = kover*nch
      allocate(charges(nover),dipstr(nover),dipvec(2,nover))
      allocate(sources(2,nover))

      do ich=1,nch
        do j=1,kover
          ipt = (ich-1)*kover + j
          sources(1,ipt) = srcover(1,j,ich)
          sources(2,ipt) = srcover(2,j,ich)
          dst = sqrt(srcover(3,j,ich)**2 + srcover(4,j,ich)**2)*wover(j)
          charges(ipt) = sigmaover(j,ich)*dst*zpars(2)
          dipstr(ipt) = sigmaover(j,ich)*dst*zpars(3)
          dipvec(1,ipt) = srcover(7,j,ich)
          dipvec(2,ipt) = srcover(8,j,ich)
        enddo
      enddo


      ifpgh = 0
      ifpghtarg = 1


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nt
        pot(i) = 0
      enddo
C$OMP END PARALLEL DO      
      
      eps = 1.0d-13
      nd = 1
      call hfmm2dpart(nd,eps,zpars(1),nover,sources,ifcharge,charges,
     1 ifdipole,dipstr,dipvec,ifpgh,tmp,tmp,tmp,nt,targs,ifpghtarg,
     2 pot,tmp,tmp)
c
c  now fix near corrections
c

      call get_fmm_thresh(2,nover,sources,2,nt,targs,thresh)
      thresh = thresh*abs(zpars(1))
      ndz = 3
      do i=1,nt
        nn = row_ptr(i+1)-row_ptr(i)
        do ii = row_ptr(i),row_ptr(i+1)-1
          ich = col_ind(ii)
          xtmp2 = 0
          call zadapquad(k,srccoefs(1,1,ich),targs(1,i),
     1       h2d_comb,2,ndd,dpars,ndz,zpars,ndi,ipars,umat,m,tsquad,
     2       wquad,xtmp2)
          pottmp = 0
          do j=1,k
            pot(i) = pot(i) + xtmp2(j)*sigma(j,ich)
          enddo

          istart = (ich-1)*kover+1
          pottmp = 0
          call h2d_directcdp_vec(nd,zpars(1),sources(1,istart),kover,
     1      charges(istart),dipstr(istart),dipvec(1,istart),
     2      targs(1,i),pottmp,thresh)
          pot(i) = pot(i) - pottmp
        enddo
      enddo
     


      
      return
      end
