c
      subroutine getdiscrete(nreg,nseg,zks,eps,k,nover,npw,rns,
     1            ilside,irside,iregionl,iregionr,irefinel,irefiner,
     1            alphas,betas,ifpr_geom,maxreg,maxc,
     1            chsmall,parsall,maxp,chunks,adjs,ders,ders2,hs,
     1            srcvals,srccoefs,whts,nch,maxseg,ipol)
      implicit real *8 (a-h,o-z)

      integer ilside(maxseg)
      integer irside(maxseg)
      integer irefinel(maxseg)
      integer irefiner(maxseg)
      integer iregionl(maxc)
      integer iregionr(maxc)
      integer adjs(2,maxc)
      real *8 chunks(2,k,maxc)
      real *8 ders(2,k,maxc)
      real *8 ders2(2,k,maxc)
      real *8 hs(maxc)
      real *8 srcvals(8,k,maxc)
      real *8 srccoefs(6,k,maxc)
      real *8 whts(k,maxc)

      complex *16 zks(0:maxreg)
      complex *16 rns(0:maxreg)
      complex *16 alphas(0:maxreg),betas(0:maxreg)
      real *8  parsall(2*maxp+3,*)
      external fcurve1
c
      ifclosed = 0
      nch = 0
      do iii = 1,nseg
         call chunkfunczk(zks(0),npw,eps,ifclosed,
     1       irefinel(iii),irefiner(iii),chsmall,
     1       parsall(2,iii),parsall(3,iii),fcurve1,parsall(1,iii),
     1       nover,k,nch1,chunks(1,1,nch+1),adjs(1,nch+1),
     2       ders(1,1,nch+1),ders2(1,1,nch+1),hs(nch+1))
         do i = 1,nch1
            iregionl(nch+i) = ilside(iii)
            iregionr(nch+i) = irside(iii)
         enddo
         nch = nch+nch1
          call prinf('nch1=*',nch1,1)
c
      enddo
      call prinf('nch=*',nch,1)
c
      if (ifpr_geom.eq.1) then
         open (unit=17,status='unknown',file='pts.m')
         write(17,*) ' xpts = ['
         do i = 1,nch
         do j = 1,k
            write(17,*) chunks(1,j,i), chunks(2,j,i)
         enddo
         enddo
         write(17,*) ' ];'
         write(17,*) ' plot(xpts(:,1),xpts(:,2),''.'')'
      endif

      call chunks_to_srcinfo(k,nch,chunks,ders,ders2,hs,srcvals,
     1  srccoefs,whts)

c
      alphas(0) = 1.0d0
      betas(0) = 1.0d0
      do ireg = 1,nreg
         alphas(ireg) = 1.0d0
         betas(ireg) = 1.0d0
         if (ipol.eq.1) betas(ireg) = 1.0d0/(rns(ireg)**2)
      enddo
c
      return
      end
c-----------------------------------------------------------
c

