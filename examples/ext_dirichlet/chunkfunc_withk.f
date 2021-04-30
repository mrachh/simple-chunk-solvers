c
c    modification of chunkfunc to make sure there are M points
c    per wavelength on a curve defined by funcurve
c
        subroutine chunkfunczk(zk, nptsperwave, eps, ifclosed, irefinel,
     1    irefiner, chsmall, ta, tb, funcurve, pars, nover, k, nch, 
     2    chunks, adjs, ders, ders2, hs)
        implicit real *8 (a-h,o-z)
        integer adjs(2,*),ifprocess(60000)
        double precision chunks(2,k,1),ders(2,k,1), ders2(2,k,1),
     1      xs(1000),ws(1000),xs2(1000),ws2(10000),u(10000),v(10000),
     2      u2(10000),v2(10000),ab(6,30000),hs(1),
     3      errs(1000),pars(1),
     4      errs0(1000)
        real *8 xcoefs(1000), ycoefs(1000), xpcoefs(1000)
        real *8 ypcoefs(1000), xpcoefs_out(1000)
        real *8 ypcoefs_out(1000)
        real *8 ch7(2,1000), der7(2,1000)
        complex *16 zk
c
        double precision, allocatable :: fvals(:,:), coefs(:,:)
c
        external funcurve
c
c       using a user-defined subroutine funcurve, split up the
c       curve into chunks. the calling sequence of funcurve
c       should be
c
c           funcurve(t,pars,x,y,dxdt,dydt,dxdt2,dydt2)
c
c       the routine will declare a chunk to be resolved when
c       chunks, ders, ders2 and dsdt = sqrt(ders(1,)**2 + ders(2,)**2)
c       have length 2k expansions with the last k coefficients
c       having relative-rmse less than eps.
c
c       NOTE: no effort whatsoever has been made to make this routine
c       efficient. it does what it does very deliberately, with some
c       redundancy.
c
c       NOTE 2: the routine returns chunks that are at most a factor
c       of two different in archlength than adjacent chunks
c
c       input:
c
c         eps - absolute precision to resolve the curve
c         ifclosed - switch telling the routine whether the curve is
c           closed or not, ifclosed=1 means a closed curve
c         irefinel - dyadically refine left (initial) end according 
c                    to tolerance chsmall
c                    1 -> refine, not otherwise
c         irefiner - dyadically refine right (terminal) end according 
c                    to tolerance chsmall
c                    1 -> refine, not otherwise
c         chsmall - if the curve is NOT closed, then this is the maximum
c           size of the intervals on either end. this parameter
c           is ignored if ifclosed=0.
c         ta,tb - assume that curve is parameterized by t \in [ta,tb)
c         funcurve - see above
c         pars - parameter array to send to funcurve, see above
c         nover - post-process oversampling factor, nover <= 1 will result
c             in no change, nover=2 will split each chunk in half (with
c             respect to arclength!!)
c         k - number of legendre nodes per chunk
c
c       output:
c
c         nch - total number of chunks created
c         chunks - points on each chunk, dimensioned chunks(2,k,nch)
c         adjs - adjacency information, adjs(1,i) is the chunk before
c             chunk i and adjs(2,i) is the chunk after chunk i
c         ders - first derivative w.r.t. t on each chunk,
c             dimensioned chunks(2,k,nch)
c         ders2 - second derivative w.r.t. t on each chunk,
c             dimensioned chunks(2,k,nch)
c         hs - weight parameter to account for the arbitrary underlying
c             parameterization, used in arclength calculation and
c             subsequent integration
c
c
        done=1
        pi=4*atan(done)
c
        maxchunks=55000
        do i=1,maxchunks
          ifprocess(i)=0
        enddo

c
c       construct legendre nodes and weights, k and 2k of them, as well
c       as the interpolation/coefficients matrices
c
        itype=2
        call legeexps(itype,k,xs,u,v,ws)
        call legeexps(itype,2*k,xs2,u2,v2,ws2)

c
c       . . . start chunking
c
        ab(1,1)=ta
        ab(2,1)=(ta+tb)/2
        ab(1,2)=(ta+tb)/2
        ab(2,2)=tb
        nch=2
        ifdone=1
        adjs(1,1)=-1
        adjs(2,1)=2
        adjs(1,2)=1
        adjs(2,2)=-1
        nchnew=nch
c
        allocate(fvals(2*k,7))
        allocate(coefs(2*k,7))
c
c       compute max length of a chunk that has nptsperwave points
c       per wavelength: 
c       k >  nwaves * nptsperwave = rmax*abs(zk)/(2*pi) * nptsperwave
c
        rmax = 2*pi*k/(abs(zk)*nptsperwave)
c
c
        maxiter=10000
        do 5000 ijk = 1,maxiter

c
c       loop through all existing chunks, if resolved store, if not split
c
          ifdone=1
          do 4600 ich=1,nchnew
c
            if (ifprocess(ich) .eq. 1) goto 4600
            ifprocess(ich)=1
c
            a=ab(1,ich)
            b=ab(2,ich)
            call chunklength(k,funcurve,pars,a,b,rlself)
c
            do i=1,2*k
              t=a+(b-a)*(xs2(i)+1)/2
              call funcurve(t,pars,fvals(i,1),fvals(i,2),
     1          fvals(i,3),fvals(i,4),fvals(i,5),fvals(i,6))
              fvals(i,7)=sqrt(fvals(i,3)**2+fvals(i,4)**2)
            enddo
c
            do i=1,7
              call chunkmatvec(2*k,2*k,u2,fvals(1,i),coefs(1,i))
              errs0(i) = 0
              errs(i) = 0
            enddo
c
            do i=1,7
              do j=1,k
                errs0(i) = errs0(i) + coefs(j,i)**2
                errs(i) = errs(i) + coefs(k+j,i)**2
              enddo
            enddo
c
            rmsemax = -1.0d0
            do i=1,7
              errs(i) = sqrt(errs(i)/errs0(i)/k)
              if (errs(i) .gt. rmsemax) rmsemax=errs(i)
            enddo

c
c       . . . mark as processed and resolved if less than eps
c
ccc            if (rmsemax .gt. eps) goto 2800
ccc            write(6,*) 'rlself = ',rlself
ccc            write(6,*) 'rmax = ',rmax
ccc            write(6,*) 'rmsemax = ',rmsemax
            if ( (rmsemax .gt. eps) .or. (rlself.gt.rmax)) goto 2800
            goto 4600
            
c
 2800 continue
ccc            write(6,*) 'not resolved = '

c
c       . . . if here, not resolved
c       divide - first update the adjacency list
c
            
            ifprocess(ich)=0
            ifdone=0
cccc            return
c
            if ((nch .eq. 1) .and. (ifclosed .gt. 0)) then
              adjs(1,nch)=2
              adjs(2,nch)=2
              adjs(1,nch+1)=1
              adjs(2,nch+1)=1
            endif
c
            if ((nch .eq. 1) .and. (ifclosed .le. 0)) then
              adjs(1,nch)=-1
              adjs(2,nch)=2
              adjs(1,nch+1)=1
              adjs(2,nch+1)=-1
            endif
c
            if (nch .gt. 1) then
              iold2=adjs(2,ich)
              adjs(2,ich)=nch+1
              if (iold2 .gt. 0) adjs(1,iold2)=nch+1
              adjs(1,nch+1)=ich
              adjs(2,nch+1)=iold2
            endif

c
c       now update the endpoints in ab
c
            ab(1,ich)=a
            ab(2,ich)=(a+b)/2
c
            nch=nch+1
            if (nch .gt. maxchunks) then
              call prinf('too many chunks in chunkfunc!*',done,0)
              call prinf('maxchunks=*',maxchunks,1)
              stop
            endif
c
            ab(1,nch)=(a+b)/2
            ab(2,nch)=b
c
 4500 continue
 4600 continue
c
          if ((ifdone .eq. 1) .and. (nchnew .eq. nch)) goto 5100
          nchnew=nch        
c
 5000 continue
 5100 continue

c
c       the curve should be resolved to precision eps now on
c       each interval ab(,i)
c       check the size of adjacent neighboring chunks - if off by a
c       factor of more than 2, split them as well. iterate until done.
c

        maxiter=1000
        do 9000 ijk=1,maxiter
c
          nchold=nch
          ifdone=1
          do 8600 i=1,nchold
c
            i1=adjs(1,i)
            i2=adjs(2,i)
ccc            write(6,*) 'nchold,i1,i2',nchold,i1,i2

c
c       calculate chunk lengths
c
            a=ab(1,i)
            b=ab(2,i)
            call chunklength(k,funcurve,pars,a,b,rlself)
c
            rl1=rlself
            rl2=rlself
c
            if (i1 .gt. 0) then
              a1=ab(1,i1)
              b1=ab(2,i1)
              call chunklength(k,funcurve,pars,a1,b1,rl1)
            endif
c
            if (i2 .gt. 0) then
              a2=ab(1,i2)
              b2=ab(2,i2)
              call chunklength(k,funcurve,pars,a2,b2,rl2)
            endif

c
c       only check if self is larger than either of adjacent blocks,
c       iterating a couple times will catch everything
c
            ifsplit=0
            sc = 2.05d0
            if (rlself .gt. sc*rl1) ifsplit=1
            if (rlself .gt. sc*rl2) ifsplit=1
ccc            write(6,*) 'ifsplit',ifsplit
            if (ifsplit .eq. 0) goto 8600
c
c       split chunk i now, and recalculate nodes, ders, etc
c
            ifdone=0
            a=ab(1,i)
            b=ab(2,i)
            ab2=(a+b)/2
c
            i1=adjs(1,i)
            i2=adjs(2,i)
c        
            adjs(1,i) = i1
            adjs(2,i) = nch+1
c
c       . . . first update nch+1
c
            adjs(1,nch+1) = i
            adjs(2,nch+1) = i2
c
c       . . . if there's an i2, update it
c
            if (i2 .gt. 0) then
              adjs(1,i2) = nch+1
            endif
c
            nch=nch+1
            if (nch .gt. maxchunks) then
                call prinf('too many chunks in chunkfunc!*',done,0)
                call prinf('maxchunks=*',maxchunks,1)
                stop
            endif
c
            ab(1,i)=a
            ab(2,i)=ab2
c
            ab(1,nch)=ab2
            ab(2,nch)=b
c
 8600 continue

c
          if (ifdone .ne. 0) goto 9100
c
 9000 continue
 9100 continue

c
c       go ahead and oversample by nover, updating
c       the adjacency information adjs along the way
c
        if (nover .le. 1) goto 6100

c
        do 6000 ijk=1,nover-1
c
        nchold=nch
        do 5600 i=1,nchold
c
        a=ab(1,i)
        b=ab(2,i)
c
c       find ab2 using newton such that 
c       len(a,ab2)=len(ab2,b)=half the chunk length
c
        call chunklength(k,funcurve,pars,a,b,rl)
c
        rlhalf=rl/2
        thresh=1.0d-8
        ifnewt=0
        ab0=(a+b)/2
c
        do 6600 iter=1,1000
c
        call chunklength(k,funcurve,pars,a,ab0,rl1)
        call funcurve(ab0,pars,x,y,dx,dy,ddx,ddy)
        dsdt=sqrt(dx**2+dy**2)
        ab1=ab0-(rl1-rlhalf)/dsdt
c
        err=rl1-rlhalf
        if (abs(err) .lt. thresh) ifnewt=ifnewt+1
c
        if (ifnewt .eq. 3) goto 6700
        ab0=ab1
c
 6600 continue
 6700 continue
c
        if (ifnewt .lt. 3) then
            call prin2('newton failed! interval not split.*',done,0)
            stop
        endif
c
        ab2=ab1
c
        i1=adjs(1,i)
        i2=adjs(2,i)
        adjs(2,i)=nch+1
        if (i2 .gt. 0) adjs(1,i2)=nch+1
c
        adjs(1,nch+1)=i
        adjs(2,nch+1)=i2
c
        ab(1,i)=a
        ab(2,i)=ab2
c
        nch=nch+1
        if (nch .gt. maxchunks) then
            call prinf('too many chunks in chunkfunc!*',done,0)
            call prinf('maxchunks=*',maxchunks,1)
            stop
        endif
c
        ab(1,nch)=ab2
        ab(2,nch)=b
c
 5600 continue
c
 6000 continue
 6100 continue
c
        if (ifclosed .eq. 1) goto 7700
c
c       if the curve is open, check the dyadic refinement at the
c       ends, first find the end segments
c
ccc        write(6,*) 'nch',nch
        do 6400 i =1,nch
        if (adjs(1,i) .lt. 0) ileft=i
        if (adjs(2,i) .lt. 0) iright=i
ccc        write(6,*) 'ileft,iright',ileft,iright
 6400 continue
c
        a1=ab(1,ileft)
        b1=ab(2,ileft)
        a2=ab(1,iright)
        b2=ab(2,iright)
c
        call chunklength(k,funcurve,pars,a1,b1,rlleft)
        call chunklength(k,funcurve,pars,a2,b2,rlright)
ccc        write(6,*) 'rlleft,rlright',rlleft,rlright
c
c       . . . dyadically split the left segment
c
        if (irefinel .eq.1) then
           do 7000 ijk=1,1000
c
           if (ijk .gt. 100) then
               call prinf('ifclosed = *', ifclosed, 1)
               call prinf('boom refining left end! ijk=*',ijk,1)
               call prinf('boom refining left end! rlleft=*',rlleft,1)
               stop
           endif
c
           if (rlleft .le. chsmall) goto 7100
ccc        write(6,*) 'refining left ',rlleft
           a=ab(1,ileft)
           b=ab(2,ileft)
           ab2=(a+b)/2
c
           i1=adjs(1,ileft)
           i2=adjs(2,ileft)
           adjs(2,ileft)=nch+1
           if (i2 .gt. 0) adjs(1,i2)=nch+1
c
           adjs(1,nch+1)=ileft
           adjs(2,nch+1)=i2
c
           ab(1,ileft)=a
           ab(2,ileft)=ab2
c
           nch=nch+1
           if (nch .gt. maxchunks) then
               call prinf('too many chunks in chunkfunc!*',done,0)
               call prinf('maxchunks=*',maxchunks,1)
               stop
           endif
c
           ab(1,nch)=ab2
           ab(2,nch)=b
c
           call chunklength(k,funcurve,pars,a,ab2,rlleft)
c        
 7000    continue
       endif
 7100 continue

c
c       . . . dyadically split the right segment
c
        if (irefiner .eq.1) then
           do 7400 ijk=1,1000
c
           if (ijk .gt. 100) then
               call prinf('boom refining right end! ijk=*',ijk,1)
               call prinf('boom right end! rlright=*',rlright,1)
               stop
           endif
c
ccc        write(6,*) 'refining right? ',rlright
           if (rlright .le. chsmall) goto 7500
ccc        write(6,*) 'yes, refining right ',rlright
           a=ab(1,iright)
           b=ab(2,iright)
           ab2=(a+b)/2
c
           i1=adjs(1,iright)
           i2=adjs(2,iright)
           adjs(2,iright)=nch+1
c
           adjs(1,nch+1)=iright
           adjs(2,nch+1)=i2
c
           ab(1,iright)=a
           ab(2,iright)=ab2
c
           nch=nch+1
           if (nch .gt. maxchunks) then
               call prinf('too many chunks in chunkfunc!*',done,0)
               call prinf('maxchunks=*',maxchunks,1)
               stop
           endif
c
           iright=nch
           ab(1,nch)=ab2
           ab(2,nch)=b
c
           call chunklength(k,funcurve,pars,ab2,b,rlright)
c
 7400    continue
      endif
 7500 continue
c
 7700 continue

c
c       up to here, everything has been done in parameter space, [ta,tb]
c       . . . finally evaluate the k nodes on each chunk, along with 
c       derivatives and chunk lengths
c
        do i = 1, nch
c
          a=ab(1,i)
          b=ab(2,i)
          hs(i)=(b-a)/2
c
          do j = 1, k
            t=a+(b-a)*(xs(j)+1)/2
            call funcurve(t, pars, x, y, dx, dy, dx2, dy2)
            chunks(1,j,i) = x
            chunks(2,j,i) = y
            ders(1,j,i) = dx
            ders(2,j,i) = dy
            ders2(1,j,i) = dx2
            ders2(2,j,i) = dy2
          enddo
        enddo
c
        return
        end
c
c
c
