      subroutine zgmres_solver(n,xmat,zid,rhs,numit,eps,niter,errs,rres,
     1   soln)
      implicit real *8 (a-h,o-z)
      integer n
      real *8 errs(numit+1)
      complex *16 xmat(n,n),zid,rhs(n),soln(n)
      complex *16, allocatable :: vmat(:,:),hmat(:,:),cs(:),sn(:)
      complex *16, allocatable :: svec(:),yvec(:),wtmp(:)
      complex *16 alpha,beta,ztmp,temp
      real *8 eps

      alpha = 1.0d0
      beta = 0.0d0
      allocate(vmat(n,numit+1),hmat(numit,numit),cs(numit),sn(numit))
      allocate(wtmp(n),svec(numit+1),yvec(numit+1))

c
c      compute norm of right hand side and initialize v
c 

      do i=1,numit
        cs(i) = 0
        sn(i) = 0
      enddo


      rb = 0
c
      do i=1,n
        rb = rb + abs(rhs(i))**2
      enddo
      rb = sqrt(rb)

      do i=1,n
        vmat(i,1) = rhs(i)/rb
      enddo

      svec(1) = rb

      do it=1,numit
        it1 = it + 1

c
c        NOTE:
c        replace this routine by appropriate layer potential
c        evaluation routine  
c

        do i=1,n
          wtmp(i) = 0
          do j=1,n
            wtmp(i) = wtmp(i) + xmat(i,j)*vmat(j,it)
          enddo
        enddo

        do k=1,it
          ztmp = 0
          do j=1,n
            ztmp = ztmp + wtmp(j)*conjg(vmat(j,k))
          enddo
          hmat(k,it) = ztmp

          do j=1,n
            wtmp(j) = wtmp(j)-hmat(k,it)*vmat(j,k)
          enddo
        enddo
          
        hmat(it,it) = hmat(it,it)+zid
        wnrm2 = 0
        do j=1,n
          wnrm2 = wnrm2 + abs(wtmp(j))**2
        enddo
        wnrm2 = sqrt(wnrm2)

        do j=1,n
          vmat(j,it1) = wtmp(j)/wnrm2
        enddo

        do k=1,it-1
          temp = cs(k)*hmat(k,it)+sn(k)*hmat(k+1,it)
          hmat(k+1,it) = -sn(k)*hmat(k,it)+cs(k)*hmat(k+1,it)
          hmat(k,it) = temp
        enddo

        ztmp = wnrm2

        call zrotmat_gmres(hmat(it,it),ztmp,cs(it),sn(it))
          
        hmat(it,it) = cs(it)*hmat(it,it)+sn(it)*wnrm2
        svec(it1) = -sn(it)*svec(it)
        svec(it) = cs(it)*svec(it)
        rmyerr = abs(svec(it1))/rb
        errs(it) = rmyerr

        if(rmyerr.le.eps.or.it.eq.numit) then

c
c            solve the linear system corresponding to
c            upper triangular part of hmat to obtain yvec
c
c            y = triu(H(1:it,1:it))\s(1:it);
c
          do j=1,it
            iind = it-j+1
            yvec(iind) = svec(iind)
            do l=iind+1,it
              yvec(iind) = yvec(iind) - hmat(iind,l)*yvec(l)
            enddo
            yvec(iind) = yvec(iind)/hmat(iind,iind)
          enddo
c
c          estimate x
c
          do j=1,n
            soln(j) = 0
            do i=1,it
              soln(j) = soln(j) + yvec(i)*vmat(j,i)
            enddo
          enddo


          rres = 0
          do i=1,n
            wtmp(i) = 0
            do j=1,n
              wtmp(i) = wtmp(i) + xmat(i,j)*soln(j)
            enddo
          enddo

          do i=1,n
            rres = rres + abs(zid*soln(i) + wtmp(i)-rhs(i))**2
          enddo
          rres = sqrt(rres)/rb
          niter = it
          return
        endif
      enddo


      return
      end
c
c
c
c
c

      subroutine zrotmat_gmres(a,b,c,s)
c-----------------------
c  Given a,b, compute sin(theta), cos(theta),
c  such that tan(theta) = b/a, note, the routine
c  implements a stabler version of computing
c  b/sqrt(a^2+b^2) and a/sqrt(b^2+a^2).
c
c  This routine is the complex version of rotmat_gmres
c
c  Input arguments:
c  
c    - a: double complex
c        cos scaling
c    - b: double complex 
c        sin scaling
c
c  Output arguments:
c    - c: double complex 
c        cos(theta)
c    - s: double complex 
c        sin(theta)
c        
c
c-----------------------
      implicit real *8 (a-h,o-z)
      complex *16 a,b,c,s,temp

      if(b.eq.0) then
        c = 1
        s = 0
      else if(abs(b).gt.abs(a)) then
        temp = a/b
        s = 1.0d0/sqrt(1.0d0+temp**2)
        c = temp*s
      else
        temp = b/a
        c = 1.0d0/sqrt(1.0d0+temp**2)
        s = temp*c
      endif

      return
      end

          


