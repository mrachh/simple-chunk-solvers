       subroutine get_cm_rad(k,nch,srcinfo,cm,rads)
       implicit real *8 (a-h,o-z)
       integer k,nch
       real *8 srcinfo(8,k,nch)
       real *8 cm(2,*),rads(*)

       do i=1,nch
         cm(1,i) = 0
         cm(2,i) = 0
         do j=1,k
           cm(1,i) = cm(1,i) + srcinfo(1,j,i)
           cm(2,i) = cm(2,i) + srcinfo(2,j,i)
         enddo
         cm(1,i) = cm(1,i)/k
         cm(2,i) = cm(2,i)/k

         rmax = 0
         do j=1,k
           dx = cm(1,i) - srcinfo(1,j,i)
           dy = cm(2,i) - srcinfo(2,j,i)
           rr = dx**2 + dy**2
           if(rr.gt.rmax) rmax = rr
         enddo
         rads(i) = sqrt(rmax)
       enddo


       return
       end
c
c
c
c


      subroutine findnearslow(xyzs,ns,rads,ndt,targets,nt,row_ptr,
     1   col_ind)
c
cc      identify targets which are contained in 
c       |xyzs(:,i)-targets(:,j)|<=rads(i)
c       This is a dense O(ns \times nt) algorithm
c
c       Calling sequence variables
c       xyzs    in: real *8 (2,ns)
c               location of the sources
c
c       ns      in: integer
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c
c       targets in: real *8(2,nt)
c               location of the targets
c
c       nt      in: integer
c               number of targets
c
c
c       OUTPUT
c       row_ptr    out: integer(nt+1)
c                itptr(i) is the starting point in the iflg
c                array for the list of sources
c                relevant for target i
c
c       col_ind  out: integer(nnz) (nnz is computed using mem routine)
c                iflg(row_ptr(i):row_ptr(i+1)) is the list
c                of sources which satisfy
c                d(s_{j},t_{i}) <= rs_{j}
c               
c-------------------------------

      implicit real *8 (a-h,o-z)
      integer ndt
      real *8 xyzs(2,ns),targets(ndt,nt),rads(ns)
      integer row_ptr(*),col_ind(*)
      integer, allocatable :: nlst(:)

      allocate(nlst(nt))

      do i=1,nt
        nlst(i) = 0
        do j=1,ns
          rr = (xyzs(1,j)-targets(1,i))**2 + 
     1          (xyzs(2,j)-targets(2,i))**2 
          if(rr.le.rads(j)**2) then
            nlst(i) = nlst(i) + 1
          endif
        enddo
      enddo

      row_ptr(1) = 1
      do i=1,nt
        row_ptr(i+1) = row_ptr(i)+nlst(i)
      enddo

      do i=1,nt
        nlst(i) = 0
      enddo

      do i=1,nt
        do j=1,ns
          rr = (xyzs(1,j)-targets(1,i))**2 + 
     1         (xyzs(2,j)-targets(2,i))**2 

          if(rr.le.rads(j)**2) then
            col_ind(row_ptr(i)+nlst(i)) = j
            nlst(i) = nlst(i) + 1
          endif
        enddo
      enddo



      return
      end
c------------------------------------------------------------------      

      subroutine findnearslowmem(xyzs,ns,rads,ndt,targets,nt,nnz)
c
cc      identify all sources which are contained in 
c       |xyzs(:,i)-targets(:,j)|<=rads(i)
c       This is a dense O(ns \times nt) algorithm.
c
c       Calling sequence variables
c       xyzs    in: real *8 (3,ns)
c               location of the sources
c
c       ns      in: integer
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c
c       targets in: real *8(3,nt)
c               location of the targets
c
c       nt      in: integer
c               number of targets
c
c
c       OUTPUT
c       nnz     out: integer
c               number of elements in the flag array
c-------------------------------

      implicit real *8 (a-h,o-z)
      integer ndt
      real *8 xyzs(2,*),targets(ndt,*),rads(*)
      integer nnz

      nnz = 0

      do i=1,nt
        do j=1,ns
          rr = (xyzs(1,j)-targets(1,i))**2 + 
     1         (xyzs(2,j)-targets(2,i))**2 

          if(rr.le.rads(j)**2) then
            nnz = nnz+1
          endif
        enddo
      enddo

      return
      end
c------------------------------------------------------------------      
      subroutine findnearmem(xyzs,ns,rads,ndt,targets,nt,nnz)
c
cc      identify all sources which are contained in 
c       |xyzs(:,i)-targets(:,j)|<=rads(i)+radt(j).
c       We use area queries to identify a nearly minimal set of
c       sources to loop over relevant for each target.
c       This is just a memory management code for flagaqt3
c       to identify length of array needed for flagging sources
c       
c
c       Calling sequence variables
c
c       xyzs    in: real *8 (3,ns)
c               location of the sources
c
c       ns      in: integer
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c
c       targets in: real *8(3,nt)
c               location of the targets
c
c       nt      in: integer
c               number of targets
c
c       radt    in: real *8(ns)
c               radii associated with the targets
c
c       targets in: real *8(3,nt)
c               location of the targets
c
c       nt      in: integer
c               number of targets
c
c       OUTPUT
c       nnz     out: integer
c               number of elements in the flag array
c
c       Note that in this code targets are treated as expansion
c       centers
c               
c-------------------------------

       implicit real *8 (a-h,o-z)
       integer ndt
       real *8 xyzs(2,ns),rads(ns),targets(ndt,nt)
       real *8 xtmp(2)
       real *8, allocatable :: rstmp(:)
       real *8, allocatable :: targs(:,:)

       integer nnz

cc       tree variables
c 
       real *8, allocatable :: centers(:,:),boxsize(:)
       integer, allocatable :: ilevel(:)
       integer, allocatable :: itree(:)
       integer ipointer(32),ltree

       logical res



       allocate(rstmp(ns))

       allocate(targs(2,nt))

       do i=1,nt
         targs(1,i) = targets(1,i)
         targs(2,i) = targets(2,i)
       enddo



       do i=1,ns
         rstmp(i) = 2*rads(i)
       enddo


       mnbors = 9
       mnlist2 = 3*mnbors

       rttmp = 0


       idivflag = 0
       ndiv = 2
       isep = 1
       nlmax = 200
       nbmax = 0

       nlevels = 0
       nboxes = 0
       mhung = 0
       mnbors = 0
       mnlist1 = 0
       mnlist2 = 0
       mnlist3 = 0
       mnlist4 = 0
       ntmp = 0

       call maketree2dmem(ier,xyzs,ns,rstmp,targs,nt,xtmp,ntmp,rttmp,
     1    idivflag,ndiv,nlmax,nbmax,nlevels,nboxes,mhung,ltree)


       allocate(centers(2,nboxes),itree(ltree),boxsize(0:nlevels))
  
       call maketree2d(xyzs,ns,rstmp,targs,nt,xtmp,ntmp,rttmp,
     1    idivflag,ndiv,mhung,nlevels,nboxes,centers,boxsize,
     2    itree,ltree,ipointer,mnlist1,mnlist2,mnlist3,mnlist4)

       allocate(ilevel(nboxes))

       do ilev=0,nlevels
         do ibox=itree(2*ilev+1),itree(2*ilev+2)
           ilevel(ibox) = ilev
         enddo
       enddo

       nnz = 0

       do ilev=0,nlevels
         do ibox=itree(2*ilev+1),itree(2*ilev+2)
           nchild = itree(ipointer(3)+ibox-1)
           if(nchild.eq.0) then 
             itstart = itree(ipointer(12)+ibox-1)
             itend = itree(ipointer(13)+ibox-1)
             do it = itstart,itend
               itarg = itree(ipointer(6)+it-1)
               nhunglistsrc = itree(ipointer(28)+ibox-1)
               do ii=1,nhunglistsrc
                 is = itree(ipointer(29)+mhung*(ibox-1)+ii-1)
                 rr = (targets(1,itarg)-xyzs(1,is))**2+ 
     1                 (targets(2,itarg)-xyzs(2,is))**2
                 if(rr.le.rads(is)**2) nnz = nnz + 1              
               enddo

               nlist1 = itree(ipointer(18)+ibox-1)
               do j=1,nlist1
                 jbox = itree(ipointer(19)+mnlist1*(ibox-1)+j-1)
                 isstart = itree(ipointer(10)+jbox-1)
                 isend = itree(ipointer(11)+jbox-1)
                 do ii = isstart,isend
                   is = itree(ipointer(5)+ii-1)
                   rr = (targets(1,itarg)-xyzs(1,is))**2+ 
     1                    (targets(2,itarg)-xyzs(2,is))**2
                   if(rr.le.rads(is)**2) nnz = nnz + 1
                 enddo
               enddo
             enddo
           endif
         enddo
       enddo


       return
       end
c-----------------------------------------------------
      subroutine findnear(xyzs,ns,rads,ndt,targets,nt,row_ptr,
     1       col_ind)  
cc      identify all sources which are contained in 
c       |xyzs(:,i)-targets(:,j)|<=rads(i).
c       We use area queries to identify a nearly minimal set of
c       sources to loop over relevant for each target.
c       
c
c       Calling sequence variables
c
c       xyzs    in: real *8 (3,ns)
c               location of the sources
c
c       ns      in: integer
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c
c       targets in: real *8(3,nt)
c               location of the targets
c
c       nt      in: integer
c               number of targets
c
c       row_ptr    out: integer(nt+1)
c                row_ptr(i) is the starting point in the col_ind
c                array for the list of sources
c                relevant for target i
c
c       col_ind     out: integer(nnz) (nnz is computed using
c                                      findnearmem)
c                col_ind(row_ptr(i):row_ptr(i+1)-1) is the list
c                of sources which satisfy
c                d(s_{j},t_{i}) <= rs_{j}
c               
c               
c-------------------------------
       implicit real *8 (a-h,o-z)
       integer ndt
       real *8 xyzs(2,*),rads(*),targets(ndt,*)
       real *8, allocatable :: rstmp(:)

       integer row_ptr(*),col_ind(*)
       integer, allocatable :: nlst(:)
c
cc       tree variables
c 
       real *8, allocatable :: centers(:,:),boxsize(:)
       real *8, allocatable :: targs(:,:)
       integer, allocatable :: itree(:)
       integer ipointer(32), ltree
       integer, allocatable :: ilevel(:)

       allocate(rstmp(ns))

       do i=1,ns
         rstmp(i) = 2*rads(i)
       enddo

       rttmp = 0.0d0

       allocate(targs(2,nt))

       do i=1,nt
         targs(1,i) = targets(1,i)
         targs(2,i) = targets(2,i)
       enddo


       mnbors = 27
       mnlist2 = 7*mnbors

       idivflag = 0
       ndiv = 2
       isep = 1
       nlmax = 200
       nbmax = 0

       nlevels = 0
       nboxes = 0
       mhung = 0
       ntmp = 0
       mnlist1 = 0
       mnlist2 = 0
       mnlist3 = 0
       mnlist4 = 0
       mnbors = 0

       call maketree2dmem(ier,xyzs,ns,rstmp,targs,nt,xtmp,ntmp,rttmp,
     1    idivflag,ndiv,nlmax,nbmax,nlevels,nboxes,mhung,ltree)


       allocate(centers(2,nboxes),itree(ltree),boxsize(0:nlevels))
  
       call maketree2d(xyzs,ns,rstmp,targs,nt,xtmp,ntmp,rttmp,
     1    idivflag,ndiv,mhung,nlevels,nboxes,centers,boxsize,
     2    itree,ltree,ipointer,mnlist1,mnlist2,mnlist3,mnlist4)


       allocate(ilevel(nboxes))

       do ilev=0,nlevels
          do ibox=itree(2*ilev+1),itree(2*ilev+2)
              ilevel(ibox) = ilev
          enddo
       enddo

       allocate(nlst(nt))


       do i=1,nt
         nlst(i) = 0
       enddo

       do ilev=0,nlevels
         do ibox=itree(2*ilev+1),itree(2*ilev+2)
           nchild = itree(ipointer(3)+ibox-1)
           if(nchild.eq.0) then
             itstart = itree(ipointer(12)+ibox-1)
             itend = itree(ipointer(13)+ibox-1)

             do it = itstart,itend
               itarg = itree(ipointer(6)+it-1)
               nhunglistsrc = itree(ipointer(28)+ibox-1)
               do ii=1,nhunglistsrc
                 is = itree(ipointer(29)+mhung*(ibox-1)+ii-1)
                 rr = (targets(1,itarg)-xyzs(1,is))**2+ 
     1                 (targets(2,itarg)-xyzs(2,is))**2
                 if(rr.le.rads(is)**2) 
     1               nlst(itarg) = nlst(itarg) + 1              
               enddo

               nlist1 = itree(ipointer(18)+ibox-1)
               do j=1,nlist1
                 jbox = itree(ipointer(19)+mnlist1*(ibox-1)+j-1)
                 isstart = itree(ipointer(10)+jbox-1)
                 isend = itree(ipointer(11)+jbox-1)
                 do ii = isstart,isend
                   is = itree(ipointer(5)+ii-1)
                   rr = (targets(1,itarg)-xyzs(1,is))**2+ 
     1                     (targets(2,itarg)-xyzs(2,is))**2
                   if(rr.le.rads(is)**2) 
     1                    nlst(itarg) = nlst(itarg) + 1
                 enddo
               enddo
             enddo
           endif
         enddo
       enddo

       row_ptr(1) = 1
       do i=1,nt
          row_ptr(i+1) = row_ptr(i) + nlst(i)
          nlst(i) = 0
       enddo



       do ilev=0,nlevels
         do ibox=itree(2*ilev+1),itree(2*ilev+2)
           nchild = itree(ipointer(3)+ibox-1)
           if(nchild.eq.0) then
             itstart = itree(ipointer(12)+ibox-1)
             itend = itree(ipointer(13)+ibox-1)

             do it = itstart,itend
               itarg = itree(ipointer(6)+it-1)
               nhunglistsrc = itree(ipointer(28)+ibox-1)
               do ii=1,nhunglistsrc
                 is = itree(ipointer(29)+mhung*(ibox-1)+ii-1)
                 rr = (targets(1,itarg)-xyzs(1,is))**2+ 
     1                 (targets(2,itarg)-xyzs(2,is))**2
                 if(rr.le.rads(is)**2) then
                   col_ind(row_ptr(itarg)+nlst(itarg)) =is
                   nlst(itarg) = nlst(itarg)+1
                 endif
               enddo

               nlist1 = itree(ipointer(18)+ibox-1)
               do j=1,nlist1
                 jbox = itree(ipointer(19)+mnlist1*(ibox-1)+j-1)
                 isstart = itree(ipointer(10)+jbox-1)
                 isend = itree(ipointer(11)+jbox-1)
                 do ii = isstart,isend
                   is = itree(ipointer(5)+ii-1)
                   rr = (targets(1,itarg)-xyzs(1,is))**2+ 
     1                  (targets(2,itarg)-xyzs(2,is))**2
                   if(rr.le.rads(is)**2) then 
                     col_ind(row_ptr(itarg)+nlst(itarg)) =is
                     nlst(itarg) = nlst(itarg)+1
                   endif
                 enddo
               enddo
             enddo
           endif
         enddo
       enddo


       return
       end
c-----------------------------------------------------
c
c
c
c
c

      subroutine rsc_to_csc(nsrc,ntarg,nnz,row_ptr,col_ind,
     1    col_ptr,row_ind,iper)
      implicit real *8 (a-h,o-z)      
c
c       convert a row-sparse compressed representation to
c       a column sparse compressed representation
c

       integer npatches,nnz,row_ptr(ntarg+1),col_ind(nnz)
       integer col_ptr(nsrc+1),row_ind(nnz),iper(nnz)
       integer, allocatable :: nslr(:)

      allocate(nslr(nsrc))
      lflg = nnz
      do i=1,nsrc
         nslr(i) = 0
      enddo

      do i=1,lflg
        nslr(col_ind(i)) = nslr(col_ind(i)) + 1
      enddo


      col_ptr(1) = 1
      do i=2,nsrc+1
         col_ptr(i) = col_ptr(i-1)+nslr(i-1)
      enddo


      do itarg=1,ntarg
         do ictr=row_ptr(itarg),row_ptr(itarg+1)-1
           jsrc = col_ind(ictr)

           iper(col_ptr(jsrc)) = ictr

           row_ind(col_ptr(jsrc)) = itarg
           col_ptr(jsrc) = col_ptr(jsrc) + 1 
         enddo
      enddo

      col_ptr(1) = 1
      do i=2,nsrc+1
         col_ptr(i) = col_ptr(i-1)+nslr(i-1)
      enddo

      return
      end
c-----------------------------------------------------------
