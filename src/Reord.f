c*****************************Svn***************************************      
c*$Date: 2011-03-16 15:32:53 -0300 (qua, 16 mar 2011) $                 
c*$Rev: 914 $                                                           
c*$Author: henrique $                                                   
c***********************************************************************          
      subroutine reord(viz,iperm,numel,nshared,alf)
c **********************************************************************
c *                                                                    *
c *   REORD: reordenador de nos Reverse Cuthill-Mckee.                 *
c *   -----                                                            *
c **********************************************************************
      use Malloc
      implicit none
      integer viz(*),iperm(*),numel,nshared
      logical alf
c ... ponteiros      
      integer*8 i0,i1,i2,i3,i4
c ......................................................................
      integer i
c ......................................................................
c
c ...    Grafo da malha:
c
c        ia(i0)=>ip(nnode+1) - ip(i) indica a posicao em ips do primeiro
c                              no vizinho ao no i.
c        ia(i1) => ips(ipos) - contem as conectividades nodais de cada no 
c        conv_graph(viz,xdaj,adjncy,numel,nshared)
      if(alf) then
        i0 = alloc_4('iaux0   ',1,numel+1)
        i1 = alloc_4('iaux1   ',nshared,numel)      
        call conv_graph(viz,numel,nshared,ia(i0),ia(i1))
        i2 = alloc_4('iaux2   ',1,numel)                 
        i3 = alloc_4('iaux3   ',1,numel)
        i4 = alloc_4('iaux4   ',1,numel)
        call mzero(iperm,numel)
        call mzero(ia(i2),numel)
        call mzero(ia(i3),numel)
        call genrcm(numel,ia(i0),ia(i1),ia(i4),ia(i2),ia(i3))
        call permInverse(ia(i4),iperm,numel)
        i4 = dealloc('iaux4   ')
        i3 = dealloc('iaux3   ')
        i2 = dealloc('iaux2   ')
        i1 = dealloc('iaux1   ')
        i0 = dealloc('iaux0   ')
      else
        do i = 1, numel
          iperm(i) = i
        enddo
      endif 
c ......................................................................      
      return
      end
      subroutine genrcm (neq,xadj,adjncy,perm,mask,xls)
c **********************************************************************
c *                                                                    *
c *   general reverse cuthill mckee                                    *
c *                                                                    *
c *   purpose - genrcm finds the reverse cuthill-mckee ordering        *
c *      for a general graph. for each connected component in          *
c *      the graph, genrcm obtains the ordering by calling the         *
c *      subroutine rcm.                                               *
c *                                                                    *
c *   input parameters -                                               *
c *      neq - number of nodes                                         *
c *      (xadj, adjncy) - array pair containing the adjacency          *
c *             structure of the graph of the matrix.                  *
c *                                                                    *
c *   output parameters -                                              *
c *      perm - vector that contains the rcm ordering.                 *
c *                                                                    *
c *   working parameters -                                             *
c *      mask - is used to mark variables that have been               *
c *             numbered during the ordering process. it               *
c *             is initialized to 1, and set to zero as each           *
c *             node is numbered.                                      *
c *      xls - the index vector for a level structure. the level       *
c *            structure is stored in the currently unused spaces      *
c *            in the permutation vector perm.                         *
c *                                                                    *
c *   program subroutines - fnroot, rcm.                               *
c *                                                                    *
c **********************************************************************
      integer xadj(*),adjncy(*),mask(*),perm(*),xls(*)
      integer ccsize,i,neq,nlvl,num,root
c
      do i = 1, neq
         mask(i) = 1
      enddo
c   
      num = 1
      do i = 1, neq
c        -----------------------------------
c        for each masked connected component
c        -----------------------------------
         if (mask(i) .ne. 0) then    
            root = i
c           -------------------------------------------------------
c           first find a pseudo peripheral node root.
c           note that the level structure found by fnroot is
c           stored starting at perm(num). then rcm is called
c           to order the component using root as the starting node.
c           -------------------------------------------------------
           call fnroot(root,xadj,adjncy,mask,nlvl,xls,perm(num))
           call rcm(root,xadj,adjncy,mask,perm(num),ccsize,xls)
           num = num + ccsize
           if (neq .lt. num) return
         endif
      enddo
c 
      return
      end
      subroutine fnroot(root,xadj,adjncy,mask,nlvl,xls,ls)
c **********************************************************************
c *                                                                    *
c *   find pseudo-peripheral node                                      *
c *                                                                    *
c *   purpose - fnroot implements a modified version of the            *
c *      scheme by gibbs, pole, and stockmeyer to find pseudo-         *
c *      peripheral nodes. it determines such a node for the           *
c *      section subgraph specified by mask and root.                  *
c *                                                                    *
c *   input parameters -                                               *
c *      (xadj, adjncy) - adjacency structure pair for the graph.      *
c *      mask - specifies a section subgraph. nodes for which          *
c *             mask is zero are ignored by fnroot.                    *
c *                                                                    *
c *   updated parameters -                                             *
c *      root - on input, it (along with mask) defines the             *
c *             component for which a pseudo-peripheral node is        *
c *             to be found. on output, it is the node obtained.       *
c *                                                                    *
c *   output parameters -                                              *
c *      nlvl - is the number of levels in the level structure         *
c *             rooted at the node root.                               *
c *      (xls,ls) - the level structure array pair containing          *
c *                 the level structure found.                         *
c *                                                                    *
c *   program subroutines - rootls.                                    *
c *                                                                    *
c **********************************************************************
      integer adjncy(*),ls(*),mask(*),xls(*),xadj(*)
      integer ccsize, j, jstrt, k, kstop, kstrt, mindeg, nabor,
     .        ndeg, nlvl, node, nunlvl, root
c
c     ---------------------------------------------
c     determine the level structure rooted at root.
c     ---------------------------------------------
      call rootls ( root, xadj, adjncy, mask, nlvl, xls, ls )
      ccsize = xls(nlvl+1) - 1
      if (nlvl .eq. 1 .or. nlvl .eq. ccsize) return
c     ----------------------------------------------------
c     pick a node with minimum degree from the last level.
c     ----------------------------------------------------
  100 continue
c
      jstrt = xls(nlvl)
      mindeg = ccsize
      root = ls(jstrt)
      if (ccsize .eq. jstrt) goto 400
c
      do j = jstrt, ccsize
        node = ls(j)
        ndeg = 0
        kstrt = xadj(node)
        kstop = xadj(node+1) - 1
c        
        do k = kstrt, kstop
          nabor = adjncy(k)
          if (mask(nabor) .gt. 0) ndeg = ndeg + 1
        enddo
        if ( ndeg .ge. mindeg ) then
          root = node
          mindeg = ndeg
        endif 
      enddo
c     ---------------------------------------
c     and generate its rooted level structure.
c     ---------------------------------------
  400 continue
c
      call rootls ( root, xadj, adjncy, mask, nunlvl, xls, ls )
      if ( nunlvl .le. nlvl ) return
      nlvl = nunlvl
      if ( nlvl .lt. ccsize ) goto 100
      return
      end
      subroutine rcm(root,xadj,adjncy,mask,perm,ccsize,deg)
c **********************************************************************
c *                                                                    *
c *   reverse cuthill-mckee ordering                                   *
c *                                                                    *
c *   purpose - rcm numbers a connected component specified            *
c *      by mask and root, using the rcm algorithm.                    *
c *      the numbering is to be started at the node root.              *
c *                                                                    *
c *   input parameters -                                               *
c *      (xadj, adjncy) - adjacency structure pair for the graph.      *
c *      root - is the node that defines the connected component       *
c *             and it is used as the starting node for the rcm        *
c *             ordering.                                              *
c *                                                                    *
c *   updated parameters -                                             *
c *      mask - only those nodes with nonzero input mask values        *
c *             are considered by the routine. the nodes numbered      *
c *             by rcm will have their mask values set to zero.        *
c *                                                                    *
c *   output parameters -                                              *
c *      perm - will contain the rcm ordering. level structure         *
c *      ccsize - is the size of the connected component that          *
c *               has been numbered by rcm.                            *
c *                                                                    *
c *   working parameters -                                             *
c *      deg - is a temporary vector used to hold the degree of        *
c *            the nodes in the section graph specified by mask        *
c *            and root.                                               *
c *                                                                    *
c *   program subroutines - degree.                                    *
c *                                                                    *
c **********************************************************************
      integer adjncy(*),deg(*),mask(*),perm(*),xadj(*)
      integer ccsize, fnbr, i, j, jstop, jstrt, k, l, lbegin,
     .        lnbr, lperm, lvlend, nbr, node,  root
c
c     ----------------------------------------------
c     find the degrees of the nodes in the component
c     specified by mask and root.
c     ----------------------------------------------
      call degree (root,xadj,adjncy,mask,deg,ccsize,perm)
      mask(root) = 0
      if ( ccsize .le. 1 ) return
      lvlend = 0
      lnbr = 1
c     ----------------------------------------------
c     lbegin and lvlend point to the beginning and
c     the end of the current level respectively.
c     ----------------------------------------------
  100 continue
      lbegin = lvlend + 1
      lvlend = lnbr
      do i = lbegin, lvlend
c        ------------------------------
c        for each node in current level
c        ------------------------------
         node = perm(i)
         jstrt = xadj(node)
         jstop = xadj(node+1) - 1
c        -----------------------------------------
c        find the unnumbered neighbors of node.
c        fnbr and lnbr point to the first and last
c        unnumbered neighbors respectively of the
c        current node in perm.
c        -----------------------------------------
         fnbr = lnbr + 1
         do j = jstrt, jstop
            nbr = adjncy(j)
            if ( mask(nbr) .ne. 0 ) then    
               lnbr = lnbr + 1
               mask(nbr) = 0
               perm(lnbr) = nbr
            endif   
         enddo   
         if ( fnbr .lt. lnbr ) then    
c           -------------------------------------------
c           sort the neighbors of node in increasing
c           order by degree. linear insertion is used.
c           -------------------------------------------
            k = fnbr
c            
  300       continue
            l = k
            k = k + 1
            nbr = perm(k)
  400       if ( l .lt. fnbr ) goto 500
            lperm = perm(l)
            if ( deg(lperm) .le. deg(nbr) ) goto 500
            perm(l+1) = lperm
            l = l - 1
            goto 400
c  
  500       continue   
c  
            perm(l+1) = nbr
            if ( k .lt. lnbr ) goto 300
         endif
      enddo
c 
      if ( lnbr .gt. lvlend ) goto 100
c     ---------------------------------------
c     we now have the cuthill mckee ordering.
c     reverse it below ...
c     ---------------------------------------
      k = ccsize/2
      l = ccsize
      do i = 1, k
         lperm = perm(l)
         perm(l) = perm(i)
         perm(i) = lperm
         l = l - 1
      enddo
      return
      end
      subroutine rootls(root,xadj,adjncy,mask,nlvl,xls,ls)
c **********************************************************************
c *                                                                    *
c *   rooted level structure                                           *
c *                                                                    *
c *   purpose - rootls generates the level structure rooted            *
c *      at the input node called root. only those nodes for           *
c *      which mask is nonzero wil be considered.                      *
c *                                                                    *
c *   input parameters -                                               *
c *      root - the node at which the level structure is to            *
c *             be rooted.                                             *
c *      (xadj, adjncy) - adjacency structure pair for the             *
c *             given graph.                                           *
c *      mask - specifies a section subgraph. nodes for which          *
c *             mask is zero are ignored.                              *
c *                                                                    *
c *   output parameters -                                              *
c *      nlvl - is the number of levels in the level structure.        *
c *      (xls,ls) - array pair for the rooted level structure.         *
c *                                                                    *
c **********************************************************************
      integer adjncy(1), ls(1), mask(1), xls(1), xadj(1)
      integer i, j, jstop, jstrt, lbegin, ccsize, lvlend, lvsize,
     .        nbr, nlvl, node, root
c
c     -------------------
c     initialization ...
c     -------------------
      mask(root) = 0
      ls(1) = root
      nlvl = 0
      lvlend = 0
      ccsize = 1
c     ------------------------------------------------------
c     lbegin is the pointer to the beginning of the current
c     level, and lvlend points to the end of this level.
c     ------------------------------------------------------
  200 lbegin = lvlend + 1
      lvlend = ccsize
      nlvl = nlvl + 1
      xls(nlvl) = lbegin
c     ------------------------------------------------------
c     generate the next level by finding all the masked
c     neighbors of nodes in the current level.
c     ------------------------------------------------------
      do 400 i = lbegin, lvlend
         node = ls(i)
         jstrt = xadj(node)
         jstop = xadj(node+1) - 1
         if ( jstop .lt. jstrt ) goto 350
            do 300 j = jstrt, jstop
               nbr = adjncy(j)
               if ( mask(nbr) .eq. 0 ) goto 250
                  ccsize = ccsize + 1
                  ls(ccsize) = nbr
                  mask(nbr) = 0
  250          continue
  300       continue
  350    continue
  400 continue
c     ------------------------------------------------------
c     compute the current level width.
c     if it is nonzero, generate the next level.
c     ------------------------------------------------------
      lvsize = ccsize - lvlend
      if ( lvsize .gt. 0 ) goto 200
c     --------------------------------------------------------
c     reset mask to one for the nodes in the level structure.
c     --------------------------------------------------------
      xls(nlvl+1) = lvlend + 1
      do 500 i = 1, ccsize
         node = ls(i)
         mask(node) = 1
  500 continue
      return
      end
      subroutine degree ( root, xadj, adjncy, mask, deg, ccsize, ls )
c **********************************************************************
c *                                                                    *
c *   degree in masked component                                       *
c *                                                                    *
c *   purpose - degree computes the degrees of the nodes in            *
c *      the connected component specified by mask and root.           *
c *      nodes for which mask is zero are ignored.                     *
c *                                                                    *
c *   input parameters -                                               *
c *      root - is the input node that defines the component.          *
c *      (xadj, adjncy) - adjacency structure pair.                    *
c *      mask - specifies a section subgraph.                          *
c *                                                                    *
c *   output parameters -                                              *
c *      deg - array containing the degrees of the nodes in            *
c *            the component.                                          *
c *      ccsize - size of the component specified by mask and root.    *
c *                                                                    *
c *   working parameters -                                             *
c *      ls - a temporary vector used to store the nodes of the        *
c *           component level by level.                                *
c *                                                                    *
c **********************************************************************
      integer adjncy(*), deg(*), ls(*), mask(*), xadj(*)
      integer i, j, jstop, jstrt, lbegin, ccsize, lvlend, lvsize,
     .        nbr, ideg, node, root
c
c     --------------------------------------------------
c     initialization ...
c     the array xadj is used as a temporary marker to
c     indicate which nodes have been considered so far.
c     --------------------------------------------------
      ls(1) = root
      xadj(root) = -xadj(root)
      lvlend = 0
      ccsize = 1
c     ------------------------------------------------------
c     lbegin is the pointer to the beginning of the current
c     level, and lvlend points to the end of this level.
c     ------------------------------------------------------
  100 continue
      lbegin = lvlend + 1
      lvlend = ccsize
c     ------------------------------------------------------
c     find the degrees of nodes in the current level,
c     and at the same time, generate the next level.
c     ------------------------------------------------------
      do i = lbegin, lvlend
         node = ls(i)
         jstrt = -xadj(node)
         jstop =  iabs(xadj(node+1)) - 1
         ideg = 0
         if ( jstop .lt. jstrt ) goto 300
            do 200 j = jstrt, jstop
               nbr = adjncy(j)
               if ( mask(nbr) .eq. 0 ) goto 200
               ideg = ideg + 1
               if ( xadj(nbr) .lt. 0 ) goto 200
               xadj(nbr) = -xadj(nbr)
               ccsize = ccsize + 1
               ls(ccsize) = nbr
  200       continue
  300    deg(node) = ideg
      enddo   
c     ------------------------------------------------------
c     compute the current level width.
c     if it is nonzero, generate the next level.
c     ------------------------------------------------------
      lvsize = ccsize - lvlend
      if ( lvsize .gt. 0 ) goto 100
c     --------------------------------------------------------
c     reset xadj to its correct sign and return.
c     --------------------------------------------------------
      do i = 1, ccsize
         node = ls(i)
         xadj(node) = -xadj(node)
      enddo
      return
      end
c *********************************************************************
      subroutine permInverse(perm,iperm,numel)
      implicit none
      integer perm(*),iperm(*)
      integer numel,i
      do i = 1, numel
        iperm(perm(i)) = i
      enddo
      return
      end