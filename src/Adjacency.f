      subroutine adjtria3(ix,nodcon,nelcon,nnode,numel,nen,nedge)
c **********************************************************************
c *                                                                    *
c *   ADJTRIA3                                                         *
c *   --------                                                         *
c *                                                                    *
c *   Determina as adjacencias dos elementos (triangulos)              *
c *                                                                    *
c *   Parametros de entrada                                            *
c *   ---------------------                                            *
c *                                                                    *
c *   ix(nen+1,numel) - conetividades nodais dos elementos             *
c *   nodcon(nnode)   - nao definido (usado como arranjo auxiliar)     *
c *   nelcon(3,numel) - nao definido                                   *
c *   nnode           - numero de nos                                  *
c *   numel           - numero de elementos                            *
c *   nedge           - nao definido                                   *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   nelcon - elementos adjacentes ao elemento j sao dados por        *
c *            nelcon(i,j), i = 1,2,3.                                 *
c *            nelcon(i,j) = -1 (face pertence ao contorno)            *
c *   nedge  - numero de arestas                                       *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nnode,numel,nen,nedge,nel,is,no,imiss,no1,no2,nel2,is2
      integer ix(nen+1,*),nodcon(*),nelcon(3,*),isnod(2,3),no21,no22
      data isnod /1,2,2,3,3,1/
c.......................................................................
      nedge = 0
      do 100 nel = 1, numel
      do 100 is = 1, 3
         nelcon(is,nel) = 0
  100 continue
      do 200 no = 1, nnode
         nodcon(no) = 0
  200 continue
c ...........................
  300 continue
      imiss = 0
      do 400 nel = 1, numel
      do 400 is = 1, 3
         if (nelcon(is,nel) .ne. 0) go to 400
         no1 = ix(isnod(1,is),nel)
         no2 = ix(isnod(2,is),nel)
         if (nodcon(no1) .eq. 0 .and. nodcon(no2) .eq. 0) then
            nodcon(no1) = nel
            nodcon(no2) = nel
            imiss = 1
         endif
  400 continue
      do 550 nel = 1, numel
      do 550 is = 1, 3
         if (nelcon(is,nel) .ne. 0) go to 550
         no1  = ix(isnod(1,is),nel)
         no2  = ix(isnod(2,is),nel)
         nel2 = nodcon(no1)
         if(nel2 .gt. 0) then
         if(nel2.eq.nodcon(no2).and. nel2.ne.nel) then
            do 510 is2 = 1, 3
               no21 = ix(isnod(1,is2),nel2)
               no22 = ix(isnod(2,is2),nel2)
               if (no21 .eq. no2 .and. no22 .eq. no1) then
                  nelcon(is,nel)  = nel2
                  nelcon(is2,nel2)= nel
                  nodcon(no1) = 0
                  nodcon(no2) = 0
                  imiss = 1
                  nedge = nedge + 1
               endif
  510       continue
         endif
         endif
  550 continue
      do 600 nel = 1, numel
      do 600 is = 1, 3
         if (nelcon(is,nel) .ne. 0) go to 600
         no1 = ix(isnod(1,is),nel)
         no2 = ix(isnod(2,is),nel)
         if (nodcon(no1).eq.nodcon(no2) .and. nodcon(no1).eq.nel) then
            nelcon(is,nel)   = -1
            nodcon(no1) = 0
            nodcon(no2) = 0
            imiss = 1
            nedge = nedge + 1
         endif
  600 continue
      if (imiss .eq. 1) go to 300
c ......................................................................      
      return
      end
      subroutine adjquad4(ix,nodcon,nelcon,nnode,numel,nen,nedge)
c **********************************************************************
c *                                                                    *
c *   ADJQUAD4                                                         *
c *   --------                                                         *
c *                                                                    *
c *   Determina as adjacencias dos elementos (quadrilateros)           *
c *                                                                    *
c *   Parametros de entrada                                            *
c *   ---------------------                                            *
c *                                                                    *
c *   ix(nen+1,numel) - conetividades nodais dos elementos             *
c *   nodcon(nnode)   - nao definido (usado como arranjo auxiliar)     *
c *   nelcon(3,numel) - nao definido                                   *
c *   nnode           - numero de nos                                  *
c *   numel           - numero de elementos                            *
c *   nen             - numero de nos por elemento                     *
c *   nedge           - nao definido                                   *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   nelcon - elementos adjacentes ao elemento j sao dados por        *
c *            nelcon(i,j), i = 1,2,3,4.                               *
c *   nedge  - numero de arestas                                       *
c *                                                                    *
c **********************************************************************
      implicit none
      integer nnode,numel,nen,nedge,nel,is,no,imiss,no1,no2,nel2,is2
      integer ix(nen+1,*),nodcon(*),nelcon(4,*),isnod(2,4),no21,no22
      data isnod /1,2,2,3,3,4,4,1/
c.......................................................................
      nedge = 0
      do 100 nel = 1, numel
      do 100 is = 1, 4
         nelcon(is,nel) = 0
  100 continue
      do 200 no = 1, nnode
         nodcon(no) = 0
  200 continue
c ...........................
  300 continue
      imiss = 0
      do 400 nel = 1, numel
      do 400 is = 1, 4
         if (nelcon(is,nel) .ne. 0) go to 400
         no1 = ix(isnod(1,is),nel)
         no2 = ix(isnod(2,is),nel)
         if (nodcon(no1) .eq. 0 .or. nodcon(no2) .eq. 0) then
            nodcon(no1) = nel
            nodcon(no2) = nel
            imiss = 1
         endif
  400 continue
      do 550 nel = 1, numel
      do 550 is = 1, 4
         if (nelcon(is,nel) .ne. 0) go to 550
         no1  = ix(isnod(1,is),nel)
         no2  = ix(isnod(2,is),nel)
         nel2 = nodcon(no1)
         if(nel2 .gt. 0) then
         if(nel2.eq.nodcon(no2).and. nel2.ne.nel) then
            do 510 is2 = 1, 4
               no21 = ix(isnod(1,is2),nel2)
               no22 = ix(isnod(2,is2),nel2)
               if (no21 .eq. no2 .and. no22 .eq. no1) then
                  nelcon(is,nel)  = nel2
                  nelcon(is2,nel2)= nel
                  nodcon(no1) = 0
                  nodcon(no2) = 0
                  imiss = 1
                  nedge = nedge + 1
               endif
  510       continue
         endif
         endif
  550 continue
      do 600 nel = 1, numel
      do 600 is = 1, 4
         if (nelcon(is,nel) .ne. 0) go to 600
         no1 = ix(isnod(1,is),nel)
         no2 = ix(isnod(2,is),nel)
         if (nodcon(no1).eq.nodcon(no2) .and. nodcon(no1).eq.nel) then
            nelcon(is,nel)   = -1
            nodcon(no1) = 0
            nodcon(no2) = 0
            imiss = 1
            nedge = nedge + 1
         endif
  600 continue
      if (imiss .eq. 1) go to 300
      nedge = nedge + 2*numel
c ......................................................................      
      return
      end
      subroutine adjhexa8(ix,nodcon,nelcon,nnode,numel)
c **********************************************************************
c *                                                                    *
c *   HEXA8FACE - determina as adjacencias do hexaedro                 *
c *   ---------                                                        *
c *                                                                    *
c *   Parametros de entrada                                            *
c *   ---------------------                                            *
c *                                                                    *
c *   ix(9,numel)- conetividades nodais dos elementos                  *
c *   nodcon(nnode)   - nao definido (usado como arranjo auxiliar)     *
c *   nelcon(n,numel) - nao definido                                   *
c *   nnode           - numero de nos                                  *
c *   numel           - numero de elementos                            *
c *   nen             - numero de nos por elemento                     *
c *   nedge           - nao definido                                   *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   nelcon - elementos adjacentes                                    *
c *   need   - numero de arestas por elemento                          *
c *   nedge  - numero de arestas                                       *
c *                                                                    *
c **********************************************************************
      implicit none
      integer numel,nnode,ipass
      integer ix(9,numel),nodcon(nnode),nelcon(6,numel)
      integer node(4),i,j,k,el,miss,hexa8face
c ============================================================================
      do 55 i = 1, numel
      do 50 j = 1, 6
         nelcon(j,i) = -2
   50 continue
   55 continue
      do 60 i = 1, nnode
         nodcon(i) = -2
   60 continue
c .........................................................................
      ipass = 0
  100 continue
      miss = 0
      do 210 i = 1, numel
      do 200 j = 1, 6
         if(nelcon(j,i) .eq. -2) then
            call hexa8fnod(i,j,ix,node)
            if ((nodcon(node(1)).eq.-2).and.(nodcon(node(2)).eq.-2).and.
     .          (nodcon(node(3)).eq.-2).and.(nodcon(node(4)).eq.-2))then
                 nodcon(node(1)) = i
                 nodcon(node(2)) = i
                 nodcon(node(3)) = i
                 nodcon(node(4)) = i
                 miss = 1
                 goto 210
            endif
         endif
  200 continue
  210 continue
c ......................................................................
      do 320 i = 1, numel
      do 310 j = 1, 6
         if(nelcon(j,i).eq.-2) then
            call hexa8fnod(i,j,ix,node)
            el = nodcon(node(1))
            if(el .ne. -2) then
             if((nodcon(node(2)).eq.el).and.(nodcon(node(3)).eq.el).and.
     .          (nodcon(node(4)).eq.el).and.(el.ne.i)) then
                 k = hexa8face(el,ix,node)
                 if (k.eq.-1) then
                    print*,'*** <ADJHEXA8> Erro na vizinhanca ***'
                    stop
                 endif
                 nelcon(k,el) = i
                 nelcon(j,i)  = el
                 nodcon(node(1)) = -2
                 nodcon(node(2)) = -2
                 nodcon(node(3)) = -2
                 nodcon(node(4)) = -2
                 miss = 1             
             endif
            endif
         endif
  310 continue
  320 continue
c ......................................................................
      do 410 i = 1, numel
      do 400 j = 1, 6
         if(nelcon(j,i).eq.-2) then
            call hexa8fnod(i,j,ix,node)
            el = nodcon(node(1))
            if((nodcon(node(2)).eq.el).and.(nodcon(node(3)).eq.el).and.
     .         (nodcon(node(4)).eq.el).and.(el.eq.i)) then
                 nelcon(j,i) = -1
                 nodcon(node(1)) = -2
                 nodcon(node(2)) = -2
                 nodcon(node(3)) = -2
                 nodcon(node(4)) = -2
                 miss = 1        
            endif
         endif
  400 continue
  410 continue
      ipass = ipass + 1
      if (miss .eq. 1) goto 100
c ......................................................................
      return
      end
      subroutine Hexa8fnod(k,j,ix,node)
c **********************************************************************
c *                                                                    *
c *   HEXA8FNOD - determina os nos da face j do elemento k             *
c *   ---------                                                        *
c *                                                                    *
c *   Parametros de entrada                                            *
c *   ---------------------                                            *
c *                                                                    *
c *   k      - numero do elemento                                      *
c *   j      - numero da face do elemento                              *
c *   ix(9,numel)- conetividades nodais dos elementos                  *
c *   node(4)- nao definido                                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   node - nos da face j (numeracao local)                           *
c *                                                                    *
c **********************************************************************
      implicit none
      integer k,j,ix(9,*),node(*),fnode(4,6),i
      data fnode /1,2,3,4,5,8,7,6,1,5,6,2,4,3,7,8,1,4,8,5,2,6,7,3/
c ======================================================================
      do 100 i = 1, 4
         node(i) = ix(fnode(i,j),k)
  100 continue
      return
      end
      integer function hexa8face(k,ix,node)
c **********************************************************************
c *                                                                    *
c *   HEXA8FACE - determina a face do hexaedro k adjacente a face j    *
c *   ---------               cujos nos estao armazenados em node(4)   *
c *                                                                    *
c *   Parametros de entrada                                            *
c *   ---------------------                                            *
c *                                                                    *
c *   k          - numero do elemento adjacente                        *
c *   node - nos da face j (numeracao local)                           *
c *   ix(9,numel - conetividades nodais dos elementos                  *
c *                                                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer k,ix(9,*),node(*),ind(3,4),no(4),i,j
      data ind /2,3,4,3,4,1,4,1,2,1,2,3/
c ======================================================================
      hexa8face = -1
      do 200 i = 1, 6
         call hexa8fnod(k,i,ix,no)
         do 100 j = 1, 4
            if(no(1) .eq. node(j)) then
               if((no(2) .eq. node(ind(3,j))).and.
     .            (no(3) .eq. node(ind(2,j))).and.
     .            (no(4) .eq. node(ind(1,j)))) then
                   hexa8face = i
                   return
               endif
            endif
  100    continue
  200 continue
c ......................................................................  
      return
      end
      integer function edg(l,k,e,n)
c **********************************************************************
c *                                                                    *
c *   EDG - determina qual lado do triangulo l e adjacente             *
c *   ---   ao triangulo k.                                            *
c *                                                                    *
c *   Parametros de entrada                                            *
c *   ---------------------                                            *
c *                                                                    *
c *   l          - numero do elemento adjacente                        *
c *   k          - numero do elemento                                  *
c *   e(n,numel) - adjacencias dos elementos                           *
c *   n          - numero de arestas do elemento                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   edg - lado do elemento l adjacente ao elemento k                 *
c *                                                                    *
c **********************************************************************
      implicit none
      integer l,k,n,e(n,*),i
c ......................................................................
      do 10 i = 1, n
         if (e(i,l) .eq. k) then
            edg = i
            return
         endif
   10 continue
c ......................................................................   
      print*, '*** Erro na funcao EDG: elementos nao adjacentes ***'
      print*, l,e(1,l),e(2,l),e(3,l)
      print*, k,e(1,k),e(2,k),e(3,k)
      stop
      end
