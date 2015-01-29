c *********************************************************************
c * CSRIA:                                                            *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ia     - nao definido                                             *
c * id     - conectividades dos vizinhos                              *
c * neq    - numero de equacoes                                       *
c * nshared- numero de arestas                                        *
c * nad    - nao definido                                             *
c * upper  - parte superio no CSR                                     *
c * lower  - parte inferior no CSR                                    *
c * diag   - diagonal CSR                                             *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * ia - arranjo apontador do csr                                     *
c * nad- numero de termos nao-nulos fora da a diagonal principal      *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine csria(ia,id,neq,nshared,nad,upper,lower,diag)
      implicit none
      integer ia(*),id(nshared,*)
      integer neq,nshared,aux,nl,nad
      integer i,j
      logical upper,lower,diag
c ... gerando o arranjo ia 
      ia(1) = 1
      do i = 1, neq
        aux = 0
        do j = 1, nshared
          nl = id(j,i)
          if(nl .ne. -1) then
c ... parte superior
            if( lower .and. nl .lt. i) then
               aux =  aux + 1
           
c ... parte inferior
            elseif( upper .and. nl .gt. i) then
              aux = aux + 1
            endif
          endif
        enddo
        if( diag ) aux = aux + 1
        ia(i+1) = ia(i) + aux
      enddo
c .....................................................................
c
c ...
      nad = ia(neq+1) - ia(1)
c .....................................................................
c
c ...      
      return
      end
c ********************************************************************
c
c *********************************************************************
c * CSRJA:                                                            *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ia - arranjo apontador do csr                                     *
c * id - conectividades dos vizinhos                                  *
c * neq    - numero de equacoes                                       *
c * nshared- numero de aresta da malha                                *
c * upper  - parte superio no CSR                                     *
c * lower  - parte inferior no CSR                                    *
c * diag   - diagonal CSR                                             *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * ia - arranjo apontador do csr                                     *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine csrja(ia,ja,id,neq,nshared,upper,lower,diag)
      implicit none
      integer ia(*),ja(*),id(nshared,*),ipont
      integer neq,nshared,aux,nl
      integer i,j
      logical upper,lower,diag
c ... gerando o arranjo ja 
      do i = 1, neq
        aux   = 0
        ipont = ia(i)
c ... diagonal principal
        if(diag) then
          ja(ipont+aux) = i
          aux           = aux + 1
        endif
        do j = 1, nshared
          nl = id(j,i)
          if(nl .ne. -1) then
c ... parte superior
            if( nl .lt. i .and. lower) then
                ja(ipont+aux) = nl
                aux           = aux + 1
c ... parte inferior
            elseif(nl .gt. i .and. upper) then
                ja(ipont+aux) = nl
                aux           = aux + 1
            endif
          endif
        enddo
      enddo
c .....................................................................
c
c ...      
      return
      end
c ********************************************************************
c
c *********************************************************************
c * CSR:  Montagem do arranjo globa A,b para sistema no formato CSR   *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ad     - nao definido                                             *
c * au     - nao definido                                             *
c * al     - nao definido                                             *
c * f      - nao definido                                             *
c * ia     - ponteiro do arranjo do csr da matrix A                   *
c * ja     - arranjo csr da matriz A                                  *
c * la     - valores locais da matriz A da celula P                   *
c * p      - valores locais do vetorde forca da celula P              *
c * ld     - numero da equacao da celula P                            *
c * nshared- numero de vizinhos P                                     *
c * cod    - 1 - CSRD( ad - diag, a superior e inferior)              *
c *          2 - CSR ( a - diag, superior e ineferior)                *
c *          3 - CSRC(a - diag, au - superior, al - ineferior)        *
c * unsysm  - matriz nao simetrica                                    *
c * forces - mantagem da vetor de forcas                              *
c * matrix - montagem da matrix                                       *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * ad   - diagonal principal da matriz A da celula P                 *
c * au,al- valores fora da matriz principal da matriz A da celula P   *
c * f    - vetor de forca global da celula P                          *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine csr(ad,au,al,f,ia,ja,la,p,ld,nshared,cod,unsysm
     .              ,forces,matrix)
      implicit none
      real*8 ad(*),au(*),al(*),f(*),la(*),p
      integer ia(*),ja(*),ld(nshared+1),nshared
      integer n,l,k,jk,iak,j,ipont,i,iaj,cod
      logical unsysm,forces,matrix
c ...
c      do i = 1 , nshared
c        la(i) = -la(i)
c      enddo
c .....................................................................
c
c ... CSRD - CSR modificado
      if(cod .eq. 1) then
c 
c ... linha da matriz
        n      = ld(nshared+1)
        if(matrix) ad(n) = la(nshared+1)
        if(forces) f(n)  = p
        ipont  = ia(n)
        iak    = ia(n+1) - ia(n)
c ... loop nos vizinhos
        if(matrix) then
          do j = 1, nshared
            l = ld(j)
            if (l .gt. 0) then
              do i = 1, iak
                k  = ipont + i - 1
                jk = ja(k)
                if( l .eq. jk) au(k) = -la(j) 
              enddo 
            endif
          enddo
        endif
c .....................................................................
c
c ... CSR
      elseif(cod .eq. 2 ) then
c 
c ... linha da matriz
        n         = ld(nshared+1)
        if(forces) f(n)      = p
        ipont     = ia(n)
        iak       = ia(n+1) - ia(n)
c ... loop nos vizinhos
        do j = 1, nshared + 1
          l = ld(j)
          if (l .gt. 0) then
            do i = 1, iak
              k  = ipont + i - 1
              jk = ja(k)
              if( l .eq. jk .and. matrix) ad(k) = -la(j) 
            enddo 
          endif
        enddo
c .....................................................................
c
c ... CSRC
      elseif(cod .eq. 3 ) then
        n      = ld(nshared+1)
        ad(n)  = la(nshared+1)
        if(forces) f(n)  = p
        ipont  = ia(n)
        iak    = ia(n+1) - ia(n)
c ... loop nos vizinhos
        do j = 1, nshared
          l   = ld(j)
          if (l .gt. 0) then
            iaj = ia(l)
            if( l .gt. n) then
c              if(unsysm) au(iaj) = la(j)
            endif
            do i = 1, iak
              k  = ipont + i - 1
              jk = ja(k)
              if( l .eq. jk .and. matrix) al(k) = -la(j)
            enddo
          endif
        enddo
      endif
c .....................................................................
c
c ...      
      return
      end
c *********************************************************************
c
c *********************************************************************
c * BANDA_CSR: calculo da banda da matriz A                           *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ia     - arranjo apontador do csr                                 *
c * id     - conectividades dos vizinhos                              *
c * neq    - numero de equacoes                                       *
c * bandm  - nao definido                                             *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * bandm  - numero maximo da banda da matriz                         *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine banda_csr(ia,ja,neq,bandMax)
      implicit none
      integer ia(*),ja(*),neq,bandMax,bandMedia
      integer i,j,col,colMax
      bandMax   = 0
      bandMedia = 0
      do i = 1, neq
        do j = ia(i), ia(i+1)-1
          col       = abs(i-ja(j))
c          print*,i,col
          colMax    = max(col,colMax)
        enddo
        bandMedia = bandMedia + colMax
        bandMax   = max(bandMax,colMax)
        colMax    = 0
      enddo
      bandMedia = bandMedia/neq 
c     do i = 1, neq+1
c       print*,i,ia(i)
c     enddo
c     do i = 1, ia(neq+1)-ia(1)
c       print*,i,ja(i)
c     enddo
c     print*,band
      write(*,'(1x,a,i9)')'Danda maxima da matriz : ',bandMax 
      write(*,'(1x,a,i9)')'Danda media  da matriz : ',bandMedia 
      return
      end
c *********************************************************************
