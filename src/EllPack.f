c *********************************************************************
c * ELPACKJA:                                                         *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ifEllPack - nao definido                                          *
c * ja        - nao definido                                          *
c * id        - conectividades dos vizinhos                           *
c * neq       - numero de equacoes                                    *
c * nad       - numero de termo nao nulos da matriz                   *
c * nshared   - numero de arestas                                     *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * ifEllpack - informacao do ellpack(nl(1)-nshared; nl(2) = neq + 1  *  
c * ja - arranjo apontador do elpack                                  *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine ellPackJa(ifEllPack,ja,id,neq,nad,nshared)
      implicit none
      integer neq,nshared,nad
      integer ja(nshared,*),id(nshared,*),ifEllPack(*)
c ... variaveis loacais
      integer i,j,nl
      ifEllPack(1) = nshared
      ifEllPack(2) = neq + 1
c ...
      nad  = 0
      do i = 1, neq
        do j = 1, nshared
          nl = id(j,i)
          if(nl .eq. -1 ) then
            ja(j,i) = 1       
          else
            nad     = nad + 1
            ja(j,i) = nl
          endif
        enddo
      enddo
c ...
      return
      end
c *********************************************************************
c
c *********************************************************************
c * ELPACK:                                                           *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ad     - nao definido                                             *
c * a      - nao definido                                             *
c * f      - nao definido                                             *
c * la     - valores locais da matriz A da celula P                   *
c * p      - valores locais do vetorde forca da celula P              *
c * ld     - numero da equacao da celula P                            *
c * nshared- numero de vizinhos P                                     *
c * unsysm  - matriz nao simetrica                                    *
c * forces - mantagem da vetor de forcas                              *
c * matrix - montagem da matrix                                       *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * ad   - diagonal principal da matriz A da celula P                 *
c * a    - valores fora da matriz principal da matriz A da celula P   *
c * f    - vetor de forca global da celula P                          *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine ellPack(ad,a,f,la,p,ld,nshared,forces,matrix)
      implicit none
      real*8 ad(*),a(nshared,*),f(*)
      real*8 la(*),p
      integer ld(*),nshared
      logical forces,matrix
c ... variaveis locais
      integer j,n,l,iCentral
      iCentral = nshared + 1
c ...
      n = ld(iCentral)
      if(forces)   f(n) =  p
      if(matrix) then 
        ad(n) = la(iCentral)
        do j = 1, nshared
          l = ld(j)
          if(l .gt. 0) then 
            a(j,n) = -la(j)
          else
            a(j,n) = 0.0d0  
          endif
        enddo
      endif
c .....................................................................
      return
      end
c *********************************************************************
