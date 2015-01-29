c *********************************************************************
c * ASSBLY : passagem dos valores locais da celula P para os arranjos *
c * globais (P sendo a um celula qualquer)                            *
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
c * cod    - 1 -    CSRD( ad - diag, a superior e inferior)           *   
c *          2 -    CSR ( a - diag, superior e ineferior)             *
c *          3 -    CSRC(a - diag, au - superior, al - ineferior)     *
c *          4 - ELLPACK(ad - diag prin, al - inferior e superior)    *
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
      subroutine assbly(ad,au,al,f,ia,ja,la,p,ld,nshared,cod,unsysm
     .                 ,forces,matrix)
      implicit none
      real*8 ad(*),au(*),al(*),f(*),la(*),p
      integer ia(*),ja(*),ld(nshared+1),nshared
      integer n,l,k,jk,iak,j,ipont,i,iaj,cod
      logical unsysm,forces,matrix
c ... sistema CSRD,CSR,CSRD
      if(cod .eq. 1 .or. cod .eq. 2 .or. cod .eq. 3) then
        call csr(ad,au,al,f,ia,ja,la,p,ld,nshared,cod,unsysm,forces
     .          ,matrix)
c .....................................................................
c
c ... ELLPACK
      elseif(cod .eq. 4) then
        call ellPack(ad,al,f,la,p,ld,nshared,forces,matrix)
      endif
c .....................................................................
c
c ...      
      return
      end
c *********************************************************************        
