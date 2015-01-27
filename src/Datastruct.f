c *********************************************************************
c * DATASTRUC: Estrutura de dados para matriz esparca                 *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ia     - nao definido                                             *
c * i_ia   - nao define                                               *
c * i_ja   - nao definido                                             *
c * neq    - numero de equacao                                        *
c * nad    - nao definido                                             *
c * nshared- numero de aresta da malha                                *
c * code   - 1 - csr                                                  *
c * banda  - nao definido                                             *
c * unsym  - matriz nao simetrica                                     *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * i_ia  - ponteiro para o arranjo apontador do csr                  *
c * i_ja  - ponteiro para o arranjo do csr                            *
c * nad   - numero de termos nao-nulos fora da adiagonaal principal   *
c * banda - numero da banda da matriz                                 *
c * ----------------------------------------------------------------- *
c *********************************************************************      
      subroutine datastruc(i_id,i_ia,i_ja,i_ad,i_au,i_al,neq,nad
     .                    ,nshared,code,banda,unsysm
     .                    ,strId,strIa,strJa,strAd,strAu,strAl)
      use Malloc
      implicit none
      integer*8 i_ia,i_ja,i_ad,i_au,i_al,i_id
      integer neq,nad,nshared,code,banda
      logical unsysm,upper,lower,diag
      character*8 strId,strIa,strJa,strAd,strAu,strAl
      i_ad = 1
      i_au = 1
      i_al = 1
c ... CSRD( ad, a)
      if(code .eq. 1) then
c ...
        lower = .true.
        if(unsysm) then
          upper = .true.
        else 
          upper = .false. 
        endif 
        diag  = .false.
c .....................................................................
c
c ...
        i_ia   = alloc_4(strIa,  1,neq+1)  
        call csria(ia(i_ia),ia(i_id),neq,nshared,nad,upper,lower,diag)
        i_ja   = alloc_4(strJa,  1,nad)
        call csrja(ia(i_ia),ia(i_ja),ia(i_id),neq,nshared,upper,lower
     .            ,diag)
c .....................................................................
c
c ... reordenando grafo em ordem crescente
        call sortgraph(ia(i_ia),ia(i_ja),neq)
c .....................................................................
c
c ...   calculo da banda maxima da matriz A
        call banda_csr(ia(i_ia),ia(i_ja),neq,banda)
c .....................................................................
c
c ... alocacao da matriz
        i_ad   = alloc_8(strAd,  1,neq)
        i_au   = alloc_8(strAu,  1,nad)
        call azero(ia(i_ad),neq)
        call azero(ia(i_au),nad)
c         i_id = dealloc(strId)
        i_ia = locate(strIa)
        i_ja = locate(strJa)
        i_ad = locate(strAd)
        i_au = locate(strAu)
        i_al = i_au
c .....................................................................
c
c ... CSR ( ad - toda a matriz)
      elseif(code .eq. 2) then
c ...
        upper = .true.
        lower = .true.
        diag  = .true.
c .....................................................................
c
c ...
        i_ia   = alloc_4(strIa,  1,neq+1)  
        call csria(ia(i_ia),ia(i_id),neq,nshared,nad,upper,lower,diag)
        i_ja   = alloc_4(strJa,  1,nad)
        call csrja(ia(i_ia),ia(i_ja),ia(i_id),neq,nshared,upper,lower
     .            ,diag)
c .....................................................................
c
c ... reordenando grafo em ordem crescente 
c     (nao ordenar o grafo para manter o termo do diagonal com primeiro
c      termo da linha i)
c        call sortgraph(ia(i_ia),ia(i_ja),neq)
c .....................................................................
c
c ... alocacao da matriz
        i_ad = alloc_8(strAd,  1,nad)
        call azero(ia(i_ad),nad)
c          i_id = dealloc(strId)
        i_ia = locate(strIa)
        i_ja = locate(strJa)
        i_ad = locate(strAd)
        i_al = 1
        i_au = 1 
c ....................................................................
c
c ... CSRC (ad-diagonal principal;au-parte superior;al-parte inferior)
      elseif(code .eq. 3) then
c ...
        upper = .false.
        lower = .true. 
        diag  = .false. 
c .....................................................................
c
c ...
        i_ia   = alloc_4(strIA,  1,neq+1)  
        call csria(ia(i_ia),ia(i_id),neq,nshared,nad,upper,lower,diag)
        i_ja   = alloc_4(strJa,  1,nad)
        call csrja(ia(i_ia),ia(i_ja),ia(i_id),neq,nshared,upper,lower
     .            ,diag)
c .....................................................................
c
c ... reordenando grafo em ordem crescente
        call sortgraph(ia(i_ia),ia(i_ja),neq)
c .....................................................................
c
c ... alocacao da matriz
        i_ad   = alloc_8(strAd,  1,neq)
        i_al   = alloc_8(strAl,  1,nad)
        if(unsysm) then 
          i_au = alloc_8(strAu,  1,nad)
          call azero(ia(i_au),nad)
        endif
        call azero(ia(i_ad),neq)
        call azero(ia(i_au),nad)
c       i_id = dealloc(strId)
        i_ia = locate(strIa)
        i_ja = locate(strJa)
        i_ad = locate(strAd)
        i_al = locate(strAl)
        if(unsysm) i_au = locate(strAu)
        i_au   = i_al
      endif       
c ....................................................................      
      return
      end
c ....................................................................
c ********************************************************************
