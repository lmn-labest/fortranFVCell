c **********************************************************************
c * VIZ: determina a vizinha dos elementos                             *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * ix(nen+1) - conectividade dos elementos                            *
c * i_nelcon  - nao definido                                           *
c * nnode     - numero de nos da malha                                 *
c * numel     - numero de elementos                                    *
c * nen       - numero de aresta por elementos                         *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * i_nelcon  - ponteiro a o arranjo que guarda as vizinhas dos elem   *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine viz(ix,i_nelcon,nnode,numel,nen)
      Use Malloc
      implicit none
      integer ix(nen+1,*)
      integer*8 i_nelcon,i_nodcon
      integer nnode,numel,nen,dum1
c ... 
      if( nen .eq. 3) then
        i_nelcon    = alloc_4('nelcon  ',nen,numel)
        i_nodcon    = alloc_4('nodcon  ',  1,nnode)
c ... identificacao dos vizinhos
        call adjtria3(ix,ia(i_nodcon),ia(i_nelcon),nnode,numel
     .               ,nen,dum1)
c .....................................................................
c
c ...
        i_nodcon    = dealloc('nodcon  ')
        i_nelcon    =  locate('nelcon  ')
c .....................................................................
c
c ... 
      elseif( nen .eq. 4) then
        i_nelcon    = alloc_4('nelcon  ',nen,numel)
        i_nodcon    = alloc_4('nodcon  ',  1,nnode)
c ... identificacao dos vizinhos
        call adjquad4(ix,ia(i_nodcon),ia(i_nelcon),nnode,numel
     .               ,nen,dum1)
c .....................................................................
c
c ...
        i_nodcon    = dealloc('nodcon  ')
        i_nelcon    =  locate('nelcon  ')
      endif      
c .....................................................................
      return
      end
