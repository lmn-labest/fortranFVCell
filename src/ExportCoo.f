c *********************************************************************
c * EXPORTCOO : exporta o sistema de equacoes no formato coo          *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * iax    - diagonal principal de A                                  *
c * ja     - coeficiantes fora da diagonal principal de A             *
c * al     - coeficiantes da matrix da parte inferior                 *
c * ad     - coeficiantes da matrix da diogonal                       *
c * au     - coeficiantes da matrix da parte superior                 *
c * neq    - numero de equacoes                                       *
c * nad    - numero de elementos nao nulos fora da diagonal principal *
c * unsym  - matriz nao simetrica                                     *
c * prename- prefixo do arquivo                                       *
c * nout   - unidade do arquivo de saida                              *
c * it     - iteracao/timestep                                        *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine exportCoo(iax,ja,al,ad,au,b,neq,nad,unsym,prename
     .                    ,nout,it)
      use Malloc
      implicit none
      character*80 prename,name
      character*90 nameout
      integer*8 i_iaCoo,i_jaCoo,i_aCoo
      integer neq,nad,nout,it
      integer iax(*),ja(*)
      real*8 ad(*),al(*),au(*),b(*)
      logical unsym
c ...
      i_iaCoo  = alloc_4('linCooPc',1,neq+nad) 
      i_jaCoo  = alloc_4('colCooPc',1,neq+nad) 
      i_aCoo = alloc_4('aCooPc  ',1,neq+nad)
c ....................................................................
c 
c ...
      call csrToCoo(ia(i_iaCoo) ,ia(i_jaCoo) ,ia(i_aCoo)
     .             ,iax         ,ja
     .             ,al          ,ad      ,au
     .             ,neq         ,nad
     .             ,unsym       ,.false.)
c ....................................................................
c 
c ...
      nameOut = name(prename,it,52) 
      call writeCoo(ia(i_iaCoo),ia(i_jaCoo),ia(i_aCoo)
     .             ,b          ,neq        ,neq+nad
     .             ,nameOut    ,nout       ,.false.)
c ....................................................................
c 
c ...
      i_aCoo  = dealloc('aCooPc  ')
      i_jaCoo = dealloc('colCooPc') 
      i_iaCoo = dealloc('linCooPc') 
c ....................................................................
      return
      end
