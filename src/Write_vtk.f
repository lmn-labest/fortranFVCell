c *********************************************************************
c * WRITE_RES: escreve o resultados no formato do vtk                 *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * el     - conectividade com material                               *
c *  x     - coordenadas                                              *
c * nnode  - numero de nos                                            *
c * numel  - numero de elementos                                      *
c * nen    - numero de nos por elementos                              *
c * ndm    - numero de dimensoes                                      *
c * bvtk   - true BINARY vtk false ASCII vtk                          *
c *    u   - reslutado reais no elmentos                              *
c *    un  - reslutado media nodal                                    *
c * filein - nome do arquivo de saida                                 *
c * nCell  - nome do campo da celula                                  *
c * nNode  - nome do campo da no                                      *      
c * nout   - arquivo de saida                                         *
c *    cod - 1 - u                                                    *
c *    cod - 2 - kDu/dn                                               *
c *    cod - 3 - Du/dt                                                *
c *    cod - 4 - uel e unode                                          *
c *    cod - 5 - gradel e gradnode                                    *
c *    cod - 6 - fluxoel e fluxnode                                   *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine write_res_vtk(el,x,u,un,nnode,numel,ndm,nen,ndf,filein
     .                        ,nCell,nNod,bvtk,cod,t,istep,nout)
      use Malloc
      implicit none
      character nCell*30,nNod*30 
      real*8 x(ndm,*)
      integer el(nen+1,*)
      integer nnode,numel,nen,ndm,ndf
      real*8 u(ndf,*),un(ndf,*)
      logical bvtk
      integer istep
      real*8 t
c ... variaveis auxiliares      
      integer idum,cod,cod1,cod2,gdl,nel,nno
c ... ponteiros      
      integer*8 i_p
c .....................................................................      
      real fdum
      real*8 ddum
      character aux*30, aux1*30
c .....................................................................      
      character*80 fileout,filein  
      integer nout
c =====================================================================
c
c ===
      fileout = filein
      if(bvtk)then
        open(unit=nout,file=fileout,access='stream'
     .      ,form='unformatted',convert='big_endian')
      else
        open(unit=nout,file=fileout)
      endif  
c =====================================================================      
c
c === cabecalho
      write(aux,'(30a)')"Volume Finitos" 
      call head_vtk(aux,bvtk,t,istep,.true.,nout)
c =====================================================================
c
c === Coordenadas
      call coor_vtk(x,nnode,ndm,bvtk,nout)
c =====================================================================
c
c === Elementos
      call elm_vtk(el,numel,nen,bvtk,nout)
c =====================================================================
c
c === cell 
      call cell_data_vtk(numel,bvtk,nout)
c ... materiais      
      i_p = alloc_4('pi      ', 1,numel)
      do nel = 1, numel
        ia(i_p+nel-1) = el(nen+1,nel)
      enddo
      write(aux1,'(30a)')"mat" 
c ... cod = 1 variaveis interias
      gdl  = 1
      cod1 = 1
      cod2 = 1
      call cell_prop_vtk(ia(i_p),fdum,ddum,numel,aux1,ndm,gdl,cod1,cod2
     .                  ,bvtk,nout)
      i_p = dealloc('pi      ')
c ... id global dos elementos      
      i_p = alloc_4('pi      ', 1,numel)
      do nel = 1, numel
        ia(i_p+nel-1) = nel
      enddo
      write(aux1,'(30a)')"Gid_elmt" 
c ... cod = 1 variaveis interias
      gdl  = 1
      cod1 = 1
      cod2 = 1 
      call cell_prop_vtk(ia(i_p),fdum,ddum,numel,aux1,ndm,gdl,cod1,cod2
     .                  ,bvtk,nout)
      i_p = dealloc('pi      ')
c ... temperatura por celula     
      if( cod .eq. 4) then
c ...        
        write(aux1,'(30a)')nCell
        gdl  = ndf
        cod1 = 1
        cod2 = 3
        call cell_prop_vtk(idum,fdum,u,numel,aux1,ndm,gdl,cod1,cod2,bvtk
     .                    ,nout)
      endif
c ... gradiente por celula      
      if( cod .eq. 5) then
c ... erro da funcao exata     
c        write(aux1,'(30a)')"Error:GradEl" 
c        i_p = alloc_8('error   ', 1,numel)
c        gdl  = 1  
c        cod1 = 1
c        cod2 = 3
c        call errosol(u,ia(i_p),x,el,numel,nen,ndm,ndm,.false.)
c        call cell_prop_vtk(idum,fdum,ia(i_p),numel,aux1,ndm,gdl,cod1
c     .                    ,cod2,bvtk,nout)
c        i_p = dealloc('error   ')
c ...        
        write(aux1,'(30a)')nCell
        gdl  = ndm
        cod1 = 2
        cod2 = 3
        call cell_prop_vtk(idum,fdum,u,numel,aux1,ndm,gdl,cod1,cod2,bvtk
     .                    ,nout)
      endif  
c ... fluxo por celula      
      if( cod .eq. 6) then
c ...        
        write(aux1,'(30a)')nCell 
        gdl  = ndm
        cod1 = 2
        cod2 = 3
        call cell_prop_vtk(idum,fdum,u,numel,aux1,ndm,gdl,cod1,cod2,bvtk
     .                    ,nout)
      endif  
c ... velocidade por celula      
      if( cod .eq. 7) then
c ...        
        write(aux1,'(30a)')nCell 
        gdl  = ndm
        cod1 = 2
        cod2 = 3
        call cell_prop_vtk(idum,fdum,u,numel,aux1,ndm,gdl,cod1,cod2,bvtk
     .                    ,nout)
      endif
c ... massa especifica por celula     
      if( cod .eq. 8) then
c ...        
        write(aux1,'(30a)')nCell 
        gdl  = 3
        cod1 = 2
        cod2 = 3
        call cell_prop_vtk(idum,fdum,u,numel,aux1,ndm,gdl,cod1,cod2,bvtk
     .                    ,nout)
      endif
c ... Residou de massa por celula     
      if( cod .eq. 9) then
c ...        
        write(aux1,'(30a)')nCell 
        gdl  = 1
        cod1 = 1
        cod2 = 3
        call cell_prop_vtk(idum,fdum,u,numel,aux1,ndm,gdl,cod1,cod2,bvtk
     .                    ,nout)
      endif              
c ======================================================================      
c
c === nos
c ... id global dos nos      
      call point_data_vtk(nnode,bvtk,nout)
      i_p = alloc_4('pi      ', 1,nnode)
      do nno = 1, nnode
        ia(i_p+nno-1) = nno
      enddo
      write(aux1,'(30a)')"Gid_no" 
c ... cod = 1 variaveis interias
      gdl  = 1
      cod1 = 1
      cod2 = 1
      call pont_prop_vtk(ia(i_p),fdum,ddum,nnode,aux1,ndm,gdl,cod1,cod2
     .                    ,bvtk,nout)
      i_p = dealloc('pi      ')
c ...      
      if ( cod .eq. 4) then
        write(aux1,'(30a)') nNod 
        gdl  = ndf
        cod1 = 1
        cod2 = 3
        call pont_prop_vtk(idum,fdum,un,nnode,aux1,ndm,gdl,cod1,cod2
     .                    ,bvtk,nout)
c .....................................................................
c
c ...
      elseif ( cod .eq. 5) then
        write(aux1,'(30a)')nNod 
        gdl  = ndm
        cod1 = 2
        cod2 = 3
        call pont_prop_vtk(idum,fdum,un,nnode,aux1,ndm,gdl,cod1,cod2
     .                    ,bvtk,nout)
c .....................................................................
c
c ...
      elseif ( cod .eq. 6) then
        write(aux1,'(30a)')nNod 
        gdl  = ndm
        cod1 = 2
        cod2 = 3
        call pont_prop_vtk(idum,fdum,un,nnode,aux1,ndm,gdl,cod1,cod2
     .                    ,bvtk,nout)
c .....................................................................
c
c ...
      elseif ( cod .eq. 7) then
        write(aux1,'(30a)')nNod 
        gdl  = ndm
        cod1 = 2
        cod2 = 3
        call pont_prop_vtk(idum,fdum,un,nnode,aux1,ndm,gdl,cod1,cod2
     .                    ,bvtk,nout)
c .....................................................................
c
c ...
      elseif ( cod .eq. 8) then
        write(aux1,'(30a)')nNod  
        gdl  = 3   
        cod1 = 2
        cod2 = 3
        call pont_prop_vtk(idum,fdum,un,nnode,aux1,ndm,gdl,cod1,cod2
     .                    ,bvtk,nout)
c .....................................................................
c
c ...     
      elseif( cod .eq. 3) then
        call point_data_vtk(nnode,bvtk,nout)
        write(aux1,'(30a)')"Derivadas_Temporais"
        gdl  = ndf
        cod1 = 1
        cod2 = 3
        call pont_prop_vtk(idum,fdum,u,nnode,aux1,ndm,gdl,cod1,cod2
     .                    ,bvtk,nout)
      endif
c .....................................................................
c
c ======================================================================      
c
c === 
      close(nout) 
      return
      end
c ======================================================================      
c *********************************************************************
c
c *********************************************************************
c * WRITE_GEO: escreve o geomarita                                    *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * el     - conectividade com material                               *
c *  x     - coordenadas                                              *
c * pedgeF - condição de contorno(fluido)                             *
c * pedgeE - condição de contorno(energia)                            *
c * pedgeT1- condição de contorno(Transporte)                         *
c * nnode  - numero de nos                                            *
c * numel  - numero de elementos                                      *
c * nen    - numero de nos por elementos                              *
c * ndfT1  - grau de liberdade(transporte)                            *
c * ndfF   - grau de liberdade(fluido)                                *
c * ndfE   - grau de liberdade(energia)                               *
c * ndm    - numero de dimensoes                                      *
c * bvtk   - true BINARY vtk false ASCII vtk                          *
c * filein - nome do arquivo de saida                                 *   
c * nout   - arquivo de saida                                         *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine write_geo_vtk(el,x,pedgeF,pedgeE,pedgeT1
     .                        ,sedgeF,sedgeE,sedgeT1  
     .                        ,nnode,numel,ndm,nen,ndfT1,ndfF,ndfE
     .                        ,filein,bvtk,nout)
      use Malloc
      implicit none
      real*8 x(ndm,*)
      integer ndfT1,ndfF,ndfE
      integer el(nen+1,*)
      integer pedgeF(nen+1,*),pedgeE(nen+1,*),pedgeT1(nen+1,*)
      real*8  sedgeF((nen+1)*ndfF,*),sedgeE((nen+1)*ndfE,*)
      real*8  sedgeT1((nen+1)*ndfT1,*)
      integer nnode,numel,nen,ndm,nno
      logical bvtk
c ... variaveis auxiliares      
      integer idum,cod1,cod2,gdl,nel
c ... ponteiros      
      integer*8 i_p
c .....................................................................      
      real fdum
      real*8 ddum
      character aux*30, aux1*15
c .....................................................................      
      character*80 fileout,filein  
      integer nout
c =====================================================================
c
c === 
      fileout = filein
      if(bvtk)then
        open(unit=nout,file=fileout,access='stream'
     .      ,form='unformatted',convert='big_endian')
      else
        open(unit=nout,file=fileout)
      endif  
c =====================================================================      
c
c === cabecalho
      write(aux,'(30a)')"Mvf geometria" 
      call head_vtk(aux,bvtk,ddum,idum,.false.,nout) 
c =====================================================================
c
c === Coordenadas
      call coor_vtk(x,nnode,ndm,bvtk,nout)
c =====================================================================
c
c === Elementos
      call elm_vtk(el,numel,nen,bvtk,nout)
c =====================================================================
c
c === cell 
      call cell_data_vtk(numel,bvtk,nout)
c ... materiais 
      i_p = alloc_4('pi      ', 1,numel)
      do nel = 1, numel
        ia(i_p+nel-1) = el(nen+1,nel)
      enddo
      write(aux1,'(15a)')"mat" 
c ... cod = 1 variaveis interias
      gdl  = 1
      cod1 = 1
      cod2 = 1
      call cell_prop_vtk(ia(i_p),fdum,ddum,numel,aux1,ndm,gdl,cod1,cod2
     .                  ,bvtk,nout)
      i_p = dealloc('pi      ')
c ... id global dos elementos      
      i_p = alloc_4('pi      ', 1,numel)
      do nel = 1, numel
        ia(i_p+nel-1) = nel
      enddo
      write(aux1,'(15a)')"Gid_elmt" 
c ... cod = 1 variaveis interias
      gdl  = 1
      cod1 = 1
      cod2 = 1 
      call cell_prop_vtk(ia(i_p),fdum,ddum,numel,aux1,ndm,gdl,cod1,cod2
     .                  ,bvtk,nout)
      i_p = dealloc('pi      ')
c ======================================================================      
c
c === 
c ... contorno de fluido      
      if(ndfF .ne. 0 ) then
        write(aux1,'(15a)')"pedgeF"   
        cod1 = 1
        cod2 = 1
        gdl  = nen+1
        call cell_prop_vtk(pedgeF,fdum,ddum,numel,aux1,ndm,gdl,cod1
     .                    ,cod2,bvtk,nout)    
      endif
c ... contorno de Energia     
      if(ndfE .ne. 0 ) then
        write(aux1,'(15a)')"pedgeE"   
        cod1 = 1
        cod2 = 1
        gdl  = nen+1
        call cell_prop_vtk(pedgeE,fdum,ddum,numel,aux1,ndm,gdl,cod1
     .                    ,cod2,bvtk,nout)    
      endif
c ... contorno de Transporte  
      if(ndfT1 .ne. 0 ) then
        write(aux1,'(15a)')"pedgeT1"   
        cod1 = 1
        cod2 = 1
        gdl  = nen+1
        call cell_prop_vtk(pedgeT1,fdum,ddum,numel,aux1,ndm,gdl,cod1
     .                    ,cod2,bvtk,nout)    
      endif      
c ======================================================================
c
c ===      
c ... contorno de fluido      
      if(ndfF .ne. 0 ) then
        write(aux1,'(15a)')"sedgeF"   
        cod1 = 1
        cod2 = 3
        gdl  = (nen+1)*ndfF
        call cell_prop_vtk(idum,ddum,sedgeF,numel,aux1,ndm,gdl,cod1
     .                    ,cod2,bvtk,nout)    
      endif
c ... contorno de Energia     
      if(ndfE .ne. 0 ) then
        write(aux1,'(15a)')"sedgeE"   
        cod1 = 1
        cod2 = 3
        gdl  = nen+1
        call cell_prop_vtk(idum,fdum,sedgeE,numel,aux1,ndm,gdl,cod1
     .                    ,cod2,bvtk,nout)    
      endif
c ... contorno de Transporte  
      if(ndfT1 .ne. 0 ) then
        write(aux1,'(15a)')"sedgeT1"   
        cod1 = 1
        cod2 = 3
        gdl  = nen+1
        call cell_prop_vtk(idum,fdum,sedgeT1,numel,aux1,ndm,gdl,cod1
     .                    ,cod2,bvtk,nout)    
      endif
c =====================================================================
c
c === nos
c ... id global dos nos      
      call point_data_vtk(nnode,bvtk,nout)
      i_p = alloc_4('pi      ', 1,nnode)
      do nno = 1, nnode
        ia(i_p+nno-1) = nno
      enddo
c      write(aux1,'(15a)')'Gid_no1'
      aux1 = 'Gid_no'
c ... cod = 1 variaveis interias
      gdl  = 1
      cod1 = 1
      cod2 = 1
      call pont_prop_vtk(ia(i_p),fdum,ddum,nnode,aux1,ndm,gdl,cod1,cod2
     .                  ,bvtk,nout)
      i_p = dealloc('pi      ')
c ======================================================================   
c
c ===
c
c ====      
      close(nout) 
      return
      end
c =====================================================================      
c *********************************************************************
c
c *********************************************************************
c * WRITEMESH: escreve a malha particionada no formato do vtk         *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * np     - nos particionados                                        *
c * ep     - elementos particionado                                   *
c * el     - conectividade com material                               *
c *  x     - coordenadas                                              *
c * nnode  - numero de nos                                            *
c * numel  - numero de elementos                                      *
c * nen    - numero de nos por elementos                              *
c * ndm    - numero de dimensoes                                      *
c * bvtk   - true BINARY vtk false ASCII vtk                          *
c * nin    - arquivo de saida                                         *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c *********************************************************************
c     subroutine writemesh(np,ep,el,x,nnode,numel,nen,ndm,filein,bvtk
c    .                    ,nout)
c ===
c     use Malloc 
c     implicit none
c ... variaveis da malha      
c     integer nnode,numel,nen,ndm
c     integer np(nnode),ep(numel)
c     integer el(nen+1,numel)
c     real*8  x(ndm,nnode)
c     integer nel,nno
c ... locais     
c     integer i_p
c     character*15 aux1
c     character*30 aux
c ... variaveis dums
c     integer dum
c     real*8 ddum
c     real*4 fdum
c ... arquivo      
c     integer nout
c     character*80 fileout,name,filein
c     logical bvtk
c     integer cod,cod2,gdl
c =====================================================================
c
c ===
c     fileout = name(filein,0,106)
c     if(bvtk)then
c       open(unit=nout,file=fileout,access='stream'
c    .      ,form='unformatted',convert='big_endian')
c     else
c       open(unit=nout,file=fileout)
c     endif  
c     print*,fileout,nout
c =====================================================================
c
c === cabecalho
c     write(aux,'(30a)')"Malha part metis" 
c     call head_vtk(aux,bvtk,nout) 
c =====================================================================
c
c === Coordenadas
c     call coor_vtk(x,nnode,ndm,bvtk,nout)
c =====================================================================
c
c === Elementos
c     call elm_vtk(el,numel,nen,bvtk,nout)
c =====================================================================
c
c === cell 
c     call cell_data_vtk(numel,bvtk,nout)
c ... materiais      
c     i_p = alloc_4("p     ", 1,numel)
c     do nel = 1, numel
c       ia(i_p+nel-1) = el(nen+1,nel)
c     enddo
c     write(aux1,'(15a)')"mat" 
c ... cod = 1 variaveis interias
c     cod = 1
c     call cell_prop_vtk(ia(i_p),fdum,ddum,numel,aux1,cod,bvtk
c    .                  ,nout)
c     i_p = dealloc("p     ")
c .....................................................................
c
c ... particionamento
c     write(aux1,'(15a)')"part_metis" 
c     cod = 1
c     call cell_prop_vtk(ep,fdum,ddum,numel,aux1,cod,bvtk,nout)
c =====================================================================
c
c === nos  
c ...       
c     call point_data_vtk(nnode,bvtk,nout)
c     i_p = alloc_4("p     ", 1,nnode)
c     do nno = 1, nnode
c       ia(i_p+nno-1) = np(nno)
c     enddo
c     write(aux1,'(15a)')"part_por_no"
c ... gdb graus de liberdade
c     cod  1 escalar
c     cod2 1 int(4bytes) 
c     gdl =  1
c     cod =  1
c     cod2 = 1
c     call pont_prop_vtk(ia(i_p),fdum,ddum,nnode,aux1,ndm,gdl,cod,cod2
c    .                  ,bvtk,nout)
c     i_p = dealloc("p     ")
c .....................................................................
c =====================================================================
c     close(nout)
c     return
c     end
c =====================================================================
c *********************************************************************
c 
c *********************************************************************
c *                                                                   *
c * WRITEMESHPART:escreve a parte da malha do processp no formato vtk *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * np     - nos particionados                                        *
c * ep     - elementos particionado                                   *
c *  x     - coordenadas                                              *
c * nnode  - numero de nos                                            *
c * numel  - numero de elementos                                      *
c * nen    - numero de nos por elementos                              *
c * ndm    - numero de dimensoes                                      *
c * elLG   - mapa local -> global de elementos                        *
c * noLG   - mapa local -> global de nos                              *
c * bvtk   - true BINARY vtk false ASCII vtk                          *
c * nin    - arquivo de saida                                         *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c *********************************************************************
c     subroutine writemeshpart(np,ep,el,x,elLG,noLG,noGL,my_nnode
c    .                        ,my_numel,nen,ndm,filein,bvtk,nout
c    .                        ,my_rank)
c ===
c     use Malloc 
c     implicit none
c     include "elementos.fi"
c     integer my_nnode,my_numel,nen,ndm
c     integer np(*),ep(*)
c     integer elLG(*),noLG(*),noGL(*) 
c     integer el(nen+1,*)
c     integer i
c     integer my_rank
c     integer nout
c     integer dum,cod
c     real*4 fdum
c     real*8 ddum
c     integer nt1,nq1,nh1,nt2,nq2,nh2,ntr1,ntr2
c     integer i_lel,i_lno,i_ix,i_x,i_p
c     real*8 x(ndm,*)
c     character*30 aux
c     character*15 aux1
c     character*80 fileout,name,filein
c     logical bvtk
c =====================================================================
c
c ===
c     fileout = name(filein,my_rank,107)
c     if(bvtk)then
c       open(unit=nout,file=fileout,access='stream'
c    .      ,form='unformatted',convert='big_endian')
c     else
c       open(unit=nout,file=fileout)
c     endif  
c =====================================================================
c
c === Caracteriscas da funcao vtk
c     nt1  = ntria3(1)
c     nt2  = ntria3(2)
c     nq1  = nquad4(1)
c     nq2  = nquad4(2)
c     ntr1 = ntetra4(1)
c     ntr2 = ntetra4(2)
c     nh1  = nhexa8(1)
c     nh2  = nhexa8(2)
c     if(nt1 .gt. 0)then
c       ntria3(1) = my_numel 
c       ntria3(2) = 1   
c     else if(ntr1 .gt. 0)then  
c       ntetra4(1) = my_numel
c       ntetra4(2) = 1       
c     else if(nq1 .gt. 0)then  
c       nquad4(1) = my_numel
c       nquad4(2) = 1       
c     else if(nh1 .gt. 0)then  
c       nhexa8(1) = my_numel
c       nquad4(2) = 1
c     endif
c =====================================================================
c
c === gerando a malha local
c     i_lel = alloc_4('lel      ',nen+1,my_numel)
c     call localel(el,ia(i_lel),elLG,noGL,my_numel,nen)
c =====================================================================
c
c === gerando as coordenados locais
c     i_lno = alloc_8('lno      ',  ndm,my_nnode)
c     call localcoor(x,ia(i_lno),noLG,my_nnode,ndm)
c =====================================================================
c
c === gerando as propriedades dos materias
c     i_p = alloc_4('p        ',    1,my_numel)
c     call copy_mat(ia(i_p),ia(i_lel),my_numel,nen)
c =====================================================================
c
c === 
c
c ... head
c     write(aux,'(30a)')'Part mesh'
c     call head_vtk(aux,bvtk,nout)
c ... escrevendo coordenadas
c     call coor_vtk(ia(i_lno),my_nnode,ndm,bvtk,nout)
c ... escrevendo elementos      
c     call elm_vtk(ia(i_lel),my_numel,nen,bvtk,nout)
c ... escrevendo propriedades dos elementos    
c     call cell_data_vtk(my_numel,bvtk,nout)
c ... materias      
c     write(aux1,'(15a)')'mat'
c     cod = 1
c     call cell_prop_vtk(ia(i_p),fdum,ddum,my_numel,aux1,cod,bvtk,nout)
c ... partcionamento      
c     write(aux1,'(15a)')'part'
c     do i = 1, my_numel
c       ia(i_p+i-1) = ep(elLG(i)) 
c     enddo  
c     call cell_prop_vtk(ia(i_p),fdum,ddum,my_numel,aux1,cod,bvtk,nout)
c =====================================================================
c
c === 
c     ntria3(1)  = nt1
c     ntria3(2)  = nt2
c     nquad4(1)  = nq1
c     nquad4(2)  = nq2
c     ntetra4(1) = ntr1
c     ntetra4(2) = ntr2
c     nhexa8(1)  = nh1
c     nhexa8(2)  = nh2
c     i_p   = dealloc('p        ')
c     i_lno = dealloc('lno      ')
c     i_lel = dealloc('lel      ')
c     close(nout)
c     return
c     end
c =====================================================================
c
c *********************************************************************
c
