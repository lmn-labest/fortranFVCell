c *********************************************************************
c * PMSHAPE: imprime o prefil da matrix A no formato vtk(diagonal     *
c * principal e parte inferior)                                       *
c * ------------------------------------------------------------------*
c * Parametros de entrada:                                            *
c * ------------------------------------------------------------------*
c * ip   - apontado da linha do csr                                   *
c * la   - colunas do csr                                             *
c * neq  - numero de equacoes                                         *
c * nout - Arquivo de saida                                           *
c * ------------------------------------------------------------------*
c * Parmetros de saida:                                               *
c * ------------------------------------------------------------------*
c * ------------------------------------------------------------------*
c *********************************************************************
      subroutine pmshape(ip,ja,neq,nout)
      use Malloc
      implicit none
      integer ip(*),ja(*)
      integer neq,nad
      integer nout
      integer*8 i_p
c ... cabecalho do arquivo vtk      
      write(nout,'(a)') '# vtk DataFile Version 3.0'
      write(nout,'(a)') 'Perfil da matriz'
      write(nout,'(a)') 'ASCII'
      write(nout,'(a)') 'DATASET UNSTRUCTURED_GRID'
c ... gerando os nos
      call makeNode(neq,nout)
      call makeConect(neq,nout)
      i_p = alloc_4('mshape  ',1,neq*neq) 
      call mzero(ia(i_p),neq*neq) 
      call makeShape(ip,ja,neq,ia(i_p),nout)
      return
      end
c *********************************************************************
c
c *********************************************************************
c * MAKENODE: gera os nos para visualizacao da matriz em um grid de   *
c * pontos no formato vtk                                             *
c * ------------------------------------------------------------------*
c * Parametros de entrada:                                            *
c * ------------------------------------------------------------------*
c * neq  - numero de equacoes                                         *
c * nout - Arquivo de saida                                           *
c * ------------------------------------------------------------------*
c * Parmetros de saida:                                               *
c * ------------------------------------------------------------------*
c *********************************************************************
      subroutine makeNode(neq,nout)
      implicit none
      integer neq,nout
      real*8 x,y
      integer line,col,i
c ... coordenadas X:dx=1/neq -> no(line,col)=((line-1)*dx,(col-1))*dy)
      i = 1
      write(nout,'(a,i10,a)') 'POINTS ',neq*(1+neq)/2,' double'
      do line = 1, neq
        y =  -(dfloat(line) -1.0) / neq
        do col = line, neq
          x =   (dfloat(col) + 1.0) /neq
          write(nout,*),x,y,0.0
          i = i + 1
        enddo  
      enddo
      return
      end
c *********************************************************************
c
c *********************************************************************
c * MAKENODE: gera os conectividades para visualizacao da matriz      *
c * em um grid de pontos no formato vtk                               *
c * ------------------------------------------------------------------*
c * Parametros de entrada:                                            *
c * ------------------------------------------------------------------*
c * neq  - numero de equacoes                                         *
c * nout - Arquivo de saida                                           *
c * ------------------------------------------------------------------*
c * Parmetros de saida:                                               *
c * ------------------------------------------------------------------*
c *********************************************************************
      subroutine makeConect(neq,nout)
      implicit none
      integer neq,nout,dl
      integer line,col,i
      write(nout,'(a,i10,i10)') 'CELLS ',neq*(neq+1)/2,neq*(neq+1)  
      i = 1
      do line = 1, neq 
        do col = 1, line 
          write(nout,*),1,i-1
          i = i + 1
        enddo  
      enddo
      write(nout,'(a,i10)') 'CELL_TYPES ',neq*(neq+1)/2
      do line = 1, neq*(neq+1)/2 
        write(nout,'(i3)') 1
      enddo
      return
      end
c *********************************************************************
c
c *********************************************************************
c * MAKESHAPE: identifica a diagonal principal e os elementos nao nulo*
c * da parte inferior da matriz                                       *
c * ------------------------------------------------------------------*
c * Parametros de entrada:                                            *
c * ------------------------------------------------------------------*
c * mshape -                                                          *
c * ip     - apontado da linha do csr                                 *
c * la     - colunas do csr                                           *
c * neq    - numero de equacoes                                       *
c * nout   - Arquivo de saida                                         *
c * ------------------------------------------------------------------*
c * Parmetros de saida:                                               *
c * ------------------------------------------------------------------*
c *********************************************************************
      subroutine makeShape(ia,ja,neq,mshape,nout)
      implicit none
      integer neq,nad,nout
      integer line,col,nd,nel,i,j
      integer mshape(*)
      integer ia(*),ja(*)
c ... diag principal
      i = 1
      do line = 1 , neq
        do col = line, neq
          if( line .eq. col ) mshape(i) = 1
          i = 1 + i 
        enddo
      enddo
c .....................................................................
c
c ... parte inferior da matriz
c
c ... CSR
      j = 0
      do line = 1, neq
        do i = ia(line), ia(line+1) -1
          nel = j + ja(i)
          if (ja(i) .gt. line) then
            mshape(nel) = 2
          endif
        enddo
        j = j + neq-line
      enddo       
c      j  = 1
c      do line = 2 , neq
c        do i = ia(line), ia(line+1) -1
c          nel = j + ja(i) 
c          mshape(nel) = 2
c        enddo
c        j = j + line
c      enddo
c .....................................................................
c
c ... saida vtk
      write(nout,'(a,i10)') 'CELL_DATA ',neq*(neq+1)/2  
      write(nout,'(a)') 'SCALARS perfil int '  
      write(nout,'(a)') 'LOOKUP_TABLE default'
      do line = 1, neq*(neq+1)/2
        write(nout,*) mshape(line) 
      enddo
c .....................................................................
c
c ...
      return
      end
c .....................................................................
c *********************************************************************

