c **********************************************************************
c Depende das informacao da malha do mefpar, utiliando as variaveis do *
c elementos.fi                                                         *  
c***********************************************************************
c
c **********************************************************************
c * COOR_VTK: escreve as coordenadas do nos                            *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * bvtk              - true formato binary false ascii                *
c * -------------------------------------------------------------------*
c * Prametros de saida:                                                *
c * -------------------------------------------------------------------*
c *--------------------------------------------------------------------*
c **********************************************************************
      subroutine coor_vtk(coor,nnode,ndm,bvtk,nfile)
      implicit none
      integer nnode,ndm,nfile
      integer i
      logical bvtk
      character buffer*1024,str1*10,lf*1
      Real*8 coor(ndm,*)
      Real*8 dum
      lf =char(10)
c ======================================================================
c
c ...
      if(bvtk)then
        write(str1(1:10),'(i10)') nnode
        buffer ='POINTS '//str1//' double'//lf
        write(nfile) trim(buffer)
c ... 
        if(ndm.eq.1)then
          dum = 0.d0
          do i = 1,nnode
            write(nfile)coor(1,i),dum,dum
          enddo
c ...
        elseif(ndm.eq.2) then
          dum = 0.d0
          do i=1,nnode
           write(nfile) coor(1,i),coor(2,i),dum
          enddo
c ...
        elseif(ndm.eq.3) then
          do i=1,nnode
             write(nfile)coor(1,i),coor(2,i),coor(3,i)
          enddo
        endif    
c ======================================================================
c .....................................................................     
      else
        write(nfile,'(a,i10,a)') 'POINTS ', nnode ,' double'
c ...    
        if(ndm.eq.1) then
          do i=1,nnode
            write(nfile,'(3e20.10)') coor(1,i),0.0,0.0
          enddo
c ======================================================================
c
c ...
        elseif(ndm.eq.2) then
          do i=1,nnode
           write(nfile,'(3e20.10)') coor(1,i),coor(2,i),0.0
          enddo
c ======================================================================
c
c ...
        elseif(ndm.eq.3) then
          do i=1,nnode
             write(nfile,'(3e20.10)') coor(1,i),coor(2,i),coor(3,i)
c             write(nfile,'(3f26.10)')   coor(1,i),coor(2,i),coor(3,i)
          enddo
c ======================================================================
        endif
      endif  
c ======================================================================
      return
      end
c ======================================================================
c
      subroutine elm_vtk(nos,numel,nen,bvtk,nfile)
c **********************************************************************
c * elm_vtk: escreve elementos nos formato vtk com os seus respectivos *
c *           materias                                                 *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * nos(nen+1,numel)  - conetividade nodai                             *
c * numel             - numero de elementos                            *
c * ndm               - numeros de dimensoes                           *
c * nen               - numero maximo de nos por elemento              *
c * nfile             - arquivo de saida                               *
c * nbar2             - numeros de barras                              *
c * ntria3            - numero de triangulos                           *
c * nquad4            - numero de quaddrilateros                       *
c * ntetra4           - numero de tetraedros                           *
c * nhexa8            - numero de hexaedros                            *
c * nquad8            - numero de quadrilateros quadraticos            *
c * bvtk              - true formato binary false ascii                *
c * -------------------------------------------------------------------*
c * Prametros de saida:                                                *
c * -------------------------------------------------------------------*
c *                                                                    *
c * -------------------------------------------------------------------*
c * OBS:                                                               *
c * -------------------------------------------------------------------*
c * Os prametors nbar2,ntria3,nquad4,ntetra4 e nhexa8 foram herdados do*
c *   mefpar                                                           *
c *--------------------------------------------------------------------*
c **********************************************************************
      implicit none
      include 'elementos.fi'
      integer numel,nen,nfile
      integer nnoel
      integer numet,nb2,nt3,nq4,nt4,nh8
      integer i,j
      integer nos(nen+1,*)
      character buffer*1024,lf*1,str1*10,str2*10 
      logical bvtk

      lf =char(10)
c ======================================================================
c
c ... Calculo do numero tipo de cada elemento 
      nb2 = nbar2(1)
      nt3 = ntria3(1) 
      nq4 = nquad4(1) 
      nt4 = ntetra4(1)
      nh8 = nhexa8(1)
      numet = 3*nb2 + 4*nt3 + 5*nq4 + 5*nt4 + 9*nh8 
c
c ... total de elemntos e tamnhanho dos dados totais
      if(bvtk)then
        write(str1(1:10),'(i10)')numel
        write(str2(1:10),'(i10)')numet
        buffer = lf//'CELLS '//str1//str2//lf
        write(nfile) trim(buffer)  
      else
        write(nfile,'(a,i10,i10)') 'CELLS ',numel,numet  
      endif 
c ... escrevendo a malha
c     
c ... nos dos elementos
c
c ......................................................................
      if (nbar2(1) .gt. 0) then
        nnoel = 2 
         do i = nbar2(2), nbar2(2) + nbar2(1)-1
           if(bvtk)then
              write(nfile) nnoel,(nos(j,i)-1,j=1,2)
           else
              write(nfile,'(10i10)') nnoel,(nos(j,i)-1,j=1,2)
           endif  
         enddo
       endif
c ......................................................................
c
c ...triangulo linear                                                   
      if (ntria3(1) .gt. 0) then
        nnoel = 3 
        do i = ntria3(2), ntria3(2) + ntria3(1)-1
          if(bvtk)then
            write(nfile) nnoel,(nos(j,i)-1,j=1,3)
          else
            write(nfile,'(10i10)') nnoel,(nos(j,i)-1,j=1,3)
          endif  
        enddo
      endif
c ......................................................................
c
c ...quadrilatero bilinear                                              
      if (nquad4(1) .gt. 0) then
        nnoel = 4 
        do i = nquad4(2), nquad4(2) + nquad4(1)-1
          if(bvtk)then
            write(nfile) nnoel,(nos(j,i)-1,j=1,4)
          else
            write(nfile,'(10i10)') nnoel,(nos(j,i)-1,j=1,4)
          endif  
        enddo
      endif
c ......................................................................
c
c ... tetraedro linear                                                  
      if (ntetra4(1) .gt. 0) then
        nnoel = 4 
        do i = ntetra4(2), ntetra4(2) + ntetra4(1)-1
          if(bvtk)then
            write(nfile) nnoel,(nos(j,i)-1,j=1,4)
          else
            write(nfile,'(10i10)') nnoel,(nos(j,i)-1,j=1,4)
          endif  
        enddo              
      endif
c ......................................................................
c
c ... hexaedro linear                                                   
      if (nhexa8(1) .gt. 0) then
        nnoel = 8 
        do i = nhexa8(2), nhexa8(2) + nhexa8(1)-1
          if(bvtk)then
            write(nfile) nnoel,(nos(j,i)-1,j=1,8)
          else
            write(nfile,'(10i10)') nnoel,(nos(j,i)-1,j=1,8)
          endif  
        enddo 
      endif
c ......................................................................
c ======================================================================
c ... tipo dos elementos
c
      if(bvtk)then
        write(str1(1:10),'(i10)')numel
        buffer = lf//'CELL_TYPES '//str1//lf
        write(nfile) trim(buffer)
      else
        write(nfile,'(a,i10)') 'CELL_TYPES ',numel
      endif  
c
c ...     
      if (nbar2(1) .gt. 0) then
        nnoel = 3
        do i = nbar2(2), nbar2(2) + nbar2(1)-1
          if(bvtk)then
            write(nfile)nnoel
          else  
            write(nfile,'(i3)')nnoel
          endif  
        enddo
      endif
c .......................................................................
c
c ...
      if (ntria3(1) .gt. 0) then
        nnoel = 5
        do i = ntria3(2), ntria3(2) + ntria3(1)-1
          if(bvtk)then
            write(nfile)nnoel
          else  
            write(nfile,'(i3)')nnoel
          endif  
        enddo
      endif
c ......................................................................
c
c ...
      if (nquad4(1) .gt. 0) then
        nnoel = 9
        do i = nquad4(2), nquad4(2) + nquad4(1)-1
          if(bvtk)then
            write(nfile)nnoel
          else  
            write(nfile,'(i3)')nnoel
          endif  
        enddo
      endif
c ......................................................................
c
c ...
      if (ntetra4(1) .gt. 0) then
        nnoel = 10
        do i = ntetra4(2), ntetra4(2) + ntetra4(1)-1
          if(bvtk)then
            write(nfile)nnoel
          else  
            write(nfile,'(i3)') nnoel
          endif  
        enddo
      endif
c ......................................................................
c
c ...
      if (nhexa8(1) .gt. 0) then
        nnoel = 12
        do i = nhexa8(2), nhexa8(2) + nhexa8(1)-1
          if(bvtk)then
            write(nfile)nnoel
          else  
            write(nfile,'(i3)')nnoel
          endif  
        enddo
      endif
c ......................................................................
      return
      end
c ======================================================================
c **********************************************************************
c
c **********************************************************************
c * res_vtk: escreve os resultados no formato vtk                      *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * res(ndf,nnode)    - campo escalar,vetorial ou tensorial            *
c * nnode             - numero de nos                                  *
c * ndm               - dimesao                                        *
c * gdl               - graus de liberdade                             *
c * nfile             - arquivo de saida                               *
c * cname             - nome da variavel                               *
c * cod1              - codigo de instrucao                            *
c *                   1 -> campo escalar                               *
c *                   2 -> campo vetorial                              *
c *                   3 -> campo tensorial                             *
c * cod2              - 1 - interio de  4 bytes                        *  
c *                   - 2 - real de 4 bytes                            *  
c *                   - 3 - real de 8 bytes                            *  
c * bvtk              - true BINARY vtk, false ASCII                   *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine pont_prop_vtk(iprop,fprop,dprop,nnode,cname,ndm,gdl
     .                        ,cod1,cod2,bvtk,nfile)
      implicit none
      integer nnode,ndm,gdl,cod1,cod2,nfile
      integer  i,j
      integer iprop(gdl,*),idum
      Real*4 fprop(gdl,*),fdum
      Real*8 dprop(gdl,*),ddum
      character buffer*1024,lf*1,cname*15
      logical bvtk
      lf = char(10)
c ======================================================================
c
c === BINARY
      if(bvtk)then
c ... campo escalar
        if(cod1.eq.1) then
c .. escalar int BINARY     
          if(cod2.eq.1)then
            buffer =' SCALARS '//cname//' int'//lf
            write(nfile)trim(buffer)
            buffer =' LOOKUP_TABLE '//'default '//lf
            write(nfile)trim(buffer)
            do i=1,nnode
              write(nfile)iprop(1,i)
            enddo
c ... escalar float BINARY         
          elseif(cod2.eq.2)then  
            buffer = lf//'SCALARS '//cname//' float'//lf
            write(nfile)trim(buffer)
            buffer = lf//'LOOKUP_TABLE '//'default '//lf
            write(nfile)trim(buffer)
            do i=1,nnode
              write(nfile)fprop(1,i)
            enddo
c ... escalardouble BINARY         
          elseif(cod2.eq.3)then  
            buffer = lf//'SCALARS '//cname//' double'//lf
            write(nfile)trim(buffer)
            buffer = lf//'LOOKUP_TABLE '//'default '//lf
            write(nfile)trim(buffer)
            do i=1,nnode
              write(nfile)dprop(1,i)
            enddo
          endif  
c ......................................................................
c
c ... campo vetorial BINARY
        elseif(cod1.eq.2) then
c .. escalar int BINARY     
          if(cod2.eq.1)then
            buffer =' VECTORS '//cname//' int'//lf
            write(nfile)trim(buffer)
            do i=1,nnode
              if (gdl .eq. 2)then
                idum = 0
                write(nfile)(iprop(j,i),j=1,gdl),idum
              endif
              if (gdl .eq. 3)then
                write(nfile)(iprop(j,i),j=1,gdl)
              endif
            enddo
c ... escalar float BINARY         
          elseif(cod2.eq.2)then  
            buffer =' VECTORS '//cname//' float'//lf
            write(nfile)trim(buffer)
            do i=1,nnode
              if (gdl .eq. 2)then
                fdum = 0.0
                write(nfile)(fprop(j,i),j=1,gdl),fdum
              endif
              if (gdl .eq. 3)then
                write(nfile)(fprop(j,i),j=1,gdl)
              endif
            enddo
c ... escalardouble BINARY         
          elseif(cod2.eq.3)then  
            buffer =' VECTORS '//cname//' double'//lf
            write(nfile)trim(buffer)
            if (gdl .eq. 2)then
              do i=1,nnode
                ddum = 0.0
                write(nfile)(dprop(j,i),j=1,gdl),ddum
              enddo  
            else if (gdl .eq. 3)then
               do i=1,nnode
                write(nfile)(dprop(j,i),j=1,gdl)
               enddo 
            endif
          endif  
c ......................................................................
c 
c ... campo tensorial BINARY
        elseif(cod1.eq.3) then
          print*,'\nnao implementado ',cod1
        endif
c ......................................................................
c
c ======================================================================
c
c === ASCII
      else
c ... campo escalar ASCII
        if(cod1.eq.1) then
c .. escalar int ASCII 
          if(cod2.eq.1)then
            write(nfile,'(a,15a,a)')'SCALARS ',cname, 'int'
            write(nfile,'(a,a)')'LOOKUP_TABLE ','DEFAULT '
            do i=1,nnode
              write(nfile,'(i10)') iprop(1,i)
            enddo
c .. escalar float ASCII 
          elseif(cod2.eq.2)then  
            write(nfile,'(a,15a,a)')'SCALARS ',cname,' float'
            write(nfile,'(a,a)')'LOOKUP_TABLE ','default '
            do i=1,nnode
              write(nfile,'(7e15.5e3)') fprop(1,i)
            enddo
c .. escalar double ASCII 
          elseif(cod2.eq.3)then  
            write(nfile,'(a,15a,a)')'SCALARS ',cname,' double'
            write(nfile,'(a,a)')'LOOKUP_TABLE ','default '
            do i=1,nnode
              write(nfile,'(7e15.5e3)') dprop(1,i)
            enddo
          endif  
c ......................................................................
c
c ... campo vetorial ASCII
        elseif(cod1.eq.2) then
c .. escalar int ASCII 
          if(cod2.eq.1)then
            write(nfile,'(a,15a,a)')'VECTORS ',cname, 'int'
            do i=1,nnode
              if(gdl.eq.2)then
                write(nfile,'(i10)') (iprop(j,i),j=1,gdl),0
              elseif(gdl.eq.3)then
                write(nfile,'(7e15.5e3)') (fprop(j,i),j=1,gdl)
              endif  
            enddo
c .. escalar float ASCII 
          elseif(cod2.eq.2)then  
            write(nfile,'(a,15a,a)')'VECTORS ',cname, 'foalt'
            do i=1,nnode
              if(gdl.eq.2)then
                write(nfile,'(7e15.5e3)') (fprop(j,i),j=1,gdl),0.d0
              elseif(gdl.eq.3)then
                write(nfile,'(7e15.5e3)') (fprop(j,i),j=1,gdl)
              endif  
            enddo
c .. escalar double ASCII 
          elseif(cod2.eq.3)then  
            write(nfile,'(a,15a,a)')'VECTORS ',cname, 'double'
            if(gdl.eq.2)then
              do i=1,nnode
                write(nfile,'(3e20.10)') (dprop(j,i),j=1,gdl),0.d0
              enddo
            elseif(gdl.eq.3)then
              do i=1,nnode
                write(nfile,'(3e20.10)') (dprop(j,i),j=1,gdl)
              enddo
            endif  
          endif  
c ......................................................................

c ... campo tensorial ASCII
        elseif(cod1.eq.3) then
          print*,'\nnao implementado '
c ......................................................................
        endif
      endif  
c ======================================================================
c
      return
      end
c ======================================================================
c **********************************************************************
c
c ======================================================================
c **********************************************************************
c * cell_prop_vtk : escreve propriedades das celulas                   *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * dprop - propriedas real 8 bytes                                    *
c * iprop - propriedas interia de 8 bytes                              *
c * fprop - propriedas real de 4 bytes                                 *
c * cod1  - 1 - escalar                                                *  
c *       - 2 - vetorial                                               *  
c *       - 3 - tensorial                                              *  
c * cod2  - 1 - interio de  4 bytes                                    *  
c *       - 2 - real de 4 bytes                                        *  
c *       - 3 - real de 8 bytes                                        *  
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine cell_prop_vtk(iprop,fprop,dprop,numel,cname,ndm,gdl
     .                        ,cod1,cod2,bvtk,nfile)
c ===
      implicit none
      integer numel,ndm,gdl
      integer iprop(gdl,*),idum
      real*4 fprop(gdl,*),fdum
      real*8 dprop(gdl,*),ddum
      integer nfile
      logical bvtk
      integer cod1,cod2
      integer i,j
      character*15 cname
      character buffer*1024,lf*1,str1*10
      lf = char(10)
c ====================================================================== 
c === DADOS          
c ... VTk BINARY
      if(bvtk)then
        if(cod1.eq.1) then  
c ... int      
          if(cod2 .eq. 1)then
            write(str1(1:10),'(i10)') gdl    
            buffer = ' SCALARS '//cname//' int '//str1//lf
            write(nfile)trim(buffer)
            buffer = 'LOOKUP_TABLE default'//lf
            write(nfile)trim(buffer)
            do i = 1, numel
              write(nfile)(iprop(j,i),j=1,gdl)
            enddo
c ... float          
          else if(cod2 .eq. 2)then
            write(str1(1:10),'(i10)') gdl
            buffer = ' SCALARS '//cname//'float '//str1//lf
            write(nfile)trim(buffer)
            buffer = 'LOOKUP_TABLE default'//lf
            write(nfile)trim(buffer)
            do i = 1, numel
              write(nfile)(fprop(j,i),j=1,gdl)
            enddo
c ... double         
          else if(cod2 .eq. 3)then  
            write(str1(1:10),'(i10)') gdl
            buffer = ' SCALARS '//cname//' double '//str1//lf
            write(nfile)trim(buffer)
            buffer = 'LOOKUP_TABLE default'//lf
            write(nfile)trim(buffer)
            do i = 1, numel
              write(nfile)(dprop(j,i),j=1,gdl)
            enddo
          endif
c ... campo vetorial binario          
        else if(cod1.eq.2) then
c ... int                
          if(cod2.eq.1) then
c ... ifloat              
          elseif(cod2 .eq. 2) then
c ... double    
          elseif(cod2 .eq. 3) then
            buffer =' VECTORS '//cname//' double'//lf
            write(nfile)trim(buffer)
            if (gdl .eq. 2)then
              do i=1,numel
                ddum = 0.0
                write(nfile)(dprop(j,i),j=1,gdl),ddum
              enddo  
            else if (gdl .eq. 3)then
              do i=1,numel
                write(nfile)(dprop(j,i),j=1,gdl),ddum
              enddo 
            endif
          endif    
        endif
c .....................................................................
c          
c ... Vtk ASCII        
      else 
c ... int
        if(cod1.eq.1) then 
          if(cod2 .eq. 1)then
            write(nfile,'(a,1x,a15,1x,a5,1x,i3)')'SCALARS',cname,'int '
     .                   ,gdl
            write(nfile,'(a)')'LOOKUP_TABLE default'
            do i = 1, numel
c              write(nfile,'(i10)')iprop(1,i)
              write(nfile,'(99i10)')(iprop(j,i),j=1,gdl)
            enddo
c ... float          
          else if(cod2 .eq. 2)then  
            write(nfile,'(a,1x,15a,a)')'SCALARS',cname,'float'
            write(nfile,'(a)')'LOOKUP_TABLE default'
            do i = 1, numel
c              write(nfile,'(e20.10)')fprop(1,i)
              write(nfile,'(99f26.10)')(fprop(j,i),j=1,gdl)
            enddo
c ... double         
          else if(cod2 .eq. 3)then  
            write(nfile,'(a,1x,a15,1x,a8,1x,i3)')'SCALARS'
     .           ,cname,'double ',gdl
            write(nfile,'(a)')'LOOKUP_TABLE default'
            do i = 1, numel
c              write(nfile,'(es20.10)')dprop(1,i)
c              write(nfile,'(f26.10)')dprop(1,i)
               write(nfile,'(99f26.10)')(dprop(j,i),j=1,gdl) 
            enddo
          endif
c ... campo vetorial ASCII
        elseif(cod1.eq.2) then
c .. escalar int ASCII 
          if(cod2.eq.1)then
c .. escalar float ASCII 
          elseif(cod2.eq.2)then  
c .. escalar double ASCII 
          elseif(cod2.eq.3)then  
            write(nfile,'(a,15a,a)')'VECTORS ',cname, 'double'
            if(gdl.eq.2)then
              do i=1,numel 
                write(nfile,'(f20.6)') (dprop(j,i),j=1,gdl),0.d0
              enddo  
            elseif(gdl.eq.3)then
              do i=1,numel 
                write(nfile,'(f20.6)') (dprop(j,i),j=1,gdl)
              enddo
            endif  
          endif
        endif  
c ......................................................................
      endif  
c ====================================================================== 
c
c ===
      return
      end
c ====================================================================== 
c **********************************************************************
c
c **********************************************************************
c * head_vtk : escreve o vabecalho do arquivo vtk                      *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * headfile - cabecalho do arquivo vtk                                *
c * bvtk     - true formato binary false ascii                         *
c * t        - tempo real                                              *
c * istep    - passo de tempo                                          *
c * nfile    - arquivo de saida                                        *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine head_vtk(headfile,bvtk,t,istep,time,nfile)
c ===
      implicit none
      character*30 headfile
      character lf*1, buffer*1024
      integer nfile
      integer istep
      real*8 t
      logical bvtk,time
      lf =char(10)
c ====================================================================== 
c
c ===
c ...Cabecalho do arquivo Vtk
      if(bvtk)then
        buffer ='# vtk DataFile Version 3.0'//lf
        write(nfile) trim(buffer)
        write(nfile) headfile
        buffer =lf//'BINARY'//lf
        write(nfile) trim(buffer)
        buffer = 'DATASET UNSTRUCTURED_GRID'//lf
        write(nfile) trim(buffer)
        buffer = 'FIELD FieldData 2'//lf
        write(nfile) trim(buffer)
        buffer = 'TIME 1 1 double'//lf
        write(nfile) trim(buffer)
        write(nfile) t           
        buffer = 'CYCLE 1 1 int'//lf
        write(nfile) trim(buffer)
        write(nfile) istep 
      else
        write(nfile,'(a)') '# vtk DataFile Version 3.0'
        write(nfile,'(a)')headfile 
        write(nfile,'(a)') 'ASCII'
        write(nfile,'(a)') 'DATASET UNSTRUCTURED_GRID'
        if(time) then
          write(nfile,'(a)') 'FIELD FieldData 2'
          write(nfile,'(a)') 'TIME 1 1 double'
          write(nfile,'(f16.6)') t
          write(nfile,'(a)') 'CYCLE 1 1 int'
          write(nfile,'(i9)') istep
        endif  
      endif  
c ======================================================================
c
c ===
      return
      end
c ======================================================================
c
c **********************************************************************
c * CELL_DATA: incia um secao CELL_DATA                                *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * numel    - numero de elementos                                     *
c * bvtk     - true formato binary false ascii                         *
c * nfile    - arquivo de saida                                        *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine cell_data_vtk(numel,bvtk,nfile)
      implicit none
      integer numel
      integer nfile
      logical bvtk
      character str1*10,buffer*1024,lf*1
      lf = char(10)
c ... VTK BINARY
      if(bvtk)then
        write(str1(1:10),'(i10)')numel
        buffer = lf//'CELL_DATA '//str1//lf
        write(nfile) trim(buffer)
c ... VTK ASCII 
      else
        write(nfile,'(a,i10)')'CELL_DATA',numel
      endif
c........................................................................
c 
c ... 
      return
      end
c **********************************************************************
c       
c **********************************************************************
c * POINT_DATA: incia um secao POINT_DATA                              *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * nnode    - numero de pontos                                        *
c * bvtk     - true formato binary false ascii                         *
c * nfile    - arquivo de saida                                        *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine point_data_vtk(nnode,bvtk,nfile)
      implicit none
      integer nnode
      integer nfile
      logical bvtk
      character str1*10,buffer*1024,lf*1
      lf = char(10)
c ... VTK BINARY
      if(bvtk)then
        write(str1(1:10),'(i10)')nnode
        buffer = lf//' POINT_DATA '//str1//lf
        write(nfile) trim(buffer)
c ... VTK ASCII 
      else
        write(nfile,'(a,i10)')'POINT_DATA',nnode
      endif
c........................................................................
c 
c ... 
      return
      end
c **********************************************************************

