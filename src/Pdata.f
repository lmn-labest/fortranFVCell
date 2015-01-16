c *********************************************************************
c * READPNODE : le o arquivo auxiliar com os nos que terao alguma     *
c * de suas grandeza impressas no tempo                               *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * fname   - nome do arquivo de entrada                              *
c * i_no    -                                                         *
c * i_nfile -                                                         *
c * num_node-                                                         *
c * nout    - arquivo de entrada                                      *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * i_no    -ponteiros para os nos                                    *
c * i_nfile -ponteiros para os arquivos                               *
c * num_node-numero de nos                                            *
c * flag    -verifica sucesso na abertura do arquivo                  *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine readpnode(fname,i_no,i_nfile,num_node,flag,nout)
      use Malloc
      implicit none
      include 'string.fi'
      integer*8 i_no,i_nfile
      integer num_node,i,j,no
      integer nout
      logical flag
      character*80 fname
      character*30 string
c .....................................................................
c
c ... 
      open(nout,file=fname,status= 'old' , err=1000 , action='read')
      flag = .true.
c ... numero total de nos      
      call readmacro(nout,.true.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*) num_node
      print*,num_node,' print set nodes'
c .....................................................................
c
c ... memoria para os nos
      i_no        = alloc_4('no      ',1,num_node)
c ... memoria para os arquivos      
      i_nfile     = alloc_4('nnodew  ',1,num_node)
c ... lendo nos
      do i = 1, num_node
         call readmacro(nout,.true.)
         write(string,'(30a)') (word(j),j=1,30)
         read(string,*) no
         ia(i_no+i-1) = no
         ia(i_nfile+i-1) = 50 + (i-1)
      enddo
c .....................................................................
      close(nout)
      return
 1000 continue
      print*, '*** Erro na abertura de arquivos: ',trim(fname)
      flag = .false.
      end
c *********************************************************************
c
c *********************************************************************
c *  PRINTNODE: imprime uma grandeza do no pelo tempo                 *
c *  ---------------------------------------------------------------- *
c *  Parametro de Entrada :                                           *
c *  ---------------------------------------------------------------- *
c *  u       - vetor com a grandeza a ser escrita                     *
c *  no      - numero do no                                           *
c *  ndf     - graus de liberdade da grandeza                         *
c *  istep   - passo de tempo                                         *
c *  dt      - intercalo de tempo                                     *
c *  nameres - nome da gradeza a ser escrita                          *
c *  prename - nome do arquivo de saida                               *
c *  nout    - numero do arquivo a ser escrito                        *
c *  open    - abertura de um novo arquivo .true.                     *
c *  ---------------------------------------------------------------- *
c *  Parametro de Saida :                                             *
c *********************************************************************
      subroutine printnode(u,no,ndf,istep,dt,nameres,prename,nout,code
     .                    ,open)
      implicit none
c ===        
      character*30 nameres
      real*8 u(ndf,*)
      real*8 dt
      integer no,istep,nout,ndf,i,code
      character*80 fileout,prename,name
      logical open
c =====================================================================        
c
c ===
c ... abre um novo arquivo
      if(open) then
        fileout = name(prename,no,code)
        open(unit=nout,file=fileout)
        write(nout,'(a,a,a,i9)')
     .  '# Variacao ',trim(nameres),' no tempo no ponto',no
      else
c ... reabre o arquivo e add uma nova linha      
        fileout = name(prename,no,code)
        open(unit=nout,file=fileout,status ='old',access='append') 
      endif
c =====================================================================        
c
c === 
      write(nout,'(i9,f15.6,9f20.10)')istep,istep*dt,(u(i,no), i=1, ndf)
c      write(nout,*)istep,',',istep*dt,',',u(no)
      close(nout)
      return
      end
c =====================================================================        
c *********************************************************************
c
c *********************************************************************
c *  SETELEV : substitui a coordenada Z pela solucao                  *
c *  ---------------------------------------------------------------- *
c *  Parametro de Entrada :                                           *
c *  ---------------------------------------------------------------- *
c *  un      - solucao nodal                                          *
c *  x       - coordenada                                             *
c *  xe      - nao atualizada                                         *
c *  nnode   - numero de nos                                          *
c *  ndm     - numero de dimensoes                                    *
c *  scale   - fator de escala                                        *
c *  ndf     - graus de liberdade                                     *
c *  pndf    - graus de liberdade a ser atualizados                   *
c *  ---------------------------------------------------------------- *
c *  Parametro de Saida :                                             *
c *  ---------------------------------------------------------------- *
c *  xe      - coordenada atualizadas                                 *
c *  ---------------------------------------------------------------- *
c *********************************************************************
      subroutine elev(un,x,xe,nnode,ndm,scale,ndf,pndf)
      implicit none
c ===        
      real*8 un(ndf,*),x(ndm,*),xe(3,*)
      real*8 scale
      integer nnode,i,ndf,pndf,ndm
c =====================================================================        
c
c ===
      do i = 1, nnode
        xe(1,i) = x(1,i)          
        xe(2,i) = x(2,i)          
        xe(3,i) = scale*un(pndf,i)
      enddo
      return
      end
c =====================================================================        
c *********************************************************************
