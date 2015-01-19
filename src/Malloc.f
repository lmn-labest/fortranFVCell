c*****************************Svn***************************************      
c*$Date: 2011-09-18 16:28:54 -0300 (Sun, 18 Sep 2011) $                 
c*$Rev: 945 $                                                           
c*$Author: henrique $                                                   
c***********************************************************************    
c **********************************************************************
c *                                                                    *
c *   MALLOC.F                                            31/08/2005   *
c *                                                                    *
c *   Este arquivo contem subrotinas para gerenciamento de memoria:    *
c *                                                                    *
c *   init_malloc                                                      *
c *   alloc_4                                                          *
c *   alloc_8                                                          *
c *   locate                                                           *
c *   dealloc                                                          *
c *   icopy                                                            *
c *   azero                                                            *
c *   mzero                                                            *
c *   maxtest                                                          *
c *   mapmalloc                                                        * 
c *   use_work_vector                                                  *
c *   avaiblemem                                                       *
c *                                                                    *
c **********************************************************************
      module Malloc
         integer, allocatable, dimension(:) :: ia
         integer*8, external :: alloc_4,alloc_8,locate,dealloc
c         integer*8, parameter :: maxmem =1800000000
         integer*8 maxmem 
         data maxmem / 1200000000/
      end module
      subroutine init_malloc()
c **********************************************************************
c *                                                                    *
c *   INIT_MALLOC: inicializa a estrutura de dados para as rotinas     *
c *   ------------ de gerenciamento de memoria                         *
c *                                                                    *
c *   Variaveis do COMMON/MALLOC/:                                     *
c *   ---------------------------                                      *
c *                                                                    *
c *   ip(maxnpts) - ponteiro do i-esimo arranjo                        *
c *   arname(maxnpts)   - nome do i-esimo arranjo                      *
c *   nalp        - numero de ponteiros alocados                       *
c *   ip(nalp+1)  - proximo ponteiro livre                             *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
      integer maxnpts
      parameter (maxnpts = 200)
      integer i,ierr
      integer nalp,align
c ... valores do ponteiros     
      integer*8 ip(maxnpts)
      character*8 arname(maxnpts)
      common /malloc_info/ arname,ip,nalp,align
c ...
      align =   16
c ...
      allocate(ia(maxmem), stat=ierr)
      if (ierr .ne. 0) then
         print *, 'error: cannot allocate memory pool in heap.'
         stop
      endif
c ......................................................................
      print*,"**************   init_malloc   ***********************" 
      print*,"Available mem in work vector.",(maxmem*4)/(1024**2),"MBs"
      print*,"******************************************************" 
c ......................................................................
      do i = 1, maxnpts
         ip(i)     = 0
         arname(i) = '        '
      enddo
      nalp = 0
      ip(nalp+1) = 1
      return
      end
      integer*8 function alloc_4(name,nl,nc)
c **********************************************************************
c *                                                                    *
c *   ALLOC_4: aloca memoria para arranjos de 4 bytes                  *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *   name - nome do arranjo a ser alocado                             *
c *   nl   - numero de linhas do arranjo                               *
c *   nc   - numero de colunas do arranjo                              *
c *                      
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   alloc_4 - ponteiro do arranjo                                    *
c *                                                                    *
c **********************************************************************
      implicit none
      integer maxnpts
      parameter (maxnpts = 200)
c ... valores do ponteiros     
      integer*8 ip(maxnpts)
      integer*8 ipi,n,locate
c ......................................................................
      character*8 arname(maxnpts),name
      integer nc,nalp,align
      integer nl
      common /malloc_info/ arname,ip,nalp,align
c ......................................................................
      n = nl*nc
c ... alinhamento com 16 bytes
c     (libera a memoria em pacotes de 16, 32 , 48, ... bytes)
      if(align .eq. 16) then
        n = 4 * n
        if(mod(n,16) .ne. 0) then
          n = (1 + n/16)*16
        endif
        n = n / 4
      endif
c ......................................................................
c     if (n .le. 0) then
c        aloc4 = ip(nalp+1)
c        return
c     endif
      if (n .le. 0) then
         print*,'*** <ALLOC_4> numero  negativo de posicoes no vetor:',
     .      '(',name,') ***'
c        call stop_mef()
         stop
      endif
c ......................................................................      
      ipi = locate(name)
      if(ipi .gt. 0) then
         print*,'*** <ALLOC_4> Nome de arranjo existente: ',
     .      '(',name,') ***'
c         call stop_mef()
         stop
      endif
c ......................................................................      
      nalp   = nalp + 1
      if (nalp+1 .gt. maxnpts) then
         print*,'*** <ALLOC_4> Max. numero de ponteiros excedido: ',
     .        '(',name,') ***'
c         call stop_mef()
         stop
      endif
c ......................................................................      
      ipi = ip(nalp)
      ip(nalp+1) = ipi + n
      call maxtest(ip(nalp+1),name)
      arname(nalp) = name
      alloc_4 = ipi
c ......................................................................                  
      return
      end
      integer*8 function alloc_8(name,nl,nc)
c **********************************************************************
c *                                                                    *
c *   ALLOC_8: aloca memoria para arranjos de 8 bytes                  *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *   name - nome do arranjo a ser alocado                             *
c *   nl   - numero de linhas do arranjo                               *
c *   nc   - numero de colunas                                         *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   alloc_8 - ponteiro do arranjo                                    *
c *                                                                    *
c **********************************************************************
      implicit none 
      integer maxnpts
      parameter  (maxnpts = 200)
c ... valores do ponteiros     
      integer*8 ip(maxnpts)
      integer*8 locate,ipi,n
c ......................................................................
      character*8 arname(maxnpts),name
      integer nl,nc,nalp,align
      common /malloc_info/ arname,ip,nalp,align
c ......................................................................
      n = nl*nc
c ... alinhamento com 16 bytes
c     (libera a memoria em pacotes de 16, 32 , 48, ... bytes)
      if(align .eq. 16) then
        n = 8 * n
        if(mod(n,16) .ne. 0) then
          n = (1 + n/16)*16
        endif
        n = n / 4
      endif
c     if (n .le. 0) then
c        aloc8 = ip(nalp+1)
c        return
c     endif      
      if (n .le. 0) then
         print*,'*** <ALLOC_8> numero negativo de posicoes no vetor:',
     .      '(',name,') ***'
c         call stop_mef()
          stop
      endif
c ......................................................................
      ipi = locate(name)
      if(ipi .gt. 0) then
         print*,'*** <ALLOC_8> Nome de arranjo existente: ',
     .        '(',name,') ***'
c         call stop_mef()
         stop
      endif
c ......................................................................            
      nalp = nalp + 1
      if (nalp+1 .gt. maxnpts) then
         print*,'*** <ALLOC_8> Max. numero de ponteiros excedido: ',
     .        '(',name,') ***'
c         call stop_mef()
         stop
      endif
c ......................................................................                  
       ipi = ip(nalp)
       if(align .eq. 0) then 
         if(mod(ipi,2) .eq. 0) then
           ipi      = ipi+1
           ip(nalp) = ipi
         endif
         n = n*2 
      endif
      ip(nalp+1) = ipi + n
      call maxtest(ip(nalp+1),name)
      arname(nalp) = name
      alloc_8 = ipi
c ......................................................................                  
      return
      end
      integer*8 function locate(name)
c **********************************************************************
c *                                                                    *
c *   LOCATE: localiza o ponteiro do arranjo 'name'                    *
c *   ------                                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *   name - nome do arranjo a ser localizado                          *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   locate - ponteiro do arranjo                                     *
c *                                                                    *
c **********************************************************************
      implicit none
      integer maxnpts
      parameter (maxnpts = 200)
c ... valores do ponteiros     
      integer*8 ip(maxnpts)
c ......................................................................
      character*8 arname(maxnpts),name
      integer i,nalp
      common /malloc_info/ arname,ip,nalp
c ......................................................................
      locate = 0      
      do i = 1, nalp
         if(name .eq. arname(i)) then
               locate = ip(i)
             return
         endif
      enddo
      return
      end
      integer*8 function dealloc(name)
c **********************************************************************
c *                                                                    *
c *   Desaloca o arranjo com o nome name, e recalcula os ponteiro      *
c *   OBS: E necessario recupera os ponteiros explicitamente           *
c *                                                                    *
c **********************************************************************
      implicit none
      integer maxnpts
      parameter (maxnpts = 200)
c ... valores do ponteiros     
      integer*8 ip(maxnpts),ip1
      integer*8 npos,npos0
c ......................................................................      
      character*8 arname(maxnpts),name
      integer nalp,i,j
      common /malloc_info/ arname,ip,nalp
c ......................................................................
c
c ... Localiza a posicao i do ponteiro a ser removido:
c
      dealloc = 0
      do 100 i = 1, nalp
         if(name .eq. arname(i)) then
            goto 200
         endif
  100 continue
      print*,'*** <DEALLOC> Ponteiro nao encontrado: ','(',name,') ***'
      stop
c      call stop_mef()
  200 continue
c
c ... Remove o arranjo da memoria:
c
      ip1   = ip(i+1)
      npos0 = ip1 - ip(i)
      npos  = npos0 - mod(npos0,2)
      do 210 j  = i, nalp-1
         ip(j)  = ip(j+1) - npos
         arname(j) = arname(j+1)
  210 continue
      call icopy(ip1,ip(nalp+1),ip(i))
      ip(nalp)   = ip(nalp+1) - npos
      ip(nalp+1) = 0
      arname(nalp)   = arname(nalp+1)
      arname(nalp+1) = '        '
      nalp = nalp - 1
      return
      end
      subroutine icopy(i1,i2,i3)
c **********************************************************************
c *                                                                    *
c *   MCOPY: desloca um vetor de inteiros na memoria.                  *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
c ... valores do ponteiros     
      integer*8 i1,i2,i3
      integer*8 i,j
c ...........................................
      j = i3
      do 100 i = i1, i2-1
         ia(j) = ia(i)
         j = j + 1
  100 continue
      return
      end
      subroutine azero(a,j)
c **********************************************************************
c *                                                                    *
c *   AZERO: zera as posicoes de 1 ate j do vetor a.                   *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8 a(*)
      integer k,j
c ......................................................................
      do 100 k = 1, j
         a(k) = 0.d0
  100 continue
      return
      end
      subroutine mzero(m,j)
c **********************************************************************
c *                                                                    *
c *   MZERO: zera as posicoes de 1 ate j do vetor m.                   *
c *                                                                    *
c **********************************************************************
      implicit none
      integer m(*)
      integer j,k
c ......................................................................
      do 100 k = 1, j
         m(k) = 0
  100 continue
      return
      end
      subroutine maxtest(ip,name)
c **********************************************************************
c *                                                                    *
c *   Verifica se ip > maxmem                                          *
c *                                                                    *
c **********************************************************************
      use Malloc
      implicit none
c ... valores do ponteiros     
      integer*8 ip
c ......................................................................      
      character*8 name
c ......................................................................
      if (ip .gt. maxmem) then
         print*,'*** <MAXTEST> Memoria insuficiente: (',name,') ***'
c
         ip = ip - maxmem
         print*, ip,' posicoes'
         ip = ip/1024/1024
         print*, ip,' MB'
c
c         call stop_mef()
         stop
      endif
      return
      end
c **********************************************************************
c
c **********************************************************************
c *                                                                    *
c *   Desaloca a memoria do sistema                                    *
c *                                                                    *
c **********************************************************************
      subroutine common_finalize()
      use Malloc
      implicit none
      integer ierr
c ...      
      if (allocated(ia)) then
         deallocate(ia, stat=ierr)
         if (ierr .ne. 0) then
            print *, 'error: cannot deallocate memory pool in heap.'
c            call stop_mef()
            stop
         endif
c .....................................................................
c
c ...
      else
         print *, 'error: trying to deallocate unallocated array.'
      endif
c .....................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
c * MAPMALLOC : gera o mapa do arranjos alocados                       *
c * ------------------------------------------------------------------ *
c * Parametro de entrada :                                             *
c * tp  - " B" bytes,"KB" kbytes,"MB" mbytes,"GB" gbytes               *
c * ------------------------------------------------------------------ *
c **********************************************************************
      subroutine mapmalloc(tp)
      use Malloc
      implicit none
      integer maxnpts
      parameter (maxnpts = 200)
c ... valores do ponteiros     
      integer*8 ip(maxnpts)
      character*8 arname(maxnpts)
c ... numero de ponteiros      
      integer nalp
c ...      
      common /malloc_info/ arname,ip,nalp
c ...      
      integer i
      integer*8 size
      real*8    bsize
      character*2 tp
c .....................................................................
c                                                                      
c ...
c ... vetor alloca para a calculo do tamanho do ultimo arranjo      
      write(*,*)"___________________________________________________"
      write(*,'(a,i15)')"***Mapa de arranjos alocados",nalp
      write(*,*)"___________________________________________________"
      write(*,100)"nome","ip","n ponteiro","Size" 
      do i = 1 , nalp
        size = ip(i+1) - ip(i)
        if(tp .eq. " B") then
          bsize = 4.0*size 
          write(*,200)arname(i),i,ip(i),bsize,"Bs"
        endif  
        if(tp .eq. "KB") then
          bsize = 4.0*size/1024.0 
          write(*,200)arname(i),i,ip(i),bsize,"KBs"
        endif  
        if(tp .eq. "MB") then
          bsize = 4.0*size/(1024.0**2) 
          write(*,200)arname(i),i,ip(i),bsize,"MBs"
        endif  
        if(tp .eq. "GB") then
          bsize = 4.0*size/(1024.0**3) 
          write(*,200)arname(i),i,ip(i),bsize,"GBs"
        endif  
      enddo
      write(*,*)"___________________________________________________"
c .....................................................................
c                                                                      
c ...
100   format('|',a4,4x,'|',3x,a2,'|',5x,a10,'|',9x,a10,'|')
200   format('|',a8,'|',i5,'|',i15,'|',f15.6,a4,'|')
      return
      end
c .....................................................................
c **********************************************************************
c
c **********************************************************************
c * USE_WORK_VECTOR: quantidade de memoria usada no vetor de trabalho  *
c * ------------------------------------------------------------------ *
c * Parametro de entrada :                                             *
c * cod - 1 Bytes, 2 KBytes, 3 MBytes, 4 GBytes                        *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ *
c **********************************************************************
      real*8 function use_work_vector(tp)
      use Malloc
      implicit none
      integer maxnpts
      parameter (maxnpts = 200)
c ... valores do ponteiros     
      integer*8 ip(maxnpts)
      character*8 arname(maxnpts)
c ... numero de ponteiros      
      integer nalp
c ...      
      common /malloc_info/ arname,ip,nalp
c ...      
      character*2 tp
c ...
      if(tp .eq.' B') then
        use_work_vector = ip(nalp + 1) * 4
      else if (tp .eq. 'KB') then
        use_work_vector = ip(nalp + 1) * 4 / 1024
      else if (tp .eq. 'MB') then
        use_work_vector = (ip(nalp + 1) * 4) / (1024**2)
      else if (tp .eq. 'GB') then
        use_work_vector = (ip(nalp + 1) * 4) / (1024**3)
      endif
c ......................................................................
c
c ...
      return
      end
c ......................................................................
c **********************************************************************
