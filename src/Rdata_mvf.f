c **********************************************************************
c *                                                                    *
c *   RDAT: leitura de dados.                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    nin     - arquivo de entrada                                    *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    nnode - numero de nos                                           *
c *    numel - numero de elementos                                     *
c *    nen   - numero max. de nos por elemento                         *
c *    ndf   - numero max. de graus de liberdade por no                *
c *    ndm   - dimensao (1, 2 ou 3)                                    *
c *    i_ix    - ponteiro para conetividades                           *
c *    i_x     - ponteiro para o arranjo x                             *
c *                                                                    *
c **********************************************************************
      subroutine rdat_mvf(i_x,i_ix,i_pedgeT1,i_sedgeT1,i_pedgeF,i_sedgeF
     .                   ,i_pedgeE,i_sedgeE,i_pnode,i_snode
     .                   ,i_fnode,i_e,i_ie,i_w,i_t1,i_u1,i_u2,i_p,i_temp
     .                   ,i_num
     .                   ,nnode,numel,ndm,ndfF,ndfE,ndfT1,sTrans,sSimple
     .                   ,sEnergy,nen,numat,nin)
      use Malloc
      implicit none
      include 'elementos.fi'
      include 'string.fi'
      include 'simple.fi'
c .....................................................................      
      integer nnode,numel,ndm,ndfT1,ndfF,ndfE,nen,numat
      logical sTrans,sSimple,sEnergy
c ... ponteiros
      integer*8 i_x,i_ix
      integer*8 i_pedgeT1,i_sedgeT1,i_pedgeF,i_sedgeF,i_pedgeE,i_sedgeE
      integer*8  i_e,i_ie,i_w,i_t1,i_u1,i_u2,i_p,i_temp
      integer*8 i_pnode,i_snode,i_fnode,i_num
c .....................................................................
      integer nin,naux,nincl,totnel
      integer nmc,nmacros,j
      character*15 macro(30),rc
      character*80 fname
c .....................................................................
      data macro/'materials      ','bar2           ','tria3          ',
     .           'quad4          ','tetra4         ','hexa8          ',
     .           '               ','               ','               ',
     .           'coordinates    ','elconstrainT1  ','edgeT1         ',
     .           'constraintemp  ','nodaltemp      ','nodalflux      ',
     .           'elveloc        ','initialtemp    ','initialPhy     ',
     .           'elconstrainCfd ','edgeCfd        ','elconstrainE   ',
     .           'edgeE          ','initialE       ','elinitialE     ',
     .           'initialT1      ','insert         ','return         ',
     .           '               ','gravity        ','end            '/
      data nmc /30/
      data nincl /4/
c ......................................................................
c
c ... 
      naux = - 1
c ... Leitura dos parametros da malha: nnode,numel,numat,nen,ndf,ndm
      call parameters(nnode,numel,nen,ndfT1,ndfF,ndfE,ndm,numat,nin)
c ......................................................................      
c     Alocacao de arranjos na memoria:
c     ---------------------------------------------------------------
c                 
c     ---------------------------------------------------------------
c
      i_num   = alloc_4('num     ',    1,numel)
      i_ix    = alloc_4('ix      ',nen+1,numel)
      i_ie    = alloc_4('ie      ',    1,numat)
      i_x     = alloc_8('x       ',  ndm,nnode)
      i_w     = alloc_8('w       ',  ndm,numel)
      i_e     = alloc_8('e       ',   10,numat)
      call mzero(ia(i_num),numel)
      call mzero(ia(i_ix),numel*(nen+1))
      call mzero(ia(i_ie),numat)
      call azero(ia(i_x),nnode*ndm)
      call azero(ia(i_e),10*numat)
      call azero(ia(i_w),ndm*numel)
c ... transporte      
      if( ndfT1 .eq. 1) then
        sTrans = .true.      
        i_pedgeT1 = alloc_4('pedgeT1 ',nen+1,numel)
        i_pnode   = alloc_4('pnode   ',ndfT1,nnode)
        i_sedgeT1 = alloc_8('sedgeT1 ',nen+1,numel)
        i_snode   = alloc_8('snode   ',ndfT1,nnode)
        i_fnode   = alloc_8('fnode   ',ndfT1,nnode)
        i_t1      = alloc_8('t1      ',ndfT1,numel)
        call mzero(ia(i_pedgeT1),numel*(nen+1))
        call mzero(ia(i_pnode),nnode*ndfT1)
        call azero(ia(i_t1),ndfT1*numel)
        call azero(ia(i_sedgeT1),numel*(nen+1)*ndfT1)
        call azero(ia(i_snode),nnode*ndfT1)
        call azero(ia(i_fnode),nnode*ndfT1)
      endif
c ... cfd     
      if( ndfF .eq. 3) then
        sSimple = .true.         
        i_pedgeF    = alloc_4('pedgeF  ',nen+1,numel)
        i_sedgeF    = alloc_8('sedgeF  ',nen+1,numel*ndfF)
        i_u1        = alloc_8('u1      ',    1,numel)
        i_u2        = alloc_8('u2      ',    1,numel)
        i_p         = alloc_8('p       ',    1,numel)
        call mzero(ia(i_pedgeF),numel*(nen+1))
        call azero(ia(i_u1),numel)
        call azero(ia(i_u2),numel)
        call azero(ia(i_p),numel)
        call azero(ia(i_sedgeF),numel*(nen+1)*ndfF)
      endif
c ... Equacao de energia
      i_pedgeE = 1
      i_sedgeE = 1
      i_temp   = 1  
      if( ndfE .eq. 1) then
        sEnergy = .true.         
        i_pedgeE= alloc_4('pedgeE  ',nen+1,numel)
        i_sedgeE= alloc_8('sedgeE  ',nen+1,numel*ndfE)
        i_temp  = alloc_8('temp    ',    1,numel)
        call mzero(ia(i_pedgeE),numel*(nen+1))
        call azero(ia(i_temp),numel)
        call azero(ia(i_sedgeE),numel*(nen+1)*ndfE)
      endif         
c .....................................................................
      nmacros = 0
      totnel  = 0
c .....................................................................
  100 continue
      call readmacro(nin,.true.)
      write(rc,'(15a)') (word(j),j=1,15)
      do 200 j = 1, nmc
         if (rc .eq. macro(j)) go to 300
  200 continue
      go to 100
  300 continue
c ......................................................................      
c
c ... Controle de Macros (PRE-processador):
c
      nmacros = nmacros + 1
c ......................................................................      
      goto(400,450,500      !materials  , bar2          , tria3        ,
     .    ,550,450,450      !quad4      , tetra4        , hexa8        ,
     .    ,450,450,450      !           ,               ,              ,
     .    ,850,900,950      !coordinates,elconstraintemp,edgesources   ,
     .    ,1000,1050,1100   
     .    ,1150,1200,1250   !elveloc    ,initialtemp    ,initialPhy    ,
     .    ,1300,1350,1400   !elconstrainCfd,edgeCfd     ,elconstrainE  ,
     .    ,1450,1500,1550   !edgeE      ,initialE       ,elinitialE    ,
     .    ,1600,1650,1700   !initialT1  ,insert         ,return        ,
     .    ,1900,1950,2000)j !           ,gravity        ,end           ,
c ......................................................................
c
c ... indefinido                    
c
  400 continue
        print*,'Loading materials ...'
        call mate(ia(i_ie),ia(i_e),numat,nin)
        print*,'load'
      go to 100
c ......................................................................      
c
c ... indefinido                    
c
  450 continue
        print*,'Macro indefinida!!!'
        stop
      go to 100
c ......................................................................      
c
c ... Conetividades quad4:          
c
  500 continue
      print*,'Loading tria3 ...'
      ntria3(1) = 0
      call elconn(ia(i_ix),nen+1,3,ntria3(1),numel,nin)
      ntria3(2) = totnel + 1
      totnel    = totnel + ntria3(1)
      print*,'load'
      go to 100 
c ......................................................................      
c
c ... Conetividades quad4:          
c
  550 continue
      print*,'Loading quad4 ...'
      nquad4(1) = 0
      call elconn(ia(i_ix),nen+1,4,nquad4(1),numel,nin)
      nquad4(2) = totnel + 1
      totnel    = totnel + nquad4(1)
      print*,'load'
      go to 100 
c ......................................................................      
c
c ... Coordenadas                   
c
  850 continue
      print*,'Loading coordinates ...'
      call coord(ia(i_x),nnode,ndm,nin)
      print*,'load'
      go to 100 
c ......................................................................      
c
c ... elconstraint1                   
c
  900 continue
      print*,'Loading elconstraint1   ...'
      call pegde(ia(i_pedgeT1),numel,nen+1,nin)
      print*,'load'
      go to 100 
c ......................................................................      
c
c ... edgeT1                     
c
  950 continue
      print*,'Loading edgeT1 ...'
      call segde(ia(i_sedgeT1),numel,nen+1,nin)
      print*,'load'
      go to 100 
c ......................................................................      
c
c ... constraintemp                 
c
 1000 continue
      print*,'Loading constraintemp ...'
      call bound(ia(i_pnode),nnode,ndfT1,nin,1)
      print*,'load'
      go to 100 
c ......................................................................      
c
c ... nodaltemp                     
c
 1050 continue
      print*,'Loading nodaltemp ...'
      call forces(ia(i_snode),nnode,ndfT1,nin,1)
      print*,'load'
      go to 100 
c ......................................................................
c
c ... nodalflux                     
c
 1100 continue
      print*,'Loading nodalflux ...'
      call forces(ia(i_fnode),nnode,ndfT1,nin,1)
      print*,'load'
      go to 100 
c ......................................................................
c
c ... elveloc                       
c
 1150 continue
      print*,'Loading elveloc ...'
      call forces(ia(i_w),numel,ndm,nin)
      print*,'load'
      go to 100 
c ......................................................................
c
c ... condicao inicial                       
c
 1200 continue
      print*,'Loading initialtemp ...'
      call forces(ia(i_t1),numel,1,nin)
      print*,'load'
      go to 100 
c ......................................................................
c
c ... condicao inicial de pressao                      
c
 1250 continue
      print*,'Loading elinitialPhy ...'
      call PressureHy(ia(i_p),ia(i_x),g,ia(i_ix),ia(i_ie),ia(i_e)
     .               ,numel,ndm,nen)
      print*,'load'
      go to 100 
c ......................................................................
c
c ... elconstrainCfd                       
c
 1300 continue
      print*,'Loading elconstrainCfd ...'
      call pegde(ia(i_pedgeF),numel,(nen+1),nin)
      print*,'load'
      go to 100 
c ......................................................................
c
c ... edgeCfd                       
c
 1350 continue
      print*,'Loading edgeCfd ...'
      call segde(ia(i_sedgeF),numel,(nen+1)*ndfF,nin)
      print*,'load'
      go to 100 
c .....................................................................
c
c ... elconstrainE                         
c
 1400 continue
      if(sEnergy) then
        print*,'Loading elconstrainE   ...'
        call pegde(ia(i_pedgeE),numel,nen+1,nin)
        print*,'load'
      endif
      go to 100 
c ......................................................................
c
c ... edgeE                         
c
 1450 continue
      if(sEnergy) then
        print*,'Loading edgeE   ...'
        call segde(ia(i_sedgeE),numel,nen+1,nin)
        print*,'load'
      endif
      go to 100 
c .....................................................................
c
c ... initialE                      
c
 1500 continue
      if(sEnergy) then
        print*,'Loading initialE ...'
        call initialP(ia(i_temp),numel,ndfE,nin)
        print*,'load'
      endif
      go to 100 
c .....................................................................
c
c ... elinitialE                      
c
 1550 continue
      if(sEnergy) then
        print*,'Loading elinitialE ...'
        call forces(ia(i_temp),numel,1,nin)
        print*,'load'
      endif
      go to 100 
c .....................................................................
c
c ... initialT1                     
c
 1600 continue
      if(sTrans) then
        print*,'Loading initialT1 ...'
        call initialP(ia(i_t1),numel,ndfT1,nin)
        print*,'load'
      endif  
      go to 100 
c .....................................................................
c
c ... (insert) Desvia leitura para arquivo auxiliar:
c
 1650 continue
      nmacros = nmacros - 1
      naux = nin
      call readmacro(nin,.false.)
      write(fname,'(80a)') (word(j),j=1,strl)
      open(nincl, file= fname,status= 'old',err=1651,action='read')
      nin = nincl
      go to 100
 1651 continue
      print*, trim(fname), ' arquivo nao existente !'
      stop
c ......................................................................      
c
c ... (return) Retorna leitura para arquivo de dados basico:
c      
 1700 continue
      nmacros = nmacros - 1
      close(nincl)
      nin = naux
      go to 100
c ......................................................................
c
c ...                               
c
 1900 continue
      go to 100 
c .....................................................................
c
c ... gravity                       
c
 1950 continue
      print*,'Loading gravity  ...'
      call gravity(g,ndm,nin)
      print*,'load'
      go to 100 
c .....................................................................
c
c ... End  
 2000 continue
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine mate(ie,e,numat,nin)
c **********************************************************************
c *                                                                    *
c *   Materiais.                                                       *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
c     include 'termprop.fi'
      integer ie(*),numat,nin,i,j,m,ma
      real*8  e(10,*)
      character*30 string
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(j),j=1,12)
c .....................................................         
      do 110 while(string .ne. 'end')
         read(string,*,err = 200,end = 200) ma
         if(ma .lt. 1 .or. ma .gt. numat) goto 200
c .....................................................            
         call readmacro(nin,.false.)
         write(string,'(12a)') (word(m),m=1,12)
         read(string,*,err = 200,end = 200) ie(ma)
c .....................................................     
         call readmacro(nin,.false.)
         write(string,'(30a)') (word(m),m=1,30)
         i = 0
c
c ... linha de comando original:
c         do 100 while(string .ne. ' ')
c
c ... formato necessário ao funcionamento em linux:
         do 100 while ( (string .ne.   ' '    ) .and.
     .                  (string .ne. CHAR(13) ) )
            i = i + 1
            if(i .gt. 10) goto 200
            read(string,*,err = 200,end = 200) e(i,ma)
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(m),m=1,30)               
  100    continue
         call readmacro(nin,.true.)
         write(string,'(12a)') (word(j),j=1,12)
  110 continue
      return
c ......................................................................                        
  200 continue
      print*,'*** Erro na leitura dos materiais !'
      stop             
      end      
c **********************************************************************
c
c **********************************************************************
      subroutine initialP(u,n,ndf,nin)
c **********************************************************************
c *                                                                    *
c *   INITIAL: codicao inicia para todo o dominio                      *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      integer n,ndf,nin
      integer i,j,m
      real*8  u(ndf,n),k
      character*30 string
c ......................................................................
c
c ...      
      call readmacro(nin,.true.)
      write(string,'(30a)') (word(m),m=1,30)
      read(string,*,err = 200,end = 200) k
c ......................................................................
c
c ...      
      do i = 1, n
        do j = 1, ndf
          u(j,i) = k
        enddo
      enddo
c ......................................................................
c
c ...      
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(j),j=1,12)
      if (string .ne. 'end') goto 200
c ......................................................................
c
c ...
      return
c ......................................................................
  200 continue
      print*,'*** Erro na leitura (initialP) !'
      stop  
      end
c **********************************************************************
c
c **********************************************************************
      subroutine gravity(g,n,nin)
c **********************************************************************
c *                                                                    *
c *   INITIAL: codicao inicia para todo o dominio                      *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      integer n
      integer nin,m,j
      real*8  g(3)
      character*30 string
c ......................................................................
c
c ... 
      read(nin,*) g(1:n)
c ......................................................................
c
c ...      
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(j),j=1,12)
      if (string .ne. 'end') goto 200
c ......................................................................
c
c ...
      return
c ......................................................................
  200 continue
      print*,'*** Erro na leitura (gravity) !'
      stop  
      end
c **********************************************************************
c 
c **********************************************************************
      subroutine coord(x,nnode,ndm,nin)
c **********************************************************************
c *                                                                    *
c *   Coordenadas.                                                     *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      integer nnode,ndm,nin,j,k,m
      real*8 x(ndm,nnode)
      character*30 string
c ......................................................................
      do j = 1, nnode
         read(nin,*) k,(x(m,k),m=1,ndm)
      enddo
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(j),j=1,12)
      if (string .ne. 'end') goto 200
      return
c ......................................................................
  200 continue
      print*,'*** Erro na leitura das coordenadas !'
      stop  
      end
c **********************************************************************
c
c **********************************************************************
      subroutine elconn(ix,nen1,nen,nel,numel,nin)
c **********************************************************************
c *                                                                    *
c *   Conetividades nodais.                                            *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'      
      integer ix(nen1,*),nen1,nen,nel,numel,nin,j,k,m
      character*12 string
c ......................................................................
       do 100 j = 1, numel
          read(nin,*) k,(ix(m,k),m=1,nen1)
  100  continue
       call readmacro(nin,.true.)
       write(string,'(12a)') (word(j),j=1,12)  
       if (string .ne. 'end') goto 200
       nel=numel
       return
      print*,'*** Erro na leitura dos elementos !'
 200  continue
c ......................................................................
      stop
      end
c **********************************************************************
c
c **********************************************************************
      subroutine readmacro(nin,newline)
c **********************************************************************
c *                                                                    *
c *   Subroutine READMACRO: le uma macro numa linha nova ou a partir   *
c *                         da coluna line_col de uma linha lida       *
c *                         anteriormente e armazenada em line(*)      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   nin - numero do arquivo de entrada                               *
c *   newline = .true.  - nova linha deve ser lida                     *
c *           = .false. - macro deve ser lida a partir da coluna       *
c *                       line_col em line(*)                          *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      integer j,k,nin
      logical newline
c ......................................................................
      if(newline) then
         line_col = 1
         read(nin,'(200a1)',err = 200,end = 200) (line(j),j=1,200)
      endif
c ......................................................................      
      do j = 1, maxstrl
         word(j) = ' '
      enddo
      strl = 0
      if(line_col .gt. maxstrl) return
c ......................................................................      
      do while (line(line_col) .eq. ' ')
         line_col = line_col + 1
         if (line_col .gt. maxstrl) return         
      end do
c ......................................................................      
      do while ( (line(line_col) .ne. ' ') .and.
     .           (line(line_col) .ne. CHAR(13)) )
         strl = strl + 1
         line_col = line_col + 1
         if (line_col .gt. maxstrl) goto 100
      end do
c ......................................................................      
  100 continue      
      k = line_col-strl
      do j = 1, strl
         write(word(j),'(a)') line(k)
         k = k + 1
      end do      
      return
c ......................................................................      
  200 continue
      print*,'*** Erro na leitura do arquivo de dados !'
      stop             
c ......................................................................      
      end 
c **********************************************************************
c
c **********************************************************************
      subroutine parameters(nnode,numel,nen,ndf,ndfF,ndfE,ndm,numat,nin)
c **********************************************************************
c *                                                                    *
c *   Parameters                                                       *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      character*12 string
      integer nnode,numel,nen,ndf,ndfF,ndfE,ndm,numat
      integer n,j
      integer nin
c ......................................................................            
      n   = 0
      call readmacro(nin,.true.)
      write(string,'(6a)') (word(j),j=1,6)
      do while (strl .ne. 0)
         if     (string .eq. 'nnode') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)            
            read(string,*,err = 100,end = 100) nnode
            n = n + 1
         elseif (string .eq. 'numel') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)
            read(string,*,err = 100,end = 100) numel
            n = n + 1
         elseif (string .eq. 'maxno') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)           
            read(string,*,err = 100,end = 100) nen
            n = n + 1
         elseif (string .eq. 'ndf') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)            
            read(string,*,err = 100,end = 100) ndf
            n = n + 1
         elseif (string .eq. 'ndfF') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)            
            read(string,*,err = 100,end = 100) ndfF
            n = n + 1
         elseif (string .eq. 'ndfE') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)            
            read(string,*,err = 100,end = 100) ndfE
            n = n + 1
         elseif (string .eq. 'dim') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)         
            read(string,*,err = 100,end = 100) ndm
            n = n + 1
         elseif (string .eq. 'numat') then
            call readmacro(nin,.false.)
            write(string,'(12a)') (word(j),j=1,12)         
            read(string,*,err = 100,end = 100) numat
            n = n + 1
         endif
         call readmacro(nin,.false.)
         write(string,'(6a)') (word(j),j=1,6)                 
      end do
      if(n .lt. 8) goto 100
      return
c ......................................................................                        
  100 continue
      print*,'*** Erro na leitura das variaveis de controle !'
      stop       
c ......................................................................                  
      end
c **********************************************************************
c
c **********************************************************************
c *                                                                    *
c * PEGDE: Leitura das aresta com restricao                            *
c *                                                                    *
c **********************************************************************
      subroutine pegde(pedge,numel,nen,nin)
      implicit none
      include 'string.fi'
      character*30 string
      integer pedge(nen,*),numel,nen
      integer i,j,k,nin
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(j),j=1,12)
c .....................................................         
      do 110 while(string .ne. 'end')
         read(string,*,err = 200,end = 200) j
         if(j .lt. 1 .or. j .gt. numel) goto 200
c .....................................................            
         call readmacro(nin,.false.)
         write(string,'(30a)') (word(k),k=1,30)
         i = 0
c
c ... linha de comando original:
c         do 100 while(string .ne. ' ')
c
c ... formato necessário ao funcionamento em linux:
         do 100 while ( (string .ne.   ' '    ) .and.
     .                  (string .ne. CHAR(13) ) )
            i = i + 1
            if(i .gt. nen) goto 200
            read(string,*,err = 200,end = 200) pedge(i,j)
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(k),k=1,30)               
  100    continue
         call readmacro(nin,.true.)
         write(string,'(12a)') (word(j),j=1,12)
  110 continue
      return
c ......................................................................                  
  200 continue
      print*,'*** Erro na leitura das aresta com restricoes !'
      print*,'Elemento ',j,'aresta ',i
      stop       
c ......................................................................                  
      end
c **********************************************************************
c
c **********************************************************************
c *                                                                    *
c * SEGDE: Leitura das restricoes por aresta                           *
c *                                                                    *
c **********************************************************************
      subroutine segde(sedge,numel,nen,nin)
      implicit none
      include 'string.fi'
      character*30 string
      real*8  sedge(nen,*)
      integer numel,nen
      integer i,j,k,nin
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(j),j=1,12)
c .....................................................         
      do 110 while(string .ne. 'end')
         read(string,*,err = 200,end = 200) j
         if(j .lt. 1 .or. j .gt. numel) goto 200
c .....................................................            
         call readmacro(nin,.false.)
         write(string,'(30a)') (word(k),k=1,30)
         i = 0
c
c ... linha de comando original:
c         do 100 while(string .ne. ' ')
c
c ... formato necessário ao funcionamento em linux:
         do 100 while ( (string .ne.   ' '    ) .and.
     .                  (string .ne. CHAR(13) ) )
            i = i + 1
            if(i .gt. nen) goto 200
            read(string,*,err = 200,end = 200) sedge(i,j)
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(k),k=1,30)               
  100    continue
         call readmacro(nin,.true.)
         write(string,'(12a)') (word(j),j=1,12)
  110 continue
      return
c ......................................................................                  
  200 continue
      print*,'*** Erro na leitura das aresta com restricoes !'
      print*,'Elemento ',j,'aresta ',i
      stop       
c ......................................................................                  
      end
c **********************************************************************
c
c **********************************************************************
      subroutine bound(id,nnode,ndf,nin,code)
c **********************************************************************
c *                                                                    *
c *   BOUND: Leitura de restricoes                                     *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'
      integer id(ndf,*),nnode,ndf,nin,code,i,j,k
      character*30 string
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(j),j=1,12)
c .....................................................         
      do 110 while(string .ne. 'end')
         read(string,*,err = 200,end = 200) j
         if(j .lt. 1 .or. j .gt. nnode) goto 200
c .....................................................            
         call readmacro(nin,.false.)
         write(string,'(30a)') (word(k),k=1,30)
         i = 0
c
c ... linha de comando original:
c         do 100 while(string .ne. ' ')
c
c ... formato necessário ao funcionamento em linux:
         do 100 while ( (string .ne.   ' '    ) .and.
     .                  (string .ne. CHAR(13) ) )
            i = i + 1
            if(i .gt. ndf) goto 200
            read(string,*,err = 200,end = 200) id(i,j)
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(k),k=1,30)               
  100    continue
         call readmacro(nin,.true.)
         write(string,'(12a)') (word(j),j=1,12)
  110 continue
      return
c ......................................................................                        
  200 continue
      if    (code .eq. 1) then
         print*,'*** Erro na leitura das restricoes !'
      elseif(code .eq. 2) then
         print*,'*** Erro na leitura das cargas nodais !'
      elseif(code .eq. 3) then
         print*,'*** Erro na leitura das cargas nos elementos !'
      endif
      stop             
      end            
c **********************************************************************
      subroutine forces(f,nnode,ndf,nin)
c **********************************************************************
c *                                                                    *
c *   Forcas nodais e valores prescritos.                              *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'string.fi'      
      integer nnode,ndf,nin,i,j,k
      real*8 f(ndf,*)
      character*30 string            
c ......................................................................
      call readmacro(nin,.true.)
      write(string,'(12a)') (word(i),i=1,12)
      do 100 while(string .ne. 'end')
         read(string,*,err = 200,end = 200) k
         if(k .lt. 1 .or. k .gt. nnode) goto 200      
         do 10 j = 1, ndf
            call readmacro(nin,.false.)
            write(string,'(30a)') (word(i),i=1,30)
            read(string,*,err = 200,end = 200) f(j,k)
   10    continue
         call readmacro(nin,.true.)
         write(string,'(12a)') (word(i),i=1,12)
  100 continue
      return
c ......................................................................      
  200 continue
      print*,'*** Erro na leitura das forcas nodais !'
      stop             
      end
c *********************************************************************
c * readConfig : leitura das configuracoes basicas de excucao         *
c * ------------------------------------------------------------------*
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * maxmem        - memoria do vetor de trabalho                      *
c * openMpCell    - openMp nas celulas                                *
c * openMpSolver  - openMp nos Solver It                              *
c * nThreadsCell  - numero de threads nas celulas                     *
c * nThreadsSolver- numero de threads nos solver                      *
c * reord         - reordanaco do sistema de equacao                  *
c * bvtk          - saida binario para o vtk legacy                   *
c * nin           - arquivo de entrada                                *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine readConfig(maxmem,openMpCell,openMpSolver
     .                     ,nThreadsCell,nThreadsSolver,reord,bvtk,nin)
      implicit none
      include 'string.fi'
      character*15 string,macro(7)
      integer maxmem,nThreadsCell,nThreadsSolver
      logical openMpCell,openMpSolver,r(7),reord,bvtk
      integer n,j,nmacro
      integer nin
      data nmacro /7/
      data macro/'memoria        ','openMpCell     ','openMpSolver   ',
     .           'nThreadsCell   ','nThreadsSolver ','reord          ',
     .           'bvtk           '/
c .....................................................................
      n      = 0
      r(1:nmacro) = .false.
      call readmacro(nin,.true.)
      write(string,'(15a)') (word(j),j=1,12)
      do while (strl .ne. 0)
         if     (string .eq. macro(1)) then
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) maxmem
c ... convertendo de Mbytes para para numero de inteiros e 4 bytes
            maxmem = maxmem*1024*1024/4
            n = n + 1
            r(1) = .true.
         elseif (string .eq. macro(2)) then
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)
            read(string,*,err = 100,end = 100) openMpCell
            n = n + 1
            r(2) = .true.
         elseif (string .eq. macro(3)) then
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)           
            read(string,*,err = 100,end = 100) openMpSolver
            n = n + 1
            r(3) = .true.
         elseif (string .eq. macro(4)) then
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) nThreadsCell
            n = n + 1
            r(4) = .true.
         elseif (string .eq. macro(5)) then
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) nThreadsSolver
            n = n + 1
            r(5) = .true.
         elseif (string .eq. macro(6)) then
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) reord          
            n = n + 1
            r(6) = .true.
         elseif (string .eq. macro(7)) then
            call readmacro(nin,.false.)
            write(string,'(15a)') (word(j),j=1,15)            
            read(string,*,err = 100,end = 100) bvtk          
            n = n + 1
            r(7) = .true.
         endif 
         call readmacro(nin,.false.)
         write(string,'(15a)') (word(j),j=1,15)                 
      end do
c ...
      call readmacro(nin,.true.)
      call readmacro(nin,.false.)
c ......................................................................
c
c ...
      if(n .lt. 7) goto 100
      return
c ......................................................................                        
  100 continue
      print*,'*** Erro na leitura das variaveis de controle !'
      do j = 1, nmacro
        if(.NOT.r(j)) print*,macro(j), ' faltando!!'
      enddo
      stop       
c ......................................................................                  
      end
c **********************************************************************
      subroutine PressureHy(p,x,g,ix,ie,e,numel,ndm,nen)
      implicit none
      include 'idealGas.fi'
      real*8 p(*),x(ndm,*),e(10,*),g(3),xc(3),ro,temp,m,tKelvin
      real*8 fIdealGas
      integer ix(nen+1,*),ie(*),numel,ndm,nen,ma
      integer nel,j,k,no
c ...
      loopCell: do nel = 1, numel
          xc(1:3) = 0.0d0
c ... calculo do centroide            
          do j = 1, nen
            no = ix(j,nel)
            do k = 1, ndm 
              xc(k) = xc(k) + x(k,no)
            enddo
          enddo
c ... Triangulo          
          if( nen .eq. 3) then
            do j = 1, ndm 
              xc(j) = xc(j)/3.0d0
            enddo
          endif
c ... Quadrilatero          
          if( nen .eq. 4) then
            do j = 1, ndm 
              xc(j) = xc(j)*0.25d0
            enddo
          endif
c ...
         ma       = ix(nen+1,nel)
         temp     =  e(5,ma)
         tKelvin  =  tConv + temp
         m  = dO2*mMoleculaO2 + dN2*mMoleculaN2 + dAr*mMoleculaAr  
         ro       = fIdealGas(tKelvin,Pa,m,R)
         p(nel)   = ro*(g(1)*xc(1)+g(2)*xc(2)+g(3)*xc(3)) 
c .....................................................................      
      enddo loopCell
c ...
      return
      end 
c **********************************************************************
