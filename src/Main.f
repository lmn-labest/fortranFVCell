      program MVF_SIMPLE
      Use Malloc
      implicit none
      include 'string.fi'
      include 'time.fi'
      include 'openmp.fi'
      include 'simple.fi'
c ----------------------------------------------------------------------
c
c ... Variaveis da estrutura interna de macro-comandos:
c
      character*8 mc,macro(60),lmacro1(50),lmacro2(50)
      character*30 string,nCell,nNod
      character*80 str
      character*80 name
      integer  iloop1,iloop2,imacro2,imacro1,nmacro1,nmacro2,nmc
      integer  loop1,loop2,iloopAux
      integer  j,i,code
c ----------------------------------------------------------------------
      character*80 prename,filein,fileout
c ... ponteiros      
      integer*8 i_x,i_ix,i_w
      integer*8 i_nelcon,i_nn,i_ls,i_lls,i_lw,i_du,i_b,i_mdf
      integer*8 i_un,i_ie,i_e,i_md,i_w0,i_wP,i_temp
      integer*8 i_num,i_ro,i_sx
c ... transporte
      integer*8 i_t10,i_t1,i_bt10,i_bt1,i_gradT1,i_rCellt1,i_fluxlT1
      integer*8 i_adT1,i_alT1,i_auT1,i_iaT1,i_jaT1,i_idT1
      integer*8 i_pedgeT1,i_sedgeT1
      integer*8 i_pnode,i_snode,i_fnode
c ... simple
      integer*8 i_u1,i_bU1,i_alU1,i_auU1,i_adU1,i_iaU1,i_jaU1,i_idU1
      integer*8 i_u2,i_bU2,i_alU2,i_auU2,i_adU2,i_iaU2,i_jaU2,i_idU2
      integer*8 i_p,i_pC,i_Pc1,i_bPc,i_alPc,i_auPc,i_adPc
      integer*8 i_iaPc,i_jaPc,i_idPc
      integer*8 i_en,i_en0,i_bE,i_alE,i_auE,i_adE,i_iaE,i_jaE,i_idE
      integer*8 i_gradU1,i_gradU2,i_gradPc,i_gradP,i_gradE,i_gradT
      integer*8 i_gradTu
      integer*8 i_sedgeF,i_pedgeF,i_sedgeE,i_pedgeE
      integer*8 i_ddU,i_iM,i_liM,i_bU10,i_bU20,i_bE0
      integer*8 i_fluxlU1,i_fluxlU2,i_fluxlPc,i_fluxlE
      integer*8 i_rCellU1,i_rCellU2,i_rCellPc,i_rCellE
      integer*8 i_div,i_mParameter
c ......................................................................
c
c ... arquivo de impressao nos nos ( temp,flux,disp,stress,...)
      character*80 pcellname,pnodename
      integer*8 i_cell,i_node,i_nfile
      integer num_pcell,num_pnode
      integer nfiles
      parameter ( nfiles = 4)
      logical new_file(nfiles),flag_pcd
c ......................................................................
c
c ... nao-linear
      real*8 ddum
      logical itPlot
c
c ... Variaveis descritivas do problema:
c
      integer nnode,numel,ndm,ndfT1,ndfF,ndfE,nen,nshared,numat
c ... sistemas de equacao

      integer maxItSol,aux
      real*8 solvtolPcg,solvtolBcg
c ... simple
      integer neqU1,neqU2,neqPc,nadU1,nadU2,nadPc,neqE,nadE
      integer bandaU1,bandaU2,bandaPc,bandaE 
      integer matrizU1,matrizU2,matrizPc,matrizE 
      integer solverU1,solverU2,solverPc,solverE  
      logical sSimple,sEnergy
      logical unsymPc
      logical iTResSimplePlot  
c ... transporte
      real*8  solvT1
      integer neqT1,nadT1
      integer bandaT1,maxItT1
      integer matrizT1,solverT1    
      logical sTrans
      logical unsymT1
c ...
      real*8 dt,t,alpha
      real*8 cfl,re,vol,prandtl,grashof 
      integer istep
      logical bs
c ... Variaveis da interface de linha de comando
      integer nargs
      character arg*80
c ----------------------------------------------------------------------      
c
c ...  
      logical freord
c ......................................................................
c
c ... Arquivos de entrada e saida:
c
      integer nin,nout,naux,ntime,noutSimple,noutSave
      logical log_nl
      logical bvtk
      real*8  setelev 
      data nin /1/, nout /2/, naux /10/, ntime /14/, noutSimple /15/
      data noutSave /16/
      data flag_pcd /.false./
      data setelev /1.0d0/
c     arquivo de impressao de nos associados aos seguintes inteiros
c     1 - entrada padrao                                            
c     2 - saida padrao                                            
c     3 - log de varaveis nos elementos                           
c     4 - arquivo de entrada auxiliar  
c     10- arquivo de entrada auxiliar 
c     13- log do solver
c     14- arquivo de log de excucao
c     15- resido da iteracao do simple
c     16- arquivo de save para o simple/simpleC
c     nfile = 50,51,52,...,num_pcell
c     100 - arquivo auxiliar de execucao
c ... Macro-comandos disponiveis:
      data nmc /60/
      data macro/'loop1   ','mesh    ','loop2   ','solvT   ','pTrans1 ',
     .           'pgrad   ','pveloc  ','pelev   ','mshape  ','simple  ',
     .           'ppres   ','pgradP  ','penergy ','vMassEsp','pmassEsp',
     .           'ptemp   ','pdiv    ','tSimpU1 ','tSimpU2 ','tSimpPc ',
     .           'tSimpEn ','tDinamic','maxItT1 ','pResMass','save    ',
     .           'load    ','itSimple','tTransT1','setpnode','pntemp  ',
     .           'underPc ','setpcell','pctemp  ','pcgrad  ','simpleC ',
     .           'underP  ','underU  ','setelev ','dt      ','config  ',
     .           'pgradU1 ','pgradU2 ','pgradU3 ','skewC   ','pgeo    ',
     .           'pgradT  ','pgradE  ','        ','        ','        ',
     .           'pvort2D ','        ','        ','        ','        ',
     .           'underRo ','tolPcg  ','tolBiPcg','maxItSol','stop    '/
c ----------------------------------------------------------------------
c
c ...
      bvtk = .true. 
c ...      
      iloop1 = 0
      iloop2 = 0
      istep  = 0
      dt     = 1.0d0
      alpha  = 1.0d0
c ...
      cfl     = 0.0d0
      prandtl = 0.0d0
      re      = 0.0d0
      vol     = 0.0d0
      grashof = 0.0d0 
c ... openmp
      openmpCell         = .false.
      openmpSolver       = .false.
      nThreadsCell       = 4
      nThreadsSolver     = 2
c ... Euller Backward de segunda ordem
      bs     = .true.
c ... reordenacao da equacoes para a diminucao da banda
      freord = .true.
c ... tipo de armazenamento ( CSRD- 1 -> ad;a | CSR -> 2 - a | CSRC -> 3 - ad;au;al)
c ... CSRD - CSR com diagonal principal separada 
c     CSR  - CSR 
c     CSRC - CSRC
      matrizT1= 1
      matrizU1= 1
      matrizU2= 1
      matrizPc= 1
      matrizE = 1 
c ... tsolver               ( PCG - 1|PBiCGSTAB - 2|PGMRES -3 )
      solverT1 = 2
      solverU1 = 2
      solverU2 = 2
      solverPc = 1
      solverE  = 2
      maxItSol = 300000
      solvtolPcg = 1.0d-14
      solvtolBcg = 1.0d-14
c ... nao linear
      log_nl   = .false.
      itPlot   = .false.
c ... transporte
      sTrans   = .false.
      unsymT1  = .true.
      maxItT1  =  100
      solvT1   =  1.0d-6
c ... simple
      skewnessCorrection = .false.
      sSimple            = .false.
      sEnergy            = .false.
      iTResSimplePlot    = .true.
      underU             = 1.0d0
      underP             = 1.0d0
      underPc            = 1.0d0
      underRo            = 1.0d0
      cPc0               = 1.0d0
      cPc                = 1.0d0  
      simpleC            =.false.
      unsymPc            =.false.
      vMass              =.false.
      tDinamico          =.false.
      closed             =.false.
      singularPressure   =.false.
      maxItSimple        = 200
      solvSimplePc       = 1.0d-6
      solvSimpleU1       = 1.0d-6
      solvSimpleU2       = 1.0d-6
      solvEnergy         = 1.0d-6
      massa              = 1.0d0
      massa0             = 1.0d0
      fluxoM             = 0.0d0
c ......................................................................
c
    5 continue 
c
c ... Abertura de arquivos:    
      nargs = iargc()
   10 continue
      if(nargs .gt. 0) then
        call getarg(1,arg)
        filein = arg
      else
        print*, 'Arquivo de dados:'
        read(*,'(a)') filein
      endif
      open(nin, file= filein, status= 'old', err=15, action= 'read')
      goto 20
   15 continue
      print*, 'Arquivo nao existente !'
      nargs = 0
      goto 10
c ......................................................................
c
c .... arquivo de saida     
   20 continue
      if(nargs .eq. 2) then
        call getarg(2,arg)
        prename = arg
      else
        print*, 'Arquivo de saida: '
        read(*,'(a)') prename
      endif
c ... log de contrelo de malha nao ortognais      
c      fileout = name(prename,0,6) 
c      open(3, file=fileout ,action= 'write')
c      write(3,'(a)')'#Controle do gradiente primario e secondudario '
c      write(3,'(a)')
c     .'#(nel) (difusao direta) (difusao secundaria) (n*e) (cv) (cvc)
c     . (pe)'
c ......................................................................
c
c ... log de controle do nao-linear
      if(log_nl) then      
        fileout = name(prename,0,10) 
        open(10, file=fileout ,action= 'write')
        write(10,'(a)')'#Controle do solver não-linear '
        write(10,'(a)')
     .  '#(it) (R(i)) (R(i)/R(i-1))'
       endif
c ......................................................................
c
c ... log de controle do solver      
      fileout = name(prename,0,50) 
      open(13, file=fileout ,action= 'write')
      write(13,'(a)')'#Controle do solver iterativo '
c ......................................................................
c
c ... log de controle do simple      
      fileout = name(prename,0,51) 
      open(15, file=fileout ,action= 'write')
      write(15,'(a)')'#Controle do simple bPc Ru1 Ru2'
c ......................................................................
c
c ... controle de tempo
      solvTime        = 0.d0
      elmTime         = 0.d0
      vecTime         = 0.d0
      dsTime          = 0.d0
      grTime          = 0.d0
      numeqTime       = 0.d0
      reordTime       = 0.d0  
      totalTime       = 0.d0
      matvecTime      = 0.d0
      dotTime         = 0.d0
      readTime        = 0.d0
      vizTime         = 0.d0
      preconTime      = 0.d0
      solvU1Time      = 0.d0
      solvU2Time      = 0.d0
      solvPctime      = 0.d0
      solvEtime       = 0.d0
      grSimpleTime    = 0.d0
      elmU1Time       = 0.d0
      elmU2Time       = 0.d0
      elmPcTime       = 0.d0
      elmETime        = 0.d0
      posVelocityTime = 0.d0
      simpleUpdateTime= 0.d0
      simpleTime      = 0.d0
      totaltime       = get_time() 
c ......................................................................
c
c ... Leitura dos macro-comandos:
c
c ......................................................................
   50 continue
      if (iloop1 .gt. 0) goto 55
      if (iloop2 .gt. 0) goto 60
      call readmacro(nin,.true.)
      write(mc,'(8a)') (word(i),i=1,8)
      do j = 1, nmc
        if (mc .eq. macro(j)) goto 70
      enddo 
      goto 50
c ... loop       
   55 continue
      if (iloop1 .eq. 0) then
         goto 60
      else
         if (imacro1 .eq. 0 .or. imacro1 .eq. nmacro1) then
            imacro1 = 1
         else
            imacro1 = imacro1 + 1
         endif
         mc     = lmacro1(imacro1)
         iloop1 = iloop1 - 1
      endif
      do j = 1, nmc
         if (mc .eq. macro(j)) goto 70
      enddo 
      goto 55
c .. loop2       
   60 continue
      if (nmacro2 .eq. imacro2) then
        iloop2  = iloop2 - 1
        if(iloop2 .gt. 0) iloop1  = iloopAux
        imacro2 = 0
        goto 50
      else
        imacro2 = imacro2 + 1
        mc      = lmacro2(imacro2)
      endif
      do j = 1, nmc
         if (mc .eq. macro(j)) goto 70
      enddo 
      goto 50
c ....................................................................
   70 continue
      goto(100,200,300,400,500,!loop1  ,mesh    ,loop2   ,solvT   ,pTrans1 
     . 600, 700, 800, 900,1000, !pgrad  ,pveloc  ,pelev   ,mshape  ,simple  
     .1100,1200,1300,1400,1500, !ppres  ,pgradP  ,penergy ,vMassEsp,pmassEsp
     .1600,1700,1800,1900,2000, !ptemp  ,pdiv    ,tSimpU1 ,tSimpU2 ,tSimpPc 
     .2100,2200,2300,2400,2500, !tSimpEn,tDinamic,maxItT1 ,pResMass,save    
     .2600,2700,2800,2900,3000, !load   ,itSimple,tTransT1,setpnode,pntemp  
     .3100,3200,3300,3400,3500, !underPc,setpcell,pctemp  ,pcgrad  ,simpleC 
     .3600,3700,3800,3900,4000, !underP ,underU  ,setelev ,dt      ,config  
     .4100,4200,4300,4400,4500, !pgradU1,pgradU2 ,pgradU3 ,skewC   ,pgeo    
     .4600,4700,4800,4900,5000, !pgradT ,pgradE  ,        ,        ,        
     .5100,5200,5300,5400,5500, !       ,        ,        ,        ,        
     .5600,5700,5800,5900,6000)j!underRo,tolPcg  ,tolBiPcg,maxItSol,stop    
c ----------------------------------------------------------------------
c
c ... Execucao dos macro-comandos:
c
c ----------------------------------------------------------------------
c
c ... Macro-comando LOOP:
c ......................................................................
  100 continue
      call readmacro(nin,.false.)
      write(string,'(12a)') (word(i),i=1,12)
      read(string,*,err = 120,end = 120) loop1
      nmacro1 = 0
      imacro1 = 0
      iloop1  = 0
  110 continue
      call readmacro(nin,.true.)
      write(mc,'(8a)') (word(i),i=1,8)        
      if (mc .eq. 'next1') then 
        if(loop2 .eq. 0 ) then 
          goto 50 
        else
          iloopAux = loop1*nmacro1 
          goto 310
        endif
      endif
      nmacro1= nmacro1 + 1
      iloop1 = loop1*nmacro1
      lmacro1(nmacro1) = mc
      goto 110
  120 continue
      print*,'Erro na leitura da macro (LOOP) !'
      goto 5000            
c ----------------------------------------------------------------------
c
c ... Macro-comando MESH:
c
c ......................................................................
  200 continue
      print*, 'Macro MESH' 
c
c ... Inicializacao da estrutura de gerenciamento de memoria:
c
      call init_malloc()
c ......................................................................
c
c ... iniciando openmp
c$    num_threads = omp_get_max_threads()
      write(*,'(1x,a,1x,i3)') 'numero threads maximo: ',num_threads
      if(openMpCell)  write(*,'(1x,a,1x,i3)') 'numero threads Cell  : '
     .              ,nThreadsCell
      if(openMpSolver)write(*,'(1x,a,1x,i3)') 'numero threads Solver: ' 
     .               ,nThreadsSolver
c ..................................................................... 
c
c ...
      if(freord) write(*,'(1x,a,1x,i3)') 'Reverse cuthill mckee.'
c .....................................................................
c
c ... leitura do arquivo de dados
      readtime = get_time()  
      call rdat_mvf(i_x,i_ix,i_pedgeT1,i_sedgeT1,i_pedgeF,i_sedgeF
     .             ,i_pedgeE,i_sedgeE,i_pnode,i_snode
     .             ,i_fnode,i_e,i_ie,i_w,i_t1,i_u1,i_u2,i_p,i_temp
     .             ,i_num,nnode,numel
     .             ,ndm,ndfF,ndfE,ndfT1,sTrans,sSimple,sEnergy,nen
     .             ,numat,nin)
      readtime = get_time() - readtime  
c       call si(ia(i_u),numel,ndf)
c       stop 
c       call wi(ia(i_w),numel)
c ......................................................................
c
c ... identificacao dos vizinhos
      viztime = get_time()
      call viz(ia(i_ix),i_nelcon,nnode,numel,nen) 
      viztime = get_time() - viztime
c .....................................................................
c
c ... transforma condicoes nodais e condicoes por aresta elemento
c      if(sTrans) then
c        call cbound(ia(i_snode),ia(i_fnode),ia(i_sedgeT1),ia(i_ix)
c     .             ,ia(i_pnode),ia(i_pedgeT1),numel,nen,nen,ndfT1)
c      endif 
c .....................................................................
c
c ... diminuicao da banda da matriz
      reordtime = get_time()
      call reord(ia(i_nelcon),ia(i_num),numel,nen,freord)
      reordtime = get_time() - reordtime
c.....................................................................
c 
c ....................................................................
c ... numeracao das equacoes
c
c ... equacao de transporte
      neqT1 = 0
      if(sTrans) then
        i_idT1 = alloc_4('idT1    ',nen,numel)
        numeqTime = get_time()
        call numeq(ia(i_idT1),ia(i_num),ia(i_nelcon),numel,nen,neqT1)
        numeqTime = get_time() - numeqTime
      endif
c .....................................................................
c
c ... equacao u,v,p
      neqPc = 0
      neqU1 = 0
      neqU2 = 0
      if(sSimple) then
        i_idPc = alloc_4('idPc    ',nen,numel)
        i_idU1 = alloc_4('idU1    ',nen,numel)
        i_idU2 = alloc_4('idU2    ',nen,numel)
        numeqTime = get_time() - numeqTime
        call numeq(ia(i_idPc),ia(i_num),ia(i_nelcon),numel,nen,neqPc)
        call numeq(ia(i_idU1),ia(i_num),ia(i_nelcon),numel,nen,neqU1)
        call numeq(ia(i_idU2),ia(i_num),ia(i_nelcon),numel,nen,neqU2)
        numeqTime = get_time() - numeqTime
      endif
c .....................................................................
c
c ... equacao Energia
      i_iaE = 1
      neqE  = 0
      if(sEnergy) then
        i_idE  = alloc_4('idE     ',nen,numel)
        numeqTime = get_time() - numeqTime
        call numeq(ia(i_idE ),ia(i_num),ia(i_nelcon),numel,nen,neqE )
        numeqTime = get_time() - numeqTime
      endif
c .....................................................................
c
c ... estrutura de dados para matriz esparça
c
c ... equacao de transporte
      if(sTrans) then
        dsTime = get_time() - dsTime
        call datastruc(i_idT1,i_iaT1,i_jaT1
     .                ,i_adT1,i_auT1,i_alT1,neqT1,nadT1,nen
     .                ,matrizT1,bandaT1,unsymT1
     .                ,'idT1    ','iaT1    ','jaT1    ','adT1    '
     .                ,'auT1    ','alT1    ')
        dsTime = get_time() - dsTime
      endif
c .....................................................................
c
c ... equacao E
      nadE  = 0
      i_iaE = 1
      i_jaE = 1
      i_adE = 1
      i_auE = 1
      i_alE = 1 
      if(sEnergy) then
        dsTime = get_time() - dsTime
        call datastruc(i_idE,i_iaE,i_jaE,i_adE,i_auE,i_alE,neqE
     .                ,nadE,nen,matrizE,bandaE,.true.
     .                ,'idE     ','iaE     ','jaE     ','adE     '
     .                ,'auE     ','alE     ')
        dsTime = get_time() - dsTime
      endif
c .....................................................................
c
c ... equacao u,v,p
      if(sSimple) then
        dsTime = get_time() - dsTime
        call datastruc(i_idU1,i_iaU1,i_jaU1,i_adU1,i_auU1,i_alU1,neqU1
     .                ,nadU1,nen,matrizU1,bandaU1,.true.
     .                ,'idU1    ','iaU1    ','jaU1    ','adU1    '
     .                ,'auU1    ','alU1    ')
        call datastruc(i_idU2,i_iaU2,i_jaU2,i_adU2,i_auU2,i_alU2,neqU2
     .                ,nadU2,nen,matrizU2,bandaU2,.true.
     .                ,'idU2    ','iaU2    ','jaU2    ','adU2    '
     .                ,'auU2    ','alU2    ')
        call datastruc(i_idPc,i_iaPc,i_jaPc,i_adPc,i_auPc,i_alPc,neqPc
     .                ,nadPc,nen,matrizPc,bandaPC,unsymPc
     .                ,'idPc    ','iaPc    ','jaPc    ','adPc    '
     .                ,'auPc    ','alPc    ')
        dsTime = get_time() - dsTime
      endif
c .....................................................................
c
c ... vetor globais
c
      i_ls        = alloc_8('ls      ',ndm*nen,numel)
      i_mdf       = alloc_8('mdf     ',  1,nnode)
      i_md        = alloc_4('md      ',  1,nnode)
      i_un        = alloc_8('un      ',  1,nnode)
      i_rO        = alloc_8('Ro      ',  3,numel)
      aux         = max(neqT1,neqPc,neqU1,neqU2,neqE)
      i_sx        = alloc_8('sx      ',  1,aux)
      call azero(ia(i_ls),  ndm*nen*numel)
      call azero(ia(i_mdf),nnode)
      call mzero(ia(i_md),nnode)
      call azero(ia(i_un),nnode)
      call azero(ia(i_rO),3*numel)
c      call azero(ia(i_sx),naux)
c ... transporte
      if(sTrans) then 
        i_du        = alloc_8('du      ',  1,neqT1)
        i_t10       = alloc_8('t10     ',  1,neqT1)
        i_gradT1    = alloc_8('gradT1  ',ndm,numel)
        i_fluxlT1   = alloc_8('fluxlT1 ',  1,numel)
        i_bT1       = alloc_8('bT1     ',  1,neqT1)
        i_bT10      = alloc_8('bT10    ',  1,neqT1)
        i_rCellT1   = alloc_8('rCellT1 ',  1,neqT1)
        call azero(ia(i_du),neqT1)
        call azero(ia(i_t10),neqT1)
        call azero(ia(i_gradT1),ndm*numel)
        call azero(ia(i_fluxlT1),numel)
        call azero(ia(i_bT1),neqT1)
        call azero(ia(i_bT10),neqT1)
        call azero(ia(i_rCellT1),neqT1)
      endif
c .....................................................................
c
c ... simple
      if(sSimple) then 
        i_gradU1    = alloc_8('gradU1  ',ndm,numel)
        i_gradU2    = alloc_8('gradU2  ',ndm,numel)
        i_gradPc    = alloc_8('gradPc  ',ndm,numel)
        i_gradTu    = alloc_8('TensorU ',ndm*ndm,numel)
        i_gradP     = alloc_8('gradP   ',ndm,numel)
        i_div       = alloc_8('div     ',  1,numel)
        i_pC        = alloc_8('presC   ',  1,neqPc)
        i_pC1       = alloc_8('presC1  ',  1,neqPc)
        i_bU1       = alloc_8('bU1     ',  1,neqU1)
        i_bU2       = alloc_8('bU2     ',  1,neqU2)
        i_bU10      = alloc_8('bU10    ',  1,neqU1)
        i_bU20      = alloc_8('bU20    ',  1,neqU2)
        i_bPc       = alloc_8('bPc     ',  1,neqPc)
        i_fluxlU1   = alloc_8('fluxlU1 ',  1,numel)
        i_fluxlU2   = alloc_8('fluxlU2 ',  1,numel)
        i_fluxlPc   = alloc_8('fluxlPc ',  1,numel)
        i_rCellU1   = alloc_8('rCellU1 ',  1,numel)
        i_rCellU2   = alloc_8('rCellU2 ',  1,numel)
        i_rCellPc   = alloc_8('rCellPc ',  1,numel)
        i_ddU       = alloc_8('ddU     ',ndfF-1,numel)
        i_iM        = alloc_8('iM      ',2*(ndfF-1),numel)
        i_wP        = alloc_8('vpseudo ',ndfF-1,numel)
        i_w0        = alloc_8('v0      ',ndfF-1,numel)
        i_mParameter= alloc_8('mPar    ',9,numel)
        call azero(ia(i_gradU1),ndm*numel)
        call azero(ia(i_gradU2),ndm*numel)
        call azero(ia(i_gradPc),ndm*numel)
        call azero(ia(i_gradP) ,ndm*numel)
        call azero(ia(i_gradTu) ,ndm*ndm*numel)
        call azero(ia(i_div  ) ,    numel)
        call azero(ia(i_bU1),neqU1)
        call azero(ia(i_bU2),neqU2)
        call azero(ia(i_bU10),neqU1)
        call azero(ia(i_bU20),neqU2)
        call azero(ia(i_bPC),neqPc)
        call azero(ia(i_fluxlU1),numel)
        call azero(ia(i_fluxlU2),numel)
        call azero(ia(i_fluxlPC),numel)
        call azero(ia(i_rCellU1),numel)
        call azero(ia(i_rCellU2),numel)
        call azero(ia(i_rCellPC),numel)
        call azero(ia(i_ddU),numel*(ndfF-1))
        call azero(ia(i_iM) ,2*numel*(ndfF-1))
        call azero(ia(i_wP) ,numel*(ndfF-1))
        call azero(ia(i_w0) ,numel*(ndfF-1))
        call azero(ia(i_mParameter),6*numel)
       endif
c ......................................................................
c
c ... Energia
      i_gradE       = 1
      i_bE          = 1
      i_bE0         = 1
      i_en0         = 1
      i_en          = 1
      i_fluxlE      = 1
      i_rCellE      = 1
      if(sEnergy) then
        i_gradE     = alloc_8('gradE   ',ndm,numel)
        i_gradT     = alloc_8('gradT   ',ndm,numel)
        i_bE        = alloc_8('bE      ',  1,neqE )
        i_bE0       = alloc_8('bE0     ',  1,neqE )
        i_en        = alloc_8('en      ',  1,neqE )
        i_en0       = alloc_8('en0     ',  1,neqE )
        i_fluxlE    = alloc_8('fluxlE  ',  1,numel)
        i_rCellE    = alloc_8('rCellE  ',  1,numel)
        call azero(ia(i_gradE),ndm*numel)
        call azero(ia(i_gradT),ndm*numel)
        call azero(ia(i_bE ),neqE )
        call azero(ia(i_bE0 ),neqE )
        call azero(ia(i_en0 ),neqE )
        call azero(ia(i_en  ),neqE )
        call azero(ia(i_fluxlE ),numel)
        call azero(ia(i_rCellE ),numel)
      endif
c ...
      if(sSimple) then
        call checkD(ia(i_pedgeF),numel,nen,singularPressure,closed)
      endif
c ......................................................................
c
c ... 
      if(sEnergy) then
        call EnthalpyForTemp(ia(i_en),ia(i_temp),ia(i_e),ia(i_w),numel
     .                      ,ndm,.false.,1)
        if(closed) then
          call closedPresure(ia(i_temp),ia(i_x),ia(i_e),ia(i_ie)
     .                      ,ia(i_nelcon),ia(i_ix),cPc0,numel,ndm,nen
     .                      ,nen,6,3)
          cPc = cPc0
        endif   
        call aequalb(ia(i_en0),ia(i_en),neqE)
        call aequalbVetor(ia(i_w0),ia(i_w),numel,ndfF-1)
c ... ro(1,*) = PM/(RT)
        call varMassEsp(ia(i_ro),ia(i_temp),numel,1,1,1)
c ... ro(2,*) = PM/(RT)
        call varMassEsp(ia(i_ro),ia(i_temp),numel,1,1,2)
c ... ro(3,*) = PM/(RT)
        call varMassEsp(ia(i_ro),ia(i_temp),numel,1,1,3)
c ... e(2,*) = PM/(RT)
c        call tRo(ia(i_ro),ia(i_e),ia(i_ie),ia(i_ix),numel,nen,.true.)
c ... parametors iniciais
        call cellParameter(ia(i_mParameter),ia(i_ro),ia(i_w),ia(i_x)
     .                    ,ia(i_e),ia(i_sedgeF),ia(i_pedgeF)
     .                    ,ia(i_ie),ia(i_nelcon),ia(i_ix)
     .                    ,numel,ndm,nen,nen,ndfF,dt,7,1)
        call calParameter(ia(i_mParameter),ia(i_div),numel,cfl,re
     .                   ,prandtl,grashof,vol,massa0,fluxoM,dt
     .                   ,.false.)
        massa = massa0
c ... incompressivel 
       else
c ...  ro(1,*) =  ro(2,*) =  ro(3,*) = e(2,*)
         call tRo(ia(i_ro),ia(i_e),ia(i_ie),ia(i_ix),numel,nen,.false.)
       endif
c ......................................................................
c
c ...
      if(sTrans) then
         if(bs)  call aequalb(ia(i_t10),ia(i_t1),neqT1)
      endif
c ......................................................................
c
c ... matriz para least square
c      if(sTrans) then
c        call lsform(ia(i_ls),ia(i_x),ia(i_sedge),ia(i_ie),ia(i_nelcon)
c     .             ,ia(i_pedge),ia(i_ix),ia(i_lsedge),ia(i_lx)
c     .             ,ia(i_lpedge),ia(i_lls),ia(i_lviz),numel,ndm,nen,nen)
c      endif
c ......................................................................
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando LOOP2:
c
c ......................................................................
  300 continue
      call readmacro(nin,.false.)
      write(string,'(12a)') (word(i),i=1,12)
      read(string,*,err = 320,end = 320) loop2
      nmacro2= 0
      imacro2= 0
      iloop2 = 0
  310 continue
      call readmacro(nin,.true.)
      write(mc,'(8a)') (word(i),i=1,8)
      if (mc .eq. 'loop1') goto 100        
      if (mc .eq. 'next2') goto 55
      nmacro2        = nmacro2+ 1
      iloop2         = loop2 
      lmacro2(nmacro2)= mc
      goto 310
  320 continue
      print*,'Erro na leitura da macro (LOOP2) !'
      goto 5000  
c .....................................................................
c
c ----------------------------------------------------------------------
c
c ... Macro-comando SOLVT:
c
c ......................................................................
  400 continue
      print*, 'Macro SOLVT'
c ...
      if(.NOT.sSimple) then
        istep = istep + 1
        t     = t + dt
      endif
      write(* ,'(1x,a,i8,a,f16.8)'),'STEP ',istep, '   Time (s) : ', t
      write(* ,'(1x,a,es20.5)'),'delta T  : ',dt 
c ......................................................................
c
c ...
      call transporte(ia(i_x)      ,ia(i_ix)     ,ia(i_e)
     .               ,ia(i_ie)     ,ia(i_nelcon) ,ia(i_pedgeT1)
     .               ,ia(i_sedgeT1),ia(i_w)      ,ia(i_num)
     .               ,ia(i_ls)     ,ia(i_sx)   
     .               ,ia(i_gradT1) ,ia(i_fluxlT1),ia(i_rCellt1)
     .               ,ia(i_adT1)   ,ia(i_auT1)   ,ia(i_alT1)
     .               ,ia(i_bT1)    ,ia(i_bT10)   
     .               ,ia(i_iaT1)   ,ia(i_jaT1)   ,ia(i_t1)
     .               ,ia(i_t10)    ,ia(i_un)     ,ia(i_ro)
     .               ,numel        ,nnode        ,ndm
     .               ,nen          ,nshared      ,ndfT1
     .               ,dt           ,t            ,matrizT1
     .               ,neqT1        ,nadT1        ,solverT1
     .               ,solvTolBcg   ,maxItSol     ,istep
     .               ,unsymT1      ,solvT1       ,maxItT1
     .               ,bs)
c ......................................................................       
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando PTRANS1:
c
c ......................................................................
  500 continue
      print*,'Macro PTRANS1'
      call uformnode(ia(i_un),ia(i_t1),ia(i_gradT1),ia(i_fluxlT1)
     .              ,ia(i_x),ia(i_mdf),ia(i_ix),ia(i_md),nnode,numel
     .              ,ndm,nen,1,2)
      fileout = name(prename,istep,1)
      nCell = 'elTrans1'
      nNod  = 'noTrans1'
      call write_res_vtk(ia(i_ix),ia(i_x),ia(i_t1),ia(i_un),nnode,numel
     .                  ,ndm,nen,ndfT1,fileout,nCell,nNod ,bvtk
     .                  ,4,t,istep,nout)
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando PGRAD:   
c
c ......................................................................
  600 continue
c      print*,'Macro PGRAD'
c      fileout     = name(prename,istep,5) 
c      i_nn        = alloc_8('nn      ',ndm,nnode)
c      call uformnode(ia(i_nn),ia(i_grad),ddum,ddum,ia(i_x),ia(i_mdf)
c     .              ,ia(i_ix),ia(i_md),nnode,numel,ndm,nen,ndm,2)
c      call write_res_vtk(ia(i_ix),ia(i_x),ia(i_grad),ia(i_nn),nnode
c     .                  ,numel,ndm,nen,ndf,fileout,.false.,5,t,istep
c     .                  ,nout)
c      i_nn        = dealloc('nn      ')
      goto 50
c ----------------------------------------------------------------------
c
c ----------------------------------------------------------------------
c
c ... Macro-comando PVELOC:   
c
c ......................................................................
  700 continue
      print*,'Macro PVELOC'
c ...
      fileout     = name(prename,istep,7) 
      i_nn        = alloc_8('nn      ',ndm,nnode)
c ......................................................................
c
c ...
      call guess(ia(i_u1),ia(i_u2),ia(i_w),numel,ndfF-1)
c
      call gform(ia(i_u1),ia(i_gradU1),ia(i_fluxlU1),ia(i_x)
     .          ,ia(i_sedgeF)
     .          ,ia(i_e),ia(i_ls),ia(i_ie),ia(i_nelcon),ia(i_pedgeF)
     .          ,ia(i_ix),numel,ndm,nen,nen,1,ndfF,2,1,1)
c
      call gform(ia(i_u2),ia(i_gradU2),ia(i_fluxlU2),ia(i_x)
     .          ,ia(i_sedgeF)
     .          ,ia(i_e),ia(i_ls),ia(i_ie),ia(i_nelcon),ia(i_pedgeF)
     .          ,ia(i_ix),numel,ndm,nen,nen,1,ndfF,2,2,1)
c ......................................................................
      call JacobionMatrixUform(ia(i_gradTu),ia(i_gradU1),ia(i_gradU2)
     .                        ,numel,ndm)
c
      call uformnode(ia(i_nn),ia(i_w),ia(i_gradTu),ia(i_fluxlU1),ia(i_x)
     .              ,ia(i_mdf),ia(i_ix),ia(i_md),nnode,numel,ndm,nen
     .              ,ndfF-1,3)
c ......................................................................
c
c ...      
      nCell = 'elVelocidade'
      nNod  = 'noVelocidade'
      call write_res_vtk(ia(i_ix),ia(i_x),ia(i_w),ia(i_nn),nnode
     .                  ,numel,ndm,nen,ndfF,fileout,nCell,nNod 
     .                  ,bvtk,7,t,istep
     .                  ,nout)
      i_nn        = dealloc('nn      ')
      goto 50
c ----------------------------------------------------------------------
c
c ----------------------------------------------------------------------
c
c ... Macro-comando PELEV:
c
c ......................................................................
  800 continue
      print*,'Macro PELEV'
      call uformnode(ia(i_un),ia(i_p),ddum,ddum,ia(i_x)
     .              ,ia(i_mdf),ia(i_ix),ia(i_md),nnode,numel,ndm,nen,1
     .              ,1)
      i_nn        = alloc_8('xe      ',  3,nnode)
      fileout = name(prename,istep,8)
      nCell = ''
      nNod = ''
      call elev(ia(i_un),ia(i_x),ia(i_nn),nnode,ndm,setelev,1,1)
      call write_res_vtk(ia(i_ix),ia(i_nn),ddum,ddum,nnode,numel
     .                  ,  3,nen,1,fileout,nCell, nNod ,bvtk,0,t
     .                  ,istep,nout)
      i_nn        = dealloc('xe      ')
      goto 50
c ----------------------------------------------------------------------
c
c ----------------------------------------------------------------------
c
c ... Macro-comando: MSHAPE: grafo da matriz de coeficientes 
c
c ......................................................................
  900 continue
      print*, 'Macro mshape'
      fileout = name(prename,0,24)
      open(naux, file= fileout)
      call pmshape(ia(i_iaT1),ia(i_jaT1),neqT1,naux)
      close(naux)
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: SIMPLE/SIMPLEC/SIMPLER 2D
c
c ......................................................................
 1000 continue
      if(simpleC) then 
        print*, 'Macro SIMPLEC'
      else
        print*, 'Macro SIMPLE'
      endif
c ... testa se existe arquivo de parada
      open(100,file='stop.mvf',status='old',err=1010)
      print*, 'Arquivo de parada achado!!'
      goto 6090     
c ...      
 1010 continue  
      istep = istep + 1
      t     = t + dt
      write(* ,'(1x,a,i8,a,f16.8)'),'STEP ',istep, '   Time (s) : ', t
      write(* ,'(1x,a,es20.5)'),'delta T  : ',dt 
      write(* ,'(1x,a,es20.5)'),'CFL      : ',cfl
      write(* ,'(1x,a,es20.5)'),'Reynalds : ', re
      write(* ,'(1x,a,es20.5)'),'Prandlt  : ', prandtl
      if(sEnergy) write(* ,'(1x,a,es20.5)'),'rMassa   : ', massa0/massa
      write(* ,'(1x,a,es20.5)'),'Fluxo M  : ', fluxoM
c      write(* ,'(1x,a,es20.5)'),'Grashof  : ', grashof
c .....................................................................
c
c ...
      simpleTime = get_time() - simpleTime
      call simple(ia(i_x)   ,ia(i_ix)     ,ia(i_e)        ,ia(i_ie)
     .     ,ia(i_nelcon)    ,ia(i_pedgeF) ,ia(i_sedgeF)   ,ia(i_pedgeE)
     .     ,ia(i_sedgeE)    ,ia(i_w)      ,ia(i_w0)       ,ia(i_wP)    
     .     ,ia(i_num)       ,ia(i_ls)     ,ia(i_gradU1)   ,ia(i_gradU2)
     .     ,ia(i_gradP)     ,ia(i_gradPc) ,ia(i_gradE)    ,ia(i_div)   
     .     ,ia(i_mParameter),ia(i_iM)     ,ia(i_fluxlU1)  ,ia(i_fluxlU2)
     .     ,ia(i_fluxlPc)   ,ia(i_fluxlE) ,ia(i_rCellU1)  ,ia(i_rCellU2)
     .     ,ia(i_rCellPc)   ,ia(i_rCellE) ,ia(i_adU1)     ,ia(i_auU1)   
     .     ,ia(i_alU1)      ,ia(i_bU1)    ,ia(i_bU10)     ,ia(i_iaU1)   
     .     ,ia(i_jaU1)      ,ia(i_adU2)   ,ia(i_auU2)     ,ia(i_alU2)   
     .     ,ia(i_bU2)       ,ia(i_bU20)   ,ia(i_iaU2)     ,ia(i_jaU2)   
     .     ,ia(i_adPc)      ,ia(i_auPc)   ,ia(i_alPc)     ,ia(i_bPC)    
     .     ,ia(i_iaPc)      ,ia(i_jaPc)   ,ia(i_adE)      ,ia(i_auE)    
     .     ,ia(i_alE)       ,ia(i_bE)     ,ia(i_bE0)      ,ia(i_iaE) 
     .     ,ia(i_jaE)       ,ia(i_u1 )    ,ia(i_u2)       ,ia(i_p)   
     .     ,ia(i_pC)        ,ia(i_Pc1)    ,ia(i_en)       ,ia(i_en0)    
     .     ,ia(i_ro)        ,ia(i_un)     ,ia(i_mdf)      ,ia(i_md)   
     .     ,ia(i_ddU)      ,ia(i_temp)    ,ia(i_sx)       ,numel
     .     ,nnode          ,ndm        
     .     ,nen            ,nen           ,ndfF           ,ndfE
     .     ,dt             ,t             ,matrizU1       ,matrizU2     
     .     ,matrizPc       ,matrizE       ,neqU1          ,neqU2
     .     ,neqPc          ,neqE          ,nadU1          ,nadU2  
     .     ,nadPc          ,nadE          ,solverU1       ,solverU2  
     .     ,solverPc       ,solverE       ,solvTolPcg     ,solvTolBcg   
     .     ,maxItSol       ,noutSimple    ,itResSimplePlot,sEnergy  
     .     ,istep          ,cfl           ,re             ,prandtl 
     .     ,grashof        ,vol           ,unsymPc        ,bs)
      simpleTime = get_time() - simpleTime 
      goto 50
c ----------------------------------------------------------------------
c
c ......................................................................
c
c ... Macro-comando: PPRES
c
c ......................................................................
 1100 continue
      print*, 'Macro PPRES'
      call gform(ia(i_p),ia(i_gradP),ia(i_fluxlPc),ia(i_x),ia(i_sedgeF)
     .            ,ia(i_e),ia(i_ls),ia(i_ie),ia(i_nelcon),ia(i_pedgeF)
     .            ,ia(i_ix),numel,ndm,nen,nen,1,ndfF,2,3,2)
      call uformnode(ia(i_un),ia(i_p),ia(i_gradPc),ia(i_fluxlPc),ia(i_x)
     .             ,ia(i_mdf),ia(i_ix),ia(i_md),nnode,numel,ndm,nen,1
     .             ,3)
      fileout = name(prename,istep,11) 
      nCell = 'elPressure'
      nNod  = 'noPressure'
      call write_res_vtk(ia(i_ix),ia(i_x),ia(i_p),ia(i_un),nnode,numel
     .                  ,ndm,nen,1,fileout,nCell,nNod ,bvtk,4,t
     .                  ,istep,nout)
      goto 50
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando: PGRADP
c
c ......................................................................
 1200 continue
      print*, 'Macro PGRADP'
      fileout     = name(prename,istep,12) 
      i_nn        = alloc_8('nn      ',ndm,nnode)
      call gform(ia(i_p),ia(i_gradP),ia(i_fluxlPc),ia(i_x),ia(i_sedgeF)
     .            ,ia(i_e),ia(i_ls),ia(i_ie),ia(i_nelcon),ia(i_pedgeF)
     .            ,ia(i_ix),numel,ndm,nen,nen,1,ndfF,2,3,2)
      call uformnode(ia(i_nn),ia(i_gradP),ddum,ddum,ia(i_x),ia(i_mdf)
     .              ,ia(i_ix),ia(i_md),nnode,numel,ndm,nen,ndm,2)
      nCell = 'elGradP'
      nNod  = 'noGradP'
      call write_res_vtk(ia(i_ix),ia(i_x),ia(i_gradP),ia(i_nn),nnode
     .                  ,numel,ndm,nen,1,fileout,nCell,nNod ,bvtk,5,t
     .                  ,istep,nout)
      i_nn        = dealloc('nn      ')
      goto 50
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando: PENERGY
c
c ......................................................................
 1300 continue
      print*, 'Macro PENERGY'
      call gform(ia(i_en),ia(i_gradE),ia(i_fluxlE),ia(i_x),ia(i_sedgeE)
     .            ,ia(i_e),ia(i_ls),ia(i_ie),ia(i_nelcon),ia(i_pedgeE)
     .            ,ia(i_ix),numel,ndm,nen,nen,1,ndfE,2,1,3)
      call uformnode(ia(i_un),ia(i_en),ia(i_gradE),ia(i_fluxlE),ia(i_x)
     .             ,ia(i_mdf),ia(i_ix),ia(i_md),nnode,numel,ndm,nen,1
     .             ,3)
      fileout = name(prename,istep,13)
      nCell = 'elEnergy'
      nNod  = 'noEnergy'
      call write_res_vtk(ia(i_ix),ia(i_x),ia(i_en),ia(i_un),nnode,numel
     .                  ,ndm,nen,1,fileout,nCell,nNod ,bvtk,4,t,istep
     .                  ,nout)
      goto 50
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando:       
c
c ......................................................................
 1400 continue
      print*, 'Macro VMASSESPE'
      vMass = .true.
      write(*,*)' Set vMass ',vMass
      goto 50
      goto 50
c ......................................................................
c
c ... Macro-comando: PMASSESP
c
c ......................................................................
 1500 continue
      print*, 'Macro PMASSEPS'
      i_nn        = alloc_8('nn      ',2,nnode)
      call uformnode(ia(i_nn),ia(i_ro),ddum,ddum,ia(i_x),ia(i_mdf)
     .              ,ia(i_ix),ia(i_md),nnode,numel,ndm,nen,3,2)
      fileout = name(prename,istep,14) 
      nCell = 'elMassEP'
      nNod  = 'noMassEP'
      call write_res_vtk(ia(i_ix),ia(i_x),ia(i_ro),ia(i_nn),nnode,numel
     .                  ,ndm,nen,3,fileout,nCell,nNod ,bvtk,8,t
     .                  ,istep,nout)
      i_nn    = dealloc('nn      ')
      goto 50
c ......................................................................
c
c ... Macro-comando:PTEMP  
c
c ......................................................................
 1600 continue
      print*, 'Macro PTEMP'
c ...
      call gform(ia(i_en),ia(i_gradE),ia(i_fluxlE),ia(i_x)
     .          ,ia(i_sedgeE),ia(i_e),ia(i_ls),ia(i_ie),ia(i_nelcon)
     .          ,ia(i_pedgeE),ia(i_ix),numel,ndm,nen,nen,1,ndfF,2,3,2)
c .......................................................................
c
c ...
      call gradEtoGradT(ia(i_gradE),ia(i_gradT),ia(i_e),ia(i_ix)
     .                 ,numel,ndm,nen)
c .......................................................................
c
c ...     
      call uformnode(ia(i_un),ia(i_temp),ia(i_gradT),ia(i_fluxlE)
     .              ,ia(i_x)
     .              ,ia(i_mdf),ia(i_ix),ia(i_md),nnode,numel,ndm,nen,1
     .              ,3)
c .......................................................................
c
c ...
      fileout = name(prename,istep,3)
      nCell = 'elTemp'
      nNod  = 'noTemp'
      call write_res_vtk(ia(i_ix),ia(i_x),ia(i_temp),ia(i_un),nnode
     .                  ,numel
     .                  ,ndm,nen,1,fileout,nCell,nNod ,bvtk,4,t,istep
     .                  ,nout)
c ........................................................................
      goto 50
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando: PDIV    
c
c ......................................................................
 1700 continue
      print*, 'Macro PDIV'
      call uformnode(ia(i_un),ia(i_div),ddum,ddum,ia(i_x),ia(i_mdf)
     .              ,ia(i_ix),ia(i_md),nnode,numel,ndm,nen,1,2)
      fileout = name(prename,istep,15)
      nCell = 'elDiv'
      nNod  = 'noDiv' 
      call write_res_vtk(ia(i_ix),ia(i_x),ia(i_div),ia(i_un),nnode,numel
     .                  ,ndm,nen,1,fileout,nCell,nNod ,bvtk,4,t
     .                  ,istep,nout)
      goto 50
c ......................................................................
c
c ... Macro-comando: TSIMPU1
c
c ......................................................................
 1800 continue
      print*, 'Macro TSIMPU1'
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =1810,end =1810) solvSimpleU1
      write(*,'(a,es10.2)')' Set solvSimpleU1 tol for ',solvSimpleU1
      goto 50
 1810 continue
      print*,'Erro na leitura da macro (TSIMPU1) !'
      goto 5000
c ......................................................................
c
c ... Macro-comando: TSIMPU2
c
c ......................................................................
 1900 continue
      print*, 'Macro TSIMPU2'
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =1910,end =1910) solvSimpleU2
      write(*,'(a,es10.2)')' Set solvSimpleU2 tol for ',solvSimpleU2
      goto 50
 1910 continue
      print*,'Erro na leitura da macro (TSIMPU2) !'
      goto 5000
c ......................................................................
c
c ... Macro-comando: TSIMPPC
c
c ......................................................................
 2000 continue
      print*, 'Macro TSIMPPC'
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =2010,end =2010) solvSimplePc
      write(*,'(a,es10.2)')' Set solvSimplePc tol for ',solvSimplePc
      goto 50
 2010 continue
      print*,'Erro na leitura da macro (TSIMPPC) !'
      goto 5000
c ......................................................................
c
c ... Macro-comando: TSIMPEN
c
c ......................................................................
 2100 continue
      print*, 'Macro TSIMPEN'
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =2110,end =2110) solvEnergy
      write(*,'(a,es10.2)')' Set  solvEnergy tol for  ', solvEnergy
      goto 50
 2110 continue
      print*,'Erro na leitura da macro (TSIMPEN) !'
      goto 5000
c ......................................................................
c
c ... Macro-comando: TDINAMIC
c
c ...................................................................... 
 2200 continue
      print*, 'Macro TDINAMIC'
      tDinamico = .true.
      write(*,*)' Set tDinamico ',tDinamico
      goto 50
c ......................................................................
c
c ... Macro-comando: TDINAMIC
c
c ......................................................................          
 2300 continue
      print*, 'Macro MAXITT1' 
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =2310,end =2310) maxItT1
      write(*,'(a,i10)')' Set max T1 Transport it for ', maxItT1
      goto 50
 2310 continue
      print*,'Erro na leitura da macro (MAXITT1) !'
      goto 5000
c ......................................................................
c
c ... Macro-comando: PRESMASS     
c
c ......................................................................
 2400 continue
      fileout = name(prename,istep,16)
      nCell = 'elResMass'
      nNod  = 'noResMass' 
      call write_res_vtk(ia(i_ix),ia(i_x),ia(i_rCellPc),ddum,nnode
     .                  ,numel,ndm,nen,1,fileout,nCell,nNod,bvtk,9,t
     .                  ,istep,nout)
      goto 50
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando: SAVE         
c
c ......................................................................
 2500 continue
      print*, 'Macro SAVE '       
      fileout = name(prename,0,30)
      call saveSimple(ia(i_w),ia(i_w0),ia(i_p),ia(i_en),ia(i_en0)
     .               ,ia(i_ro),nnode,numel,ndm,istep,t
     .               ,fileout,.false.,sEnergy,noutSave)
      goto 50
c ......................................................................
c
c ... Macro-comando: load         
c
c ......................................................................
 2600 continue
      print*, 'Macro LOAD '       
      fileout = name(prename,0,30)
      call saveSimple(ia(i_w),ia(i_w0),ia(i_p),ia(i_en),ia(i_en0)
     .               ,ia(i_ro),nnode,numel,ndm,istep,t
     .               ,fileout,.true.,sEnergy,noutSave)
      goto 50
c ----------------------------------------------------------------------
c/ro
c ... Macro-comando:
c
c ......................................................................
 2700 continue
      print*, 'Macro ITSIMPLE    '
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =2710,end =2710) maxItSimple
      write(*,'(a,i10)')' Set max simple/simpleC it for '
     .                  , maxItSimple
      goto 50
 2710 continue
      print*,'Erro na leitura da macro (ITSIMPLE) !'
      goto 5000
c ----------------------------------------------------------------------
c
c ... Macro-comando:
c
c ......................................................................
 2800 continue
      print*, 'Macro TTRANST1'
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =2810,end =2810) solvT1
      write(*,'(a,es10.2)')' Set noliner trans T1 tol for  ',solvT1
      goto 50
 2810 continue
      print*,'Erro na leitura da macro (TTRANST1) !'
      goto 5000
c ----------------------------------------------------------------------
c
c ... Macro-comando: SETPNODE impressao de grandezas por no no tempo
c
c ......................................................................
 2900 continue
      print*, 'Macro SETPNODE'
      call readmacro(nin,.false.)
      write(str,'(80a)') (word(i),i=1,80)
      read(str,*,err=2910,end = 2910) pnodename
      goto 2920
c ... problema no arquivo auxiliar        
 2910 continue
      print*,'Erro na leitura da macro (SETPNODE)'
      flag_pcd = .false.
      goto 2930
c ... leitura normal 
 2920 continue     
      call readpnode(pnodename,i_node,i_nfile,num_pnode,flag_pcd,nout)
      new_file(1:nfiles) = .true.
 2930 continue
c ... erro na letura do nome do arquivo auxiliar      
      if( flag_pcd .eqv. .false.) stop
      goto 50
c .....................................................................
c 
c .....................................................................
c
c ... Macro-comando: PNTEMP impressao da temperatura por no no tempo
c     (SETPNODE)                                                   
c ......................................................................
 3000 continue
      print*, 'Macro PNTEMP    '
      if(flag_pcd.eqv..false.) then
        print*,'Nemhum no de impressao para PNTEMP!'
        stop 
      endif  
c ... codigo para o arquivo _temp_node.txt      
      code = 25
c .....................................................................
      string = 'Temperatura'
      call uformnode(ia(i_un),ia(i_t1),ia(i_gradT1),ia(i_x),ia(i_mdf)
     .              ,ia(i_ix),ia(i_md),nnode,numel,ndm,nen,ndfT1,2)
      do j = 1, num_pnode
         call printnode(ia(i_un),ia(i_node+j-1),ndfT1,istep,t
     .                ,string,prename,ia(i_nfile+j-1),code,new_file(1))
      enddo
      new_file(1) = .false.
      goto 50
c ----------------------------------------------------------------------
c
c .....................................................................
c
c ... Macro-comando: UNDERPC                                         
c                                                                  
c ...................................................................... 
 3100 continue
      print*, 'Macro UNDERPC'
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =2810,end =2810) underPc
      write(*,'(a,d10.2)')' Set underPc for ',underPc
      goto 50
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: SETPCELL impressao de grandezas por no no tempo
c
c ......................................................................
 3200 continue
      print*, 'Macro SETPCELL'
      call readmacro(nin,.false.)
      write(str,'(80a)') (word(i),i=1,80)
      read(str,*,err=3310,end = 3310) pcellname
      goto 3320
c ... problema no arquivo auxiliar        
 3310 continue
      print*,'Erro na leitura da macro (SETPCEEL)'
      flag_pcd = .false.
      goto 3330
c ... leitura normal 
 3320 continue     
      call readpnode(pcellname,i_cell,i_nfile,num_pcell,flag_pcd,nout)
      new_file(1:nfiles) = .true.
 3330 continue
c ... erro na letura do nome do arquivo auxiliar      
      if( flag_pcd .eqv. .false.) stop
      goto 50
c .....................................................................
c
c ... Macro-comando: PCTEMP impressao da temperatura por no no tempo
c     (SETPCELL)                                                   
c ......................................................................
 3300 continue
      print*, 'Macro PCTEMP    '
      if(flag_pcd.eqv..false.) then
        print*,'Nemhum no de impressao para PCTEMP!'
        stop 
      endif  
c ... codigo para o arquivo _temp.txt      
      code = 27
c .....................................................................
      string = 'Temperatura'
      do j = 1, num_pcell
        call printnode(ia(i_t1),ia(i_cell+j-1),ndfT1,0,0.0d0
     .                ,string,prename,ia(i_nfile+j-1),code,new_file(1))
      enddo
      new_file(1) = .false.
      goto 50
c ----------------------------------------------------------------------
c
c ... Macro-comando: PCGRAD impressao da temperatura por no no tempo
c     (SETPCELL)                                                   
c ......................................................................
 3400 continue
      print*, 'Macro PCGRAD    '
      if(flag_pcd.eqv..false.) then
        print*,'Nemhum no de impressao para PCGRAD!'
        stop 
      endif  
c ... codigo para o arquivo _grad.txt      
      code = 28
c .....................................................................
      string = 'Gradiente'
      do j = 1, num_pcell
        call printnode(ia(i_gradT1),ia(i_cell+j-1),ndm,0,0.0d0
     .                ,string,prename,ia(i_nfile+j-1),code,new_file(2))
      enddo
      new_file(2) = .false.
      goto 50
c .....................................................................
c
c......................................................................
c
c ... Macro-comando: SIMPLEV
c
c ......................................................................
 3500 continue
      print*, 'Macro SIMPLEC'
      simpleC = .true.
      write(*,*)' Set SimpleC ',simpleC
      goto 50
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando: UNDERP 
c
c ......................................................................
 3600 continue
      print*, 'Macro UNDERP'
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =2810,end =2810) underP
      write(*,'(a,e10.2)')' Set underP for ',underP
      goto 50
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando: UNDERU 
c
c ......................................................................
 3700 continue
      print*, 'Macro UNDERU'
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =2810,end =2810) underU
      write(*,'(a,e10.2)')' Set underU for  ',underU
      goto 50
c ......................................................................
c
c ... Macro-comando: SETELEV
c
c ......................................................................
 3800 continue
      print*, 'Macro SETELEV   '
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err = 3810,end = 3810) setelev
      goto 50
 3810 continue
      print*,'Erro na leitura da macro (SETELEV) !'
      goto 5000      
c ......................................................................
c
c ... Macro-comando: DT
c
c ......................................................................
 3900 continue
      print*, 'Macro DT   '
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err = 3910,end = 3910) dt
      write(*,'(a,es10.2)')' Set DT for ',dt
      goto 50
 3910 continue
      print*,'Erro na leitura da macro (DT) !'
      goto 5000      
c ----------------------------------------------------------------------
c
c ......................................................................
c
c ... Macro-comando: CONFIG
c
c ......................................................................
 4000 continue
      print*, 'Macro CONFIG   '
      call readConfig(maxmem,openMpCell,openMpSolver
     .              ,nThreadsCell,nThreadsSolver,freord,bvtk,nin)
      goto 50
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando: PGRADU1
c
c ......................................................................
 4100 continue
      print*, 'Macro PGRADU1'
      fileout     = name(prename,istep,31) 
      i_nn        = alloc_8('nn      ',ndm,nnode)
      call uformnode(ia(i_nn),ia(i_gradU1),ddum,ddum,ia(i_x),ia(i_mdf)
     .              ,ia(i_ix),ia(i_md),nnode,numel,ndm,nen,ndm,2)
      nCell = 'elGradU1'
      nNod  = 'noGradU1'
      call write_res_vtk(ia(i_ix),ia(i_x),ia(i_gradU1),ia(i_nn),nnode
     .                  ,numel,ndm,nen,1,fileout,nCell,nNod ,bvtk,5,t
     .                  ,istep,nout)
      i_nn        = dealloc('nn      ')
      goto 50
 4110 continue
      print*,'Erro na leitura da macro (PGRADU1) !'
      goto 5000
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando: PGRADU2
c
c ......................................................................
 4200 continue
      print*, 'Macro PGRADU2'
      fileout     = name(prename,istep,32) 
      i_nn        = alloc_8('nn      ',ndm,nnode)
      call uformnode(ia(i_nn),ia(i_gradU2),ddum,ddum,ia(i_x),ia(i_mdf)
     .              ,ia(i_ix),ia(i_md),nnode,numel,ndm,nen,ndm,2)
      nCell = 'elGradU2'
      nNod  = 'noGradU2'
      call write_res_vtk(ia(i_ix),ia(i_x),ia(i_gradU2),ia(i_nn),nnode
     .                  ,numel,ndm,nen,1,fileout,nCell,nNod ,bvtk,5,t
     .                  ,istep,nout)
      i_nn        = dealloc('nn      ')
      goto 50
 4210 continue
      print*,'Erro na leitura da macro (PGRADU2) !'
      goto 5000
c ......................................................................
c
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando: PGRADU3
c
c ......................................................................
 4300 continue
      print*, 'Macro '
      goto 50
 4310 continue
      print*,'Erro na leitura da macro () !'
      goto 5000
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando: SKEWNESSCORRECTION
c
c ......................................................................
 4400 continue
      print*, 'Macro SKEWNESSCORRECTION'
      skewnessCorrection = .true.
      write(*,*)' Set SKEWNESSCORRECTION ',skewnessCorrection
      goto 50
 4410 continue
      print*,'Erro na leitura da macro (SKEWC) !'
      goto 5000
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando: PGEO
c
c ......................................................................
 4500 continue
      print*, 'Macro PGEO'
      fileout = name(prename,istep,33)
      i_pedgeT1 = 1
      i_sedgeT1 = 1
      call write_geo_vtk(ia(i_ix),ia(i_x)
     .                  ,ia(i_pedgeF),ia(i_pedgeE),ia(i_pedgeT1)
     .                  ,ia(i_sedgeF),ia(i_sedgeE),ia(i_sedgeT1)
     .                  ,nnode,numel,ndm,nen,ndfT1,ndfF,ndfE
     .                  ,fileout,bvtk,nout)
      goto 50
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando: PGRADT 
c
c ......................................................................
 4600 continue
      print*, 'Macro PGRADT'
c ...
      fileout     = name(prename,istep,34) 
      i_nn        = alloc_8('nn      ',ndm,nnode)
c .....................................................................
c
c ...
      call gform(ia(i_en),ia(i_gradE),ia(i_fluxlE),ia(i_x)
     .          ,ia(i_sedgeE),ia(i_e),ia(i_ls),ia(i_ie),ia(i_nelcon)
     .          ,ia(i_pedgeE),ia(i_ix),numel,ndm,nen,nen,1,ndfF,2,3,2)
c .......................................................................
c
c ...
       call gradEtoGradT(ia(i_gradE),ia(i_gradT),ia(i_e),ia(i_ix)
     .                 ,numel,ndm,nen)
c .......................................................................
c
c ...
       call uformnode(ia(i_nn),ia(i_gradT),ddum,ddum,ia(i_x),ia(i_mdf)
     .               ,ia(i_ix),ia(i_md),nnode,numel,ndm,nen,ndm,2)
c .......................................................................
c
c ...
      nCell = 'elGradT'
      nNod  = 'noGradT'
      call write_res_vtk(ia(i_ix),ia(i_x),ia(i_gradT),ia(i_nn),nnode
     .                  ,numel,ndm,nen,1,fileout,nCell,nNod ,bvtk,5,t
     .                  ,istep,nout)
      i_nn        = dealloc('nn      ')
      goto 50
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando:        
c
c ......................................................................
 4700 continue
      print*, 'Macro PGRADE'
      fileout     = name(prename,istep,35) 
      i_nn        = alloc_8('nn      ',ndm,nnode)
      call gform(ia(i_en),ia(i_gradE),ia(i_fluxlE),ia(i_x)
     .          ,ia(i_sedgeE),ia(i_e),ia(i_ls),ia(i_ie),ia(i_nelcon)
     .          ,ia(i_pedgeE),ia(i_ix),numel,ndm,nen,nen,1,ndfF,2,3,2)
      call uformnode(ia(i_nn),ia(i_gradE),ddum,ddum,ia(i_x),ia(i_mdf)
     .              ,ia(i_ix),ia(i_md),nnode,numel,ndm,nen,ndm,2)
      nCell = 'elGradE'
      nNod  = 'noGradE'
      call write_res_vtk(ia(i_ix),ia(i_x),ia(i_gradE),ia(i_nn),nnode
     .                  ,numel,ndm,nen,1,fileout,nCell,nNod ,bvtk,5,t
     .                  ,istep,nout)
      i_nn        = dealloc('nn      ')
      goto 50
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando:        
c
c ......................................................................
 4800 continue
      print*, 'Macro 4800'
      goto 50
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando:        
c
c ......................................................................
 4900 continue
      print*, 'Macro 4900'
      goto 50
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando:        
c
c ......................................................................
 5000 continue
      print*, 'Macro 5000'
      goto 50
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando: pvort2D
c
c ......................................................................
 5100 continue
      print*, 'Macro PVORT2D'
c ...
      call guess(ia(i_u1),ia(i_u2),ia(i_w),numel,ndfF-1)
c ......................................................................
c
c ...
      call gform(ia(i_u1),ia(i_gradU1),ia(i_fluxlU1),ia(i_x)
     .          ,ia(i_sedgeF)
     .          ,ia(i_e),ia(i_ls),ia(i_ie),ia(i_nelcon),ia(i_pedgeF)
     .          ,ia(i_ix),numel,ndm,nen,nen,1,ndfF,2,1,1)
c
      call gform(ia(i_u2),ia(i_gradU2),ia(i_fluxlU2),ia(i_x)
     .          ,ia(i_sedgeF)
     .          ,ia(i_e),ia(i_ls),ia(i_ie),ia(i_nelcon),ia(i_pedgeF)
     .          ,ia(i_ix),numel,ndm,nen,nen,1,ndfF,2,2,1)
c ......................................................................
c
c ...
      i_nn        = alloc_8('nn      ',1,numel)
      call vort2D(ia(i_gradU1),ia(i_gradU2),ia(i_nn),numel,ndm)
c ...
c
c ...
      call uformnode(ia(i_un),ia(i_nn),ddum,ddum,ia(i_x),ia(i_mdf)
     .              ,ia(i_ix),ia(i_md),nnode,numel,ndm,nen,1,2)
c ......................................................................
c
c ...      
      fileout     = name(prename,istep,36)
      nCell = 'elVort'
      nNod  = 'noVort'
      call write_res_vtk(ia(i_ix),ia(i_x),ia(i_nn),ia(i_un),nnode
     .                  ,numel
     .                  ,ndm,nen,1,fileout,nCell,nNod ,bvtk,4,t,istep
     .                  ,nout) 
      i_nn        = dealloc('nn      ')
c ......................................................................
c
c ...
      goto 50
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando:        
c
c ......................................................................
 5200 continue
      print*, 'Macro 5200'
      goto 50
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando:        
c
c ......................................................................
 5300 continue
      print*, 'Macro 5300'
      goto 50
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando:        
c
c ......................................................................
 5400 continue
      print*, 'Macro 5400'
      goto 50
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando:        
c
c ......................................................................
 5500 continue
      print*, 'Macro 5500'
      goto 50
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando: UNDERRO
c
c ......................................................................
 5600 continue
      print*, 'Macro UNDERRO'
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =5610,end =5610) underRo
      write(*,'(a,e10.2)')' Set underU for  ',underRo
      goto 50
 5610 continue
      print*,'Erro na leitura da macro (UNDERRO) !'
      goto 5000
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando: TOLPCG
c
c ......................................................................
 5700 continue
      print*, 'Macro TOLPCG'
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =5710,end =5710) solvtolPcg
      write(*,'(a,e10.2)')' Set tol PCG for   ',solvtolPcg
      goto 50
 5710 continue
      print*,'Erro na leitura da macro (TOLPCG) !'
      goto 5000
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando: TOLBIPCG
c
c ......................................................................
 5800 continue
      print*, 'Macro TOLBIPCG'
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =5810,end =5810) solvtolBcg
      write(*,'(a,e10.2)')' Set tol BIPCG for ',solvtolBcg
      goto 50
 5810 continue
      print*,'Erro na leitura da macro (TOLBIPCG) !'
      goto 5000
c ......................................................................
c
c ......................................................................
c
c ... Macro-comando: MAXITSOL
c
c ......................................................................
 5900 continue
      print*, 'Macro MAXITSOL'
      call readmacro(nin,.false.)
      write(string,'(30a)') (word(i),i=1,30)
      read(string,*,err =5910,end =5910) maxItSol
      write(*,'(a,i10)')' Set MaxIt solver(s) for ',maxItSol
      goto 50
 5910 continue
      print*,'Erro na leitura da macro (MAXITSOL) !'
      goto 5000
c ......................................................................
c
c ......................................................................                                                                        
c
c ... Macro-comando: STOP   
c
c ......................................................................
 6090 continue
      print*, 'Salvando ...'       
      fileout = name(prename,0,30)
      call saveSimple(ia(i_w),ia(i_w0),ia(i_p),ia(i_en),ia(i_en0)
     .               ,ia(i_ro),nnode,numel,ndm,istep,t
     .               ,fileout,.false.,sEnergy,noutSave)
      print*, 'Salvo.'
c ......................................................................        
 6000 continue
      close(3)
      close(10)
      close(nin)
      totaltime = get_time() - totaltime
      call write_log(ntime,prename,neqU1,neqU2,neqPc,neqE
     .              ,nadU1,nadU2,nadPc,nadE,bandaU1,bandaU2,bandaPc
     .              ,bandaE,cfl,re,prandtl,grashof,Massa0,Massa,FluxoM
     .              ,vol)
      call common_finalize()
      end
c **********************************************************************
c
c **********************************************************************
      subroutine s(n,m)
      implicit none
      integer n(4,*),i,m
      do i=1, m
        print*,n(1,i),n(2,i),n(3,i),n(4,i),i
      enddo
      return
      end
c **********************************************************************
c **********************************************************************
      subroutine si(n,m,nl,nout)
      implicit none
      real*8 n(nl,*)
      integer i,j,m,nl
      integer nout
      do i=1, m
        write(nout,'(i5,99f20.8)'),i,(n(j,i),j=1,nl)
      enddo
      return
      end
c **********************************************************************
      subroutine wi(v,nl)
      implicit none
      real*8 v(2,*)
      integer i,nl
      do i=1, nl
        v(1,i) = 1.0d0
        v(2,i) = 1.0d0
      enddo
      return
      end
c *********************************************************************
c   
c **********************************************************************
      subroutine sistema(ad,a,b,ia,ja,neq,nout)
      implicit none
      real*8 ad(*),a(*),b(*)
      integer ia(*),ja(*)
      integer i,j,neq
      integer nout
      do i=1, neq
        write(nout,'(i,f16.8,4f16.8,f16.8)')i,ad(i)
     .       ,(a(j),j=ia(i),ia(i+1)-1),b(i)
      enddo
      return
      end
c *********************************************************************
