c *********************************************************************
c * SIMPLE: Algoritmo SIMPLE/SIMPLEC para o acoplamento da pressao/   *
c * velocidade                                                        *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * x        - coordenadas nodais                                     *
c * ix       - conectividade dos elementos                            *
c * e        - propriedade                                            *
c * ie       - tipo do elemento do material                           *
c * nelcon   - vizinho por celula                                     *
c * pedge    - tipo de condicao de contorno u1,u2,p                   *
c * sedge    - condicoes de contorno nas arestas u1,u2,p              *
c * pedgeE   - tipo de condicao de contorno E                         *
c * sedgeE   - condicoes de contorno nas arestas  E                   *
c * w        - campo de velocidade por celula (n-1)                   *
c * w0       - campo de velocidade por celula (n-2)                   *
c * wP       - campo de velocidade por celula auxiliar                *
c * num      - numeracao dos celulas                                  *
c * ls       - reconstrucao de gradiente least square                 *
c * gradU1   - gradiente u1 da celula                                 *
c * gradU2   - gradiente u2 da celula                                 *
c * gradP    - gradiente p da celula                                  *
c * gradPc   - gradiente pC da celula                                 *
c * gradE    - gradiente E da celula                                  *
c * div      - divergente da velocidade                               *
c * mP       - parametros da equacao de momentos(cfl,Reynalds)        *
c * iM       - interpolaco do momentos                                *
c * fluxlU1  - limitador do fluxo                                     *
c * fluxlU2  - limitador do fluxo                                     *
c * fluxlPc  - limitador do fluxo                                     *
c * fluxlE   - limitador do fluxo                                     *
c * rCu1     - residuo o celula u1                                    *
c * rCu2     - residuo o celula u2                                    *
c * rCPc     - residuo o celula Pc                                    *
c * rCE      - residuo o celula E                                     *
c * adU1     - matriz Au                                              *
c * auU1     - matriz Au                                              *
c * alU1     - matriz Au                                              *
c * bU1      - vetor de forcas bU1                                    *
c * bU10     - vetor de forcas bU1                                    *
c * iaU1     - ponteiro do para o csr                                 *
c * jaU1     - estrutura do csr                                       *
c * adU2     - matriz Av                                              *
c * auU2     - matriz Av                                              *
c * alU2     - matriz Av                                              *
c * bU2      - vetor de forcas bU2                                    *
c * bU20     - vetor de forcas bU2                                    *
c * iaU2     - ponteiro do para o csr                                 *
c * jaU2     - estrutura do csr                                       *
c * adPc     - matriz ApC                                             *
c * auPc     - matriz ApC                                             *
c * alPc     - matriz ApC                                             *
c * bPc      - vetor de forcas bPc                                    *
c * iaPc     - ponteiro do para o csr                                 *
c * jaPc     - estrutura do csr                                       *
c * adE      - matriz Ae                                              *
c * auE      - matriz Ae                                              *
c * alE      - matriz Ae                                              *
c * bE       - vetor de forcas bE                                     *
c * bE0      - vetor de forcas bE                                     *
c * iaE      - ponteiro do para o csr                                 *
c * jaE      - estrutura do csr                                       *
c * u1       - velocidade u1                                          *
c * u2       - velocidade u2                                          *
c * p        - pressao                                                *
c * pC       - pressao de correcao                                    *
c * pC1      - segunda pressao de correcao                            *
c * en       - energia                                                *
c * en0      - energia 0                                              *
c * r0       - massa especifica variavel                              *
c * un       - nao definido variavel qualquer por no                  *
c * mdf      - usado para interpolacap simples das variaveis          *
c * md       - usado para interpolacap simples das variaveis          *
c * ddu      - campo d usado no simple                                *
c * temp     - temperatura                                            *
c * sx       - arronjo para guardar a solucao do sistemas lineares    *
c * eddyVisc - viscoside turbulenta                                   *
c * cs       - coeficiente smagorinsk dinamico                        *
c * yPlus    - distancia adiminensional a parede                      *
c * numel  - numero de celula                                         *
c * ndm    - numero de dimensoes                                      *
c * nen    - numero de nos por celula                                 *
c * nshared- numero de celula compartilhados pela celula central      *
c * ndf    - graus de liberdade                                       *
c * ndfE   - graus de liberdade da equacao de energia                 *
c * dt     - passo de tempo                                           *
c * matrizU1- estrutura da matriz A ( 1 - CSR )                       *
c * matrizU2- estrutura da matriz A ( 1 - CSR )                       *
c * matrizPc- estrutura da matriz A ( 2 - CSR )                       *
c * matrizE - estrutura da matriz A ( 2 - CSR )                       *
c * solverU1-tipo do solver         ( 2 - PBiCGSTAB)                  *
c * solverU2-tipo do solver         ( 2 - PBiCGSTAB)                  *
c * solverPc-tipo do solver         ( 1 - PCG      )                  *
c * solverE -tipo do solver         ( 2 - PBiCGSTAB)                  *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * w      - campo de velocidades                                     *
c * p      - campo de pressão                                         *
c * ----------------------------------------------------------------- *
c *********************************************************************
       subroutine simple(x         ,ix             ,e        ,ie
     .                  ,nelcon    ,pedge          ,sedge    ,pedgeE
     .                  ,sedgeE    ,w              ,w0       ,wP       
     .                  ,num       ,ls             ,gradU1   ,gradU2   
     .                  ,gradP     ,gradPc         ,gradE    ,div      
     .                  ,mP        ,iM             ,fluxlU1  ,fluxlU2  
     .                  ,fluxlPc   ,fluxlE         ,rCu1     ,rCu2
     .                  ,rCpC      ,rCe            ,adU1     ,auU1 
     .                  ,alU1      ,bU1       ,bU10          ,iaU1    
     .                  ,jaU1      ,adU2      ,auU2          ,alU2     
     .                  ,bU2       ,bU20      ,iaU2          ,jaU2     
     .                  ,adPC      ,auPc      ,alPc          ,bPc      
     .                  ,iaPc      ,jaPc      ,adE           ,auE    
     .                  ,alE       ,bE        ,bE0           ,iaE     
     .                  ,jaE       ,u1        ,u2            ,p   
     .                  ,pC        ,pC1       ,en            ,en0
     .                  ,rO        ,un        ,mdf           ,md
     .                  ,ddU       ,temp      ,sx            ,eddyVisc
     .                  ,cs        ,yPlus     
     .                  ,numel     ,nnode     ,ndm      
     .                  ,nen       ,nshared   ,ndf           ,ndfE  
     .                  ,dt        ,t         ,matrizU1      ,matrizU2 
     .                  ,matrizPc  ,matrizE   ,neqU1         ,neqU2   
     .                  ,neqPc     ,neqE      ,nadU1         ,nadU2  
     .                  ,nadPc     ,nadE      ,solverU1      ,solverU2 
     .                  ,solverPc  ,solverE   ,solvTolPcg    ,solvTolBcg
     .                  ,maxIt    ,noutSimple,itResSimplePlot,sEnergy   
     .                  ,istep    ,cfl       ,reynalds       ,kEnergy 
     .                  ,prandtl  ,disfiltro ,transSubGrid
     .                  ,grashof  ,vol       ,unsymPc        ,bs
     .                  ,prename  ,noutCoo   ,flagTurbulence ,flagCoo)
c...
      implicit none
      include 'time.fi'
      include 'simple.fi'
      include 'openmp.fi'
      real(8)  x(*),e(*),sedge(*),sedgeE(*),w(ndm,*),w0(ndm,*)
      real(8)  ro(*),temp(*),un(*),mdf(*),wP(ndm,*)
      integer ix(*),nelcon(*),ie(*),pedge(*),pedgeE(*),num(*),ls(*)
      integer md(*)
      logical unsymPc
      logical bs
c ... turbulencia 
      real*8 eddyVisc(*),cs(*),yPlus(*)
c ...
      character*80 prename
c ... gradientes
      real(8) gradU1(ndm,*),gradU2(ndm,*),gradP(*),gradPc(*),gradE(*)
      real(8) div(*)
      real(8) fluxlU1(*),fluxlU2(*),fluxlPc(*),fluxlE(*),iM(*),mP(*)
c ... sistema
      real(8) rCu1(*),rCu2(*),rCPc(*),rCe(*)
      real(8) adU1(*),auU1(*),alU1(*),bU1(*),bU10(*)
      real(8) adU2(*),auU2(*),alU2(*),bU2(*),bU20(*)
      real(8) adE(*),auE(*),alE(*),bE(*),bE0(*)
      real(8) adPc(*),auPc(*),alPc(*),bPc(*)
      integer iaU1(*),jaU1(*),iaU2(*),jaU2(*)
      integer iaPc(*),jaPc(*),iaE(*),jaE(*)
      integer matrizU1,matrizU2,matrizPc,matrizE
      integer neqU1,neqU2,neqPc,neqE,nadU1,nadU2,nadPc,nadE
c ... variaveis
      real(8)  u1(*),u2(*),p(*),pC(*),pC1(*),en(*),en0(*),ddU(*),dt,t
      integer numel,ndm,nen,nshared,ndf,ndfE,istep,nnode
c ... Solver      
      real(8)  solvtolPcg,solvtolBcg 
      integer solverU1,solverU2,solverPc,solverE,maxIt
      real(8)  sx(*)
c ...
      real(8) dot,rU1,rU2,rPc,rE,rU10,rU20,rE0,rPc0
      integer itSimple,conv
      integer noutSimple,noutCoo
      logical itResSimplePLot,solvU1,solvU2,solvE,flagCoo
      logical ResAbs,sEnergy,xMomentum,yMomentum,flagTurbulence
      integer i
c ... variaveis
      real(8) cfl,reynalds,prandtl,grashof,vol,kEnergy
      real(8) disfiltro,transSubGrid
      real(8) ZERO,kZero
      parameter (ZERO=1.0d-32)
c ...
      kZero         = 2
      ResAbs        = .false.
c ... tecnica de interpolação de velocidade nas faces (checkerboard problem) 
      intVel        = 5
c ... interpolacao da pressao nas faces de segunda ordem (momentum)    
      sPressure     = .true.
c .....................................................................
c
c ...       
      solvU1       = .true.
      solvU2       = .true.
      solvE        = .true.
c .....................................................................
c
c ...
       call azero(Pc,neqPc)
       call azero(gradPc,numel*ndm)
c .....................................................................
c
c ...
c      open(19, file= 'matriz2.txt')
      rPc0 = 1.d0
      rU10 = 1.d0
      rU20 = 1.d0 
      rE0  = 1.d0
      conv = 0
c .....................................................................
c
c ... calculo de parametro por celula
      call cellParameter(mP      ,rO     ,w      ,x  ,eddyVisc
     .                  ,gradU1  ,gradU2
     .                  ,e       ,sedge  ,pedge  ,ie ,nelcon
     .                  ,ix      ,numel  ,ndm    ,nen,nshared 
     .                  ,ndf     ,dt     ,      7,  1)
      call calParameter(mP     ,div       ,numel         ,cfl ,reynalds
     .                 ,kEnergy,disfiltro,transSubGrid
     .                 ,prandtl,grashof  ,vol   ,Massa   ,fluxoM
     .                 ,dt     ,tDinamico)
c .....................................................................
c
c ...
      write(noutSimple,'(a,i6,1x,a,es20.8)')'Istep :',istep,'Time  :',t
      write(noutSimple,'(a,1x,es20.8)')'Time  :',dt
c .....................................................................
c
c ... u1 do passo n e n-1
      u1(1:neqU1) =  w(1,1:neqU1)
      if(bs) u2(1:neqU1) = w0(1,1:neqU1)
c ... F = fontes no volume + Vp*u0
      elmU1Time = get_time() - elmU1Time
      call pform(adU1,auU1,alU1,iaU1,jaU1,num,un,u1,u2,u1,gradU1,gradU2
     .           ,gradP,div,fluxlU1,iM,ro
     .           ,w,bU10,rCu1,x,sedge,e,ddU,eddyVisc
     .           ,ie,nelcon,pedge,ix
     .           ,numel,ndm,nen,nshared,1,ndf,dt,0.0d0
     .           ,matrizU1,.true.,3,1,1,.true.,.false.,.false.
     .           ,.false.,.false.,bs)
      elmU1Time = get_time() - elmU1Time
c .....................................................................
c
c ... u2 do passo n e n-1
      u1(1:neqU1) =  w(2,1:neqU2)
      if(bs) u2(1:neqU1) = w0(2,1:neqU2)
c ... F = fontes no volume + Vp*v0
      elmU2Time = get_time() - elmU2Time
      call pform(adU2,auU2,alU2,iaU2,jaU2,num,un,u1,u2,u1,gradU1,gradU2
     .           ,gradP,div,fluxlU2,iM,ro
     .           ,w,bU20,rCu2,x,sedge,e,ddU,eddyVisc
     .           ,ie,nelcon,pedge,ix
     .           ,numel,ndm,nen,nshared,1,ndf,dt,0.0d0
     .           ,matrizU2,.true.,3,2,1,.true.,.false.,.false.
     .           ,.false.,.false.,bs)
      elmU2Time = get_time() - elmU2Time
c .....................................................................
c
c ... w0 = w(n-1)
      call aequalbVetor(w0,w,numel,ndf-1)
c .....................................................................
c
c ... en = e(n) e en0 = e(n-1)
      if(sEnergy) then
c ...
c         call EnthalpyForTemp(en,temp,e,numel,.false.,1)
c .....................................................................
c
c ... F = fontes no volume + Vp*en0
        elmETime = get_time() - elmETime    
        call pform(adE,auE,alE,iaE,jaE,num,un,en,en0,en0,gradE,gradE
     .            ,gradE,div,fluxlE,iM,rO
     .            ,w,bE0,rCe,x,sedgeE,e,ddU,eddyVisc
     .            ,ie,nelcon,pedgeE,ix
     .            ,numel,ndm,nen,nshared,1,ndfE,dt,0.0d0
     .            ,matrizE,.true.,4,1,3,.true.,.false.,.false.
     .            ,.false.,.false.,bs)
        elmETime = get_time() - elmETime
c ... en0 = e(n)
        call aequalb(en0,en,neqE)
      endif 
c ......................................................................
c
c ...
c     if(.NOT. simpleR) then
c       call pressureGuess(p,pedge,sedge,numel,nshared,ndf)
c     endif
c .....................................................................
c
c ...
      simpleLoop: do itSimple = 1 ,  maxItSimple
        simpleItTime = get_time()
c ...     
        call guess(u1,u2,w,numel,ndf-1)
c .....................................................................
c       
c ... recontrucao do gradiente de u1
        grSimpleTime = get_time() - grSimpleTime
c       call uformnode(un,u1,ddum,ddum,x,mdf,ix,md,nnode,numel,ndm
c   .                ,nen,1,2)
        call gform(u1,gradU1,fluxlU1,x,sedge,e,ls,ie,nelcon,pedge,ix
     .            ,numel,ndm,nen,nshared,1,ndf,2,1,1)
        grSimpleTime = get_time() - grSimpleTime
c .....................................................................
c
c ... recontrucao do gradiente de u2
        grSimpleTime = get_time() - grSimpleTime
c       call uformnode(un,u2,ddum,ddum,x,mdf,ix,md,nnode,numel,ndm
c   .                ,nen,1,2)
        call gform(u2,gradU2,fluxlU2,x,sedge,e,ls,ie,nelcon,pedge,ix
     .            ,numel,ndm,nen,nshared,1,ndf,2,2,1)
        grSimpleTime = get_time() - grSimpleTime
c .....................................................................
c
c ... calculo de divergente
        call divergente(div,gradU1,gradU2,numel,ndf-1)
c .....................................................................
c
c ... recontrucao do gradiente da pressao
        grSimpleTime = get_time() - grSimpleTime 
c       call uformnode(un,p,ddum,ddum,x,mdf,ix,md,nnode,numel,ndm
c   .                 ,nen,1,2)
        call gform(p,gradP,fluxlPc,x,sedge,e,ls,ie,nelcon,pedge,ix
     .            ,numel,ndm,nen,nshared,1,ndf,2,3,2)
        grSimpleTime = get_time() - grSimpleTime
c .....................................................................
c
c ...   
        if(flagTurbulence) then
          turbLes = get_time() - turbLes 
          call turbulenceLes(w  ,gradU1,gradU2  , x   ,eddyVisc
     .                      ,cs ,yPlus                               
     .                      ,sedge,e  ,     ie, nelcon ,pedge ,ix
     .                      ,numel,ndm,   nen,nshared,ndf
     .                      ,nshared,8,0,1)
          turbLes = get_time() - turbLes 
        endif
c .....................................................................
c ... SIMPLER 
c        if(simpleR) then
c ... pseudo-velocidades
c ... u1
c          elmU1Time = get_time() - elmU1Time
c          call pform(adU1,auU1,alU1,iaU1,jaU1,num,un,u1,p,p,gradU1
c     .              ,gradU2,gradP,div,fluxlU1,iM,rO
c     .              ,w,bU1,bU1,x,sedge,e,ddU,ie,nelcon,pedge,ix
c     .              ,numel,ndm,nen,nshared,1,ndf,dt,0.0d0
c     .              ,matrizU1,.true.,1,1,1,.false.,.false.,.true.,bs)
c          elmU1Time = get_time() - elmU1Time
c .....................................................................
c
c ... u2
c          elmU2Time = get_time() - elmU2Time
c          call pform(adU2,auU2,alU2,iaU2,jaU2,num,un,u2,p,p,gradU2
c     .              ,gradU1,gradP,div,fluxlU2,iM,rO
c     .              ,w,bU2,bU1,x,sedge,e,ddU,ie,nelcon,pedge,ix
c     .              ,numel,ndm,nen,nshared,1,ndf,dt,0.0d0
c     .              ,matrizU2,.true.,1,2,1,.false.,.false.,.true.,bs)
c          elmU2Time = get_time() - elmU2Time
c .....................................................................
c
c ... 
c          call pseudoVel(wP,iM,numel,ndm)
c .....................................................................
c
c ... equacao da presssao
c          print*,'pressao'  
c          elmPcTime = get_time() - elmPCTime
c          call pform(adPc,auPc,alPc,iaPc,jaPc,num,un,p,p,p,gradP
c     .              ,gradP,gradP,div,fluxlPc,iM,rO
c     .              ,wP,bPc,rCPc,x,sedge,e,ddU,ie,nelcon,pedge,ix
c     .              ,numel,ndm,nen,nshared,1,ndf,dt,0.0d0
c     .              ,matrizPc,unsymPc,4,2,2,.true.,.true.,.false.,bs)
c          elmPcTime = get_time() - elmPCTime
c .....................................................................
c
c ... Solver Ap = bp
c          call getRes(p,sx,num,neqPc,numel,.true.)
c          solvPcTime = get_time() - solvPCTime
c          call solver(adPc,auPc,alPc,sx,bPc,iaPc,jaPc,neqPc,nadPc
c     .               ,unsymPc,solvTolPcg,maxIt,solverPc,matrizPc,2)
c          solvPcTime = get_time() - solvPCTime
c          call getRes(p,sx,num,neqPc,numel,.false.)
c .....................................................................
c
c ... recontrucao do gradiente da pressao
c          grSimpleTime = get_time() - grSimpleTime 
c          call gform(p,gradP,fluxlPc,x,sedge,e,ls,ie,nelcon,pedge,ix
c     .              ,numel,ndm,nen,nshared,1,ndf,3,3,2)
c          grSimpleTime = get_time() - grSimpleTime
c .....................................................................
c        endif
c .....................................................................
c
c ... quantidade de movimento u1
c 
c ... montagem do sistema u1
        print*,'quantidade de movimento u1'
        elmU1Time = get_time() - elmU1Time
        call pform(adU1,auU1,alU1,iaU1,jaU1,num,un,u1,p,p,gradU1,gradU2
     .            ,gradP,div,fluxlU1,iM,rO
     .            ,w,bU1,rCu1,x,sedge,e,ddU,eddyVisc
     .            ,ie,nelcon,pedge,ix
     .            ,numel,ndm,nen,nshared,1,ndf,dt,0.0d0
     .            ,matrizU1,.true.,1,1,1,.true.,.true.,.true.
     .            ,.true.,.true.,bs)
        elmU1Time = get_time() - elmU1Time
c .....................................................................
c
c ...
        call asumb(bU1,bU10,neqU1,bU1)
        call getRes(sx,bU10,num,neqU2,numel,.false.)
        call asumb(rCu1,sx,neqU1,rCu1)
c .....................................................................
c
c ...
        rU1 = dot(bU1,bU1,neqU1)
        xMomentum = .true.
        if( rU1 .lt. ZERO) xMomentum = .false.
        if(itSimple .eq. kZero .AND. xMomentum) then
          rU1 = dsqrt(dot(rCU1,rCU1,neqU1)) 
          rU10 = rU1
        endif
        solvU1 = .false. 
c .....................................................................
c
c ... montagem do sistema u2
        print*,'quantidade de movimento u2'
        elmU2Time = get_time() - elmU2Time    
        call pform(adU2,auU2,alU2,iaU2,jaU2,num,un,u2,p,p,gradU2,gradU1
     .            ,gradP,div,fluxlU2,iM,rO
     .            ,w,bU2,rCu2,x,sedge,e,ddU,eddyVisc
     .            ,ie,nelcon,pedge,ix
     .            ,numel,ndm,nen,nshared,1,ndf,dt,0.0d0
     .            ,matrizU2,.true.,1,2,1,.true.,.true.,.true.
     .            ,.true.,.true.,bs)
        elmU2Time = get_time() - elmU2Time
c .....................................................................
c
c ...
        call asumb(bU2,bU20,neqU2,bU2)
        call getRes(sx,bU20,num,neqU2,numel,.false.)
        call asumb(rCu2,sx,neqU2,rCu2)
c .....................................................................
c
c ...
        rU2 = dot(bU2,bU2,neqU2)
        yMomentum = .true.
        if( rU2 .lt. ZERO) yMomentum = .false.
        if( itSimple .eq. kZero .AND. YMomentum) then
          rU2  = dsqrt(dot(rCU2,rCU2,neqU2)) 
          rU20 = rU2
        endif
        solvU2 = .false. 
c .....................................................................
c
c ... Solver Au1 = b1
        if(xMomentum) then
          call getRes(u1,sx,num,neqU1,numel,.true.)
          solvU1Time = get_time() - solvU1Time 
          call solver(adU1,auU1,alU1,sx,bU1,iaU1,jaU1,neqU1,nadU1
     .               ,.true.,solvTolBcg,maxIt,solverU1,matrizU1,2)
          solvU1Time = get_time() - solvU1Time
          call getRes(u1,sx,num,neqU1,numel,.false.)
        endif   
c .....................................................................
c
c ... Solver Au2 = b2
        if(yMomentum) then
          call getRes(u2,sx,num,neqU2,numel,.true.)  
          solvU2Time = get_time() - solvU2Time 
          call solver(adU2,auU2,alU2,sx,bU2,iaU2,jaU2,neqU2,nadU2
     .               ,.true.,solvTolBcg,maxIt,solverU2,matrizU2,2)
          solvU2Time = get_time() - solvU2Time
          call getRes(u2,sx,num,neqU2,numel,.false.)
        endif
c .....................................................................
c
c ... atualizando o campo de velocidades estimadas
c ... interpolacao dos momentos         
c        if(intVel .eq. 3 ) then
c ... recontrucao do gradiente de u1
c          grSimpleTime = get_time() - grSimpleTime
c          call gform(u1,gradU1,fluxlU1,x,sedge,e,ls,ie,nelcon,pedge,ix
c     .              ,numel,ndm,nen,nshared,1,ndf,2,1,1)
c          grSimpleTime = get_time() - grSimpleTime
c .....................................................................
c
c ... recontrucao do gradiente de u2
c          grSimpleTime = get_time() - grSimpleTime
c          call gform(u2,gradU2,fluxlU2,x,sedge,e,ls,ie,nelcon,pedge,ix
c     .              ,numel,ndm,nen,nshared,1,ndf,2,2,1)
c          grSimpleTime = get_time() - grSimpleTime
c .....................................................................            
c
c ... calculo de divergente
c         call divergente(div,gradU1,gradU2,numel,ndf-1)
c .....................................................................
c
c ... u1
c          elmU1Time = get_time() - elmU1Time
c          call pform(adU1,auU1,alU1,iaU1,jaU1,num,un,u1,p,p,gradU1
c     .            ,gradU2,gradP,div,fluxlU1,iM,rO
c     .            ,wP,bU1,bU1,x,sedge,e,ddU,ie,nelcon,pedge,ix
c     .            ,numel,ndm,nen,nshared,1,ndf,dt,0.0d0
c     .            ,matrizU1,.true.,1,1,1,.false.,.false.,.true.
c     .            ,.false.,.false.,bs)
c          elmU1Time = get_time() - elmU1Time
c ..................................................................
c          
c ... u2
c          elmU2Time = get_time() - elmU2Time
c          call pform(adU2,auU2,alU2,iaU2,jaU2,num,un,u2,p,p,gradU2
c     .            ,gradU1,gradP,div,fluxlU2,iM,rO
c     .            ,wP,bU2,bU1,x,sedge,e,ddU,ie,nelcon,pedge,ix
c     .            ,numel,ndm,nen,nshared,1,ndf,dt,0.0d0
c     .            ,matrizU2,.true.,1,2,1,.false.,.false.,.true.
c     .            ,.false.,.false.,bs)
c          elmU2Time = get_time() - elmU2Time
c        endif
c .....................................................................
c
c ... atualizando o campo de velocidades estimadas
        simpleUpdateTime = get_time() - simpleUpdateTime 
        call updateVel(wP,u1,u2,numel,ndm)
        simpleUpdateTime = get_time() - simpleUpdateTime 
c .....................................................................
c
c      call uformnode(un,Pc,ddum,ddum,x,mdf,ix,md,nnode,numel,ndm
c   .                ,nen,1,1)
c ... equacao da correcao da pressao
        print*,'correcao de pressao'  
        elmPcTime = get_time() - elmPCTime
        call pform(adPc,auPc,alPc,iaPc,jaPc,num,un,pC,p,p,gradPc
     .            ,gradP,gradPc,div,fluxlPc,iM,rO
     .            ,wP,bPc,rCPc,x,sedge,e,ddU,eddyVisc
     .            ,ie,nelcon,pedge,ix
     .            ,numel,ndm,nen,nshared,1,ndf,dt,0.0d0
     .            ,matrizPc,unsymPc,1,2,2,.true.,.true.,.false.
     .            ,.false.,.true.,bs)
        elmPcTime = get_time() - elmPCTime
c .....................................................................
c 
c ... export coo format
        if(flagCoo) then
          call exportCoo(iaPc,jaPc,alPc   ,adPc    ,auPc,bPc,neqPc,nadPc
     .                  ,unsymPc  ,prename,noutCoo,itSimple ,matrizPc)
        endif
c .....................................................................
c 
c ...
        rPc = dsqrt(dot(rCPc,rCPc,neqPc))
        if( rPc .lt. ZERO) goto 100  
        if(itSimple .eq. kZero) rPc0 = rPc
        rU1 = dsqrt(dot(rCU1,rCU1,neqU1))
        rU2 = dsqrt(dot(rCU2,rCU2,neqU2))
        conv = 0
c .....................................................................       
c
c ... Solver ApC = pC
        call getRes(pC,sx,num,neqPc,numel,.true.)
        solvPcTime = get_time() - solvPCTime 
        call solver(adPc,auPc,alPc,sx,bPc,iaPc,jaPc,neqPc,nadPc
     .             ,unsymPc,solvTolPcg,maxIt,solverPc,matrizPc,2)
        solvPcTime = get_time() - solvPCTime
        call getRes(pC,sx,num,neqPc,numel,.false.)
c .....................................................................
c
c ... segunda equacao da correcao da pressao com malhas nao-ortognais
        if(skewnessCorrection) then
          print*,'skewness Correction: '
c ...     
          grSimpleTime = get_time() - grSimpleTime 
          call gform(pC,gradPc,fluxlPc,x,sedge,e,ls,ie,nelcon,pedge,ix
     .              ,numel,ndm,nen,nshared,1,ndf,2,2,2)
          grSimpleTime = get_time() - grSimpleTime 
c .....................................................................
c
c ...          
          elmPcTime = get_time() - elmPCTime
          call pform(adPc,auPc,alPc,iaPc,jaPc,num,un,pC,p,p,gradPc
     .              ,gradP,gradPc,div,fluxlPc,iM,rO
     .              ,wP,bPc,rCPc,x,sedge,e,ddU,eddyVisc
     .              ,ie,nelcon,pedge,ix
     .              ,numel,ndm,nen,nshared,1,ndf,dt,0.0d0
     .              ,matrizPc,unsymPc,4,2,2,.true.,.false.,.false.
     .              ,.false.,.false.,bs)
          elmPcTime = get_time() - elmPCTime
c .....................................................................
c
c ... Solver ApC = pC
          call getRes(pC1,sx,num,neqPc,numel,.true.)
          solvPcTime = get_time() - solvPCTime
          call solver(adPc,auPc,alPc,sx,bPc,iaPc,jaPc,neqPc,nadPc
     .              ,unsymPc,solvTolPcg,maxIt,solverPc,matrizPc,2)
          solvPcTime = get_time() - solvPCTime
          call getRes(pC1,sx,num,neqPc,numel,.false.)
c .....................................................................
c
c ...
          call asumb(pC,pC1,neqPc,pC)
c .....................................................................          
        endif
c .....................................................................
c
c ... recontrucao do gradiente da pressao de correcao
        grSimpleTime = get_time() - grSimpleTime 
        call gform(pC,gradPc,fluxlPc,x,sedge,e,ls,ie,nelcon,pedge,ix
     .            ,numel,ndm,nen,nshared,1,ndf,2,2,2)
        grSimpleTime = get_time() - grSimpleTime 
c .....................................................................
c
c ... atualizacao de u,v e p
        simpleUpdateTime = get_time() - simpleUpdateTime 
        call simpleUpdate(w,p,u1,u2,pC,gradPc,ddU,numel,ndf-1
     .                   ,underU,underPc,simpleR)
        simpleUpdateTime = get_time() - simpleUpdateTime
c .....................................................................
c
c ... equacao de transporte da energia
  100   continue
        if(sEnergy) then
c ... posprocessamento do campo de velocidade
c         posVelocityTime = get_time() - posVelocityTime               
c         call posProVelocityField(w,x,e,ie,nelcon,ix,wP,pedgeE
c    .                           ,numel,ndm,nen,nshared,5,3)
c         posVelocityTime = get_time() - posVelocityTime
c ......................................................................
c
c ...
          print*,'Transporte de energia'
c ... recontrucao do gradiente de en
          grSimpleTime = get_time() - grSimpleTime 
          call gform(en,gradE,fluxlE,x,sedgeE,e,ls,ie,nelcon,pedgeE,ix
     .              ,numel,ndm,nen,nshared,1,1,2,1,3)
          grSimpleTime = get_time() - grSimpleTime
c .....................................................................
c
c ...
          elmETime = get_time() - elmETime  
          call pform(adE,auE,alE,iaE,jaE,num,un,en,en,en,gradE,gradE
     .             ,gradE,div,fluxlE,iM,rO
     .             ,w,bE,rCe,x,sedgeE,e,ddU,eddyVisc
     .             ,ie,nelcon,pedgeE,ix
     .             ,numel,ndm,nen,nshared,1,1,dt,0.0d0
     .             ,matrizE,.true.,1,1,3,.true.,.true.,.false.
     .             ,.false.,.true.,bs)     
          elmETime = get_time() - elmETime
c ...  
          call asumb(bE,bE0,neqE,bE)
          call getRes(sx,bE0,num,neqE,numel,.false.)
          call asumb(rCe,sx,neqE,rCe)
c .....................................................................
c
c ...
          rE = dsqrt(dot(bE,bE,neqE))
          if( rE .lt. ZERO) goto 250
          rE  = dsqrt(dot(rCe,rCe,neqE)) 
          if(solvE) rE0 = rE
          solvE = .false.
c .....................................................................
c
c ... Solver Ae = bE
          call getRes(en,sx,num,neqE,numel,.true.)
          solvETime = get_time() - solvETime 
          call solver(adE,auE,alE,sx,bE,iaE,jaE,neqE,nadE
     .               ,.true.,solvTolBcg,maxIt,solverE,matrizE,2)
          solvETime = get_time() - solvETime
          call getRes(en,sx,num,neqE,numel,.false.)
c .....................................................................
c
c ... campo de massa especifica
          call EnthalpyForTemp(en,temp,e,w,numel,ndm,.true.,1)
          if(closed) then 
            call closedPresure(temp,x,e,ie,nelcon,ix,cPc,numel,ndm,nen
     .                        ,nshared,6,3)
          endif
          call  varMassEsp(ro,temp,numel,1,1,3)  
c ......................................................................
        endif
  250   continue
c ......................................................................
c
c ...
        if(ResAbs) then
          if( rPc .lt. solvSimplePc ) conv = conv + 1
          if(.NOT.solvU1) then
            if( rU1 .lt. solvSimpleU1 ) conv = conv + 1
          endif
          if(.NOT.solvU2) then
            if( rU2 .lt. solvSimpleU2 ) conv = conv + 1
          endif
          if(sEnergy) then
            if( rE  .lt. solvEnergy   ) conv = conv + 1
          endif 
c .....................................................................       
c
c ...
        else 
         if( rPc/rPc0   .lt. solvSimplePc 
     .               .or.
     .            rPc   .lt. ZERO) conv = conv + 1
         if(.NOT.solvU1) then
           if( rU1/rU10 .lt. solvSimpleU1 
     .               .or.
     .            rU1   .lt. ZERO) conv = conv + 1
         endif
         if(.NOT.solvU2) then
           if( rU2/rU20 .lt. solvSimpleU2 
     .               .or.
     .            rU2   .lt. ZERO) conv = conv + 1
         endif
         if(sEnergy) then
            if( rE/rE0   .lt. solvEnergy   ) conv = conv + 1
         endif
       endif
c .....................................................................       
c
c ...      
        if(sEnergy) then  
          if(conv .eq. 4) goto 300
        else
          if(conv .eq. 3) goto 300
        endif 
c .....................................................................       
c
c ...
        write(*,'(1x,a,i5)')    'It simple:                      ' 
     .                     ,itSimple
c .....................................................................       
c
c ...
        if(ResAbs) then
          write(*,'(1x,a,es20.8)')'residuo da conservacao de massa:'
     .                           ,rPc
          write(*,'(1x,a,es20.8)')'residuo do mometum x1          :'
     .                           ,rU1
          write(*,'(1x,a,es20.8)')'residuo do mometum x2          :'
     .                           ,rU2
          if(sEnergy) then
            write(*,'(1x,a,es20.8)')'residuo da trans de energia    :'
     .                          ,rE
          endif    
c .....................................................................       
c
c ...
        else
          write(*,'(1x,a,es20.8)')'residuo da conservacao de massa:'
     .                           ,rPc/rPc0
          write(*,'(1x,a,es20.8)')'residuo do mometum x1          :'
     .                           ,rU1/rU10
          write(*,'(1x,a,es20.8)')'residuo do mometum x2          :'
     .                           ,rU2/rU20
          if(sEnergy) then
            write(*,'(1x,a,es20.8)')'residuo da trans de energia    :'
     .                             ,rE/rE0
          endif 
        endif
c .....................................................................
        if(sEnergy) then
          write(noutSimple,'(i5,8es20.6)')itSimple
     .         ,rPc/rPc0,rPc,rU1/rU10,rU1,rU2/rU20,rU2,rE/rE0,rE
        else
           write(noutSimple,'(i5,6es20.6)')itSimple
     .                    ,rPc/rPc0,rPc,rU1/rU10,rU1,rU2/rU20,rU2
        endif
c .....................................................................
c
c ...
        simpleItTime = get_time() - simpleItTime
        write(*,'(1x,a,12x,a,1x,f19.5)')
     .       'Simple CPU time (s)',':',simpleItTime
c .....................................................................
      enddo simpleLoop
c .....................................................................
c
c ...
  300 continue
c ...
c     if(sEnergy) then  
c       if(conv .ne. 4) then
c         print*,"Simple nao convergiu!!"
c         stop
c       endif
c     else
c       if(conv .ne. 3) then
c         print*,"Simple nao convergiu!!"
c         stop
c       endif
c     endif 
      write(*,'(1x,a,i5)')    'it simple',itSimple
      write(*,'(1x,a,es20.8)')'residuo da conservacao de massa:'
     .                       ,rPc
      write(*,'(1x,a,es20.8)')'residuo do mometum x1          :'
     .                       ,rU1
      write(*,'(1x,a,es20.8)')'residuo do mometum x2          :'
     .                       ,rU2    
      if(sEnergy) then
        write(*,'(1x,a,es20.8)')'residuo da trans de energia    :'
     .                         ,rE
      endif    
c .....................................................................
c
c ... calculo de parametro por celula
      call cellParameter(mP   ,rO     ,w      ,x  ,eddyVisc
     .                  ,gradU1  ,gradU2
     .                  ,e       ,sedge  ,pedge  ,ie ,nelcon
     .                  ,ix      ,numel  ,ndm    ,nen,nshared 
     .                  ,ndf     ,dt     ,      7,  1)
      call calParameter(mP     ,div       ,numel         ,cfl ,reynalds
     .                 ,kEnergy,disfiltro,transSubGrid
     .                 ,prandtl,grashof  ,vol   ,Massa   ,fluxoM
     .                 ,dt     ,tDinamico)
c .....................................................................
c
c ...
      if(sEnergy) call massUpdate(ro,numel)
c .....................................................................
      return
      end
c *********************************************************************
c * SIMPLEUPDATE : atualizacao das variaveis do metodo simple         *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * w       - campo de velocidade desatualizado                       *
c * pre     - campo de pressao desatualizado                          *
c * u1      - campo de velocidade u estimado                          *
c * u3      - campo de velocidade w estimado                          *
c * preC    - pressao de correcao                                     *
c * gradP   - gradiente da pressao de correcao                        *
c * d       -                                                         *
c * numel   - numero de elementos                                     *
c * underU  - fator de sobrerrelaxamento para as velocidades          *
c * underP  - fator de sobrerrelaxamento para a pressao               *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * w       -vetor atualizado                                         *
c * pre     -vetor atualizado                                         *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine simpleUpdate(w,pre,u1,u2,preC,gradPc,d,numel,ndm
     .                       ,underU,underPc,simpleR)
      implicit none
      real(8) w(ndm,*),u1(*),u2(*),pre(*),preC(*),d(ndm,*)
      real(8) gradPc(ndm,*)
      real(8) underU,underPc
      integer numel,ndm
      integer i
      logical simpleR
c ...
      if(simpleR) then
        do i = 1, numel
          w(1,i) = u1(i)  - d(1,i)*gradPc(1,i)
          w(2,i) = u2(i)  - d(2,i)*gradPc(2,i)
        enddo
c .....................................................................
c
c ...
      else
        do i = 1, numel
          w(1,i) = u1(i)   - d(1,i)*gradPc(1,i)
          w(2,i) = u2(i)   - d(2,i)*gradPc(2,i)
          pre(i) = pre(i)  + underPc*preC(i)
        enddo
      endif
c .....................................................................
      return
      end
c .....................................................................
c *********************************************************************
c
c *********************************************************************
c * updateVEl:                                          X,Y,Z         *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine updateVel(w,u1,u2,numel,ndm)
      implicit none
      real(8) w(ndm,*),u1(*),u2(*)
      integer ndm,numel
      integer i
      do i = 1, numel
        w(1,i) = u1(i)
        w(2,i) = u2(i)
      enddo
      return
      end
c *********************************************************************
c
c *********************************************************************
c * PSEUDOVEl:                                                        *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine pseudoVel(w,iM,numel,ndm)
      implicit none
      real(8) w(ndm,*),iM(2*ndm,numel)
      integer ndm,numel
      integer i
      do i = 1, numel
        w(1,i) = iM(1,i)/iM(3,i)
        w(2,i) = iM(2,i)/iM(4,i)
      enddo
      return
      end
c *********************************************************************
c
c *********************************************************************
c * GUESS:   :                                                        *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c *********************************************************************    
      subroutine guess(u1,u2,w,numel,n)
      implicit none
      real(8) u1(*),u2(*),w(n,*)
      integer numel,n
      integer i
      do i = 1, numel
        u1(i) = w(1,i)
        u2(i) = w(2,i)
      enddo
      return
      end
c *********************************************************************
c
c *********************************************************************
c * GUESS:   :                                                        *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c *********************************************************************      
      subroutine pressureGuess(p,pedge,sedge,numel,nshared,ndf)
      implicit none
      real(8) p(*),sedge(ndf,nshared+1,numel)
      integer numel,pedge(nshared+1,numel)
      integer i,nshared,ndf,idCell
      idCell = nshared + 1
      do i = 1, numel
        if(pedge(idCell,i) .eq. -1) p(i) = sedge(ndf,idCell,i)
      enddo
      return
      end
c *********************************************************************
c
c *********************************************************************
c * MASSUPDATE:                                                       *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c *********************************************************************    
      subroutine massUpdate(ro,numel)
      implicit none
      real(8) ro(3,numel)
      integer numel,i
      do i = 1, numel
        ro(1,i) = ro(2,i)
        ro(2,i) = ro(3,i)
      enddo
      return
      end
c *********************************************************************
c
c *********************************************************************
c * DIVERGENTE: Calculo do divergente da velocidade                   *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * div    - nao definido                                             *
c * gradU1 - gradiente do u1                                          *
c * gradU2 - gradiente do u2                                          *
c * numel  - numero de elementos                                      *
c * n      - dimensao                                                 *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * div    - divergente                                               *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine divergente(div,gradU1,gradU2,numel,n)
      implicit none
      real(8) div(numel),gradU1(n,numel),gradU2(n,numel)
      integer numel,i,n
      do i = 1, numel
        div(i) = gradU1(1,i) + gradU2(2,i)
      enddo
      return
      end
c *********************************************************************
c
c *********************************************************************
c * calParameter : calculo de parametros importantes no dominio       *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * mP    - parametros por celula                                     *
c *            1 - cfl                                                *
c *            2 -Reynalds                                            *
c *            3 - volume                                             *
c *            4 - velocidade Vm                                      *
c *            5 - prandlt numer local                                *
c *            6 - grahosf local                                      *
c *            7 - massa do volume                                    * 
c *            8 - fluxo de massa entrando no dominio                 * 
c *            9 - fluxo de massa saindo do dominio                   * 
c *           10 - energia cinetica total                             * 
c *           11 - dissipacao de energia resolvivel                   * 
c *           12 - transferencia de energia resolvivel para o sub grid* 
c * numel - numero de elementos                                       *
c * massa - indefinido                                                *
c * fM    - indefinido                                                *
c * vol   - indefinido                                                *
c * dt    - intervalo de tempo                                        *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * cfl          - numero de courant                                  *
c * re           - numero de reynalds                                 *
c * kEnrgy       - energia cinetica total                             *
c * disfiltro    - dissipaco de energia resolvivel pelas tensoes      *
c * viscosas                                                          *   
c * transSubGrid - tranferencia de enrgia para o subgrid              *
c * pr           - Prandtl number                                     *
c * gr           - Grashof number                                     *        
c * vol          - volume total do dominio                            *
c * massa        - massa total do dominio                             *
c * fm           - fluxo de massa atraves do domino aberto            *
c * dt           - intervalo de tempo                                 *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine calParameter(mP       ,div          ,numel,cfl,re
     .                       ,kEnergy  ,disfiltro    ,transSubGrid
     .                       ,pr       ,gr           ,vol 
     .                       ,massa   ,fM            ,dt
     .                       ,deltaDinamico)
      implicit none
      real(8) mP(20,numel),div(*),vol,dt,dtm,dtm1
      real(8) cfl,re,kEnergy,disfiltro,transSubGrid
      real(8) pr,gr
      real(8) massa,fM,massInput,massOut
      logical deltaDinamico
      integer numel
      integer i
c ...
      cfl          =  mP(1,1)
      re           =  mP(2,1)
      vol          =  mP(3,1)
      pr           =  mP(5,1)
      gr           =  mP(6,1)
      massa        =  mP(7,1)
      massInput    =  mP(8,1)
      massOut      =  mP(9,1)
      kEnergy      = mP(10,1)
      disfiltro    = mP(11,1)
      transSubGrid = mP(12,1)
      dtm  = min(dt,0.9d0*(dsqrt(mP(3,1))/mP(4,1)))
      dtm1  = min(dt,dabs(1/div(1)))
c .....................................................................
c
c ...      
      do i = 2, numel
c ... cfl maximo da malha
        cfl        = max(mP(1,i),cfl)
c ... reynalds maximo da malha
        re         = max(mP(2,i),re)
c ... volume da dominio
        vol        = vol + mP(3,i)
c ... prandlt maximo da malha
        pr         = max(mP(5,i),pr)
c ... grasholf maximo da malha
        gr         = max(mP(6,i),gr)
c ... massa total
        massa      = massa + mP(7,i)
c ... massa total input
        massInput  = massInput + mP(8,i)
c ... massa total output
        massOut    = massOut   + mP(9,i)
c ...            
        kEnergy    = kEnergy   + mP(10,i)
c ...            
        disfiltro    = disfiltro    + mP(11,i)
c ...
        transSubGrid = transSubGrid + mP(12,i)
c ... deltaT que satisfaca o cfl
        dtm    = min(dtm,0.9d0*(dsqrt(mP(3,i))/mP(4,i)))
c ... continuidade
        dtm1   = min(dtm1,dabs(1/div(i)))
      enddo
c .....................................................................
c
c ...
      dtm = min(dtm,dtm1)
      if(dt .gt. dtm .and. deltaDinamico) then
        dt = dtm
        write(*,'(1x,a,es20.8)')'Delta t madificado para',dt
      endif
c .....................................................................
c
c ...
      kEnergy      = kEnergy/(massa)
      disfiltro    = disfiltro/vol
      transSubGrid = transSubGrid/vol
c ...
      fM = massInput + massOut
c .....................................................................
      return
      end
c *********************************************************************
c
c *********************************************************************
c * VARMASSESP : variacao da mass especifica                          *
c ------------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ro   F            - massa especifica (n+1,n,n-1)                  *
c * temp              - temperatura                                   *
c * typeF             - 1 - Ar; 2 - Agua                              *
c * typeT             - 1 - °C; 0 - Kelvin                            *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * ro(3,i)  - massa especifica (n+1)                                 *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine  varMassEsp(ro,temp,numel,typeF,typeT,pos)
      implicit none
      real(8) ro(*),temp(*)
      integer numel
      integer typeT,typeF,pos
c ... gas ideal Ar
      if(typeF .eq. 1) then
        call idealGas(ro,temp,numel,typeT,pos)
c ... agua
      else if(typeF .eq. 2) then
        call water(ro,temp,numel,typeT)
      endif
      return
      end
c *********************************************************************
c
c *********************************************************************
c * CHECKD : checa caracteriscas do domino                            *
c ------------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * pedgeF            - tipo de condicao de contorno                  *
c * numel             - numero de elementos                           *
c * singularPressure  - nao definido                                  *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * singularPressure  - (true|false)                                  *
c * cloded            - (true|false)                                  *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine checkD(pedgeF,numel,nshared,singularPressure
     .                 ,closed)
      implicit none
      integer pedgeF(nshared+1,*)
      integer numel,nshared
      logical singularPressure,closed,checkS,checkC
      integer i,j
c ...      
      checkS           = .false.
      singularPressure = .false.
      checkC           = .true. 
      do i = 1, numel
        if( pedgeF(nshared+1,i) .eq. -1) checkS = .true.
        do j = 1, nshared
           if( pedgeF(j,i) .eq. 3) checkS = .true.
           if( pedgeF(j,i) .eq. 2 .or. pedgeF(j,i) .eq. 3) then 
             checkC = .false.
          endif
        enddo
      enddo
      if(.not.checkS) singularPressure = .true.
      if(checkC)      closed           = .true.
c ......................................................................      
      return
      end
c *********************************************************************
c
c *********************************************************************
c * SAVESIMPLE : salva o estado das principais variaveis no metodo    *
c * simple                                                            *
c ------------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * w                 - campo de velocidade no tempo (n)              *
c * w0                - campo de velocidade no tempo (n-1)            *
c * p                 - campo de pressao                              *
c * en                - campo de temperatura (n-1)                    *
c * en0               - campo de temperatura (n-2)                    *
c * ro                - campo de massa especifica (n-1,n-1,n-2)       *
c * numel             - numero de elementos                           *
c * ndm               - numero de dimensoes                           *
c * istep             - passo de tempo                                *
c * t                 - tempo de simulacao                            *
c * load              - true|false                                    *
c * sEnergy           - equacao de conservacao de energia             *
c * nout              -                                               *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c *********************************************************************      
      subroutine saveSimple(w,w0,p,en,en0,ro,nnode,numel,ndm
     .                     ,istep,t,fileout,load,sEnergy,nout)
      implicit none
      integer nnode,numel,ndm,istep
      integer nout
      integer i,j
      real(8) w(ndm,*),w0(ndm,*),p(*),en(*),en0(*)
      real(8) ro(3,*),t
      character*80 fileout
      character string*80
      logical load,sEnergy
c ...
      if(load) then
        open(unit=nout,file=fileout,status='old')
      else  
        open(unit=nout,file=fileout)
      endif
c .....................................................................
c
      if(load) then
        read(nout,'(80a)') string
        read(nout,'(i9)') istep
        read(nout,'(80a)') string
        read(nout,'(f20.8)') t 
      else  
        write(nout,'(a)')'passo de tempo'
        write(nout,'(i9)') istep
        write(nout,'(a)')'tempo '
        write(nout,'(f20.8)') t     
      endif
c ... velocidade      
      if(load) then
        read(nout,'(80a)') string
        do i=1, numel
          read(nout,'(3f20.8)') (w(j,i),j = 1,ndm)
        enddo
      else
        write(nout,'(a)')'Velocidade (n)'
        do i=1, numel
          write(nout,'(3f20.8)') (w(j,i),j = 1,ndm)
        enddo
      endif
c .....................................................................
c
c ... velocidade
      if(load) then
        read(nout,'(80a)') string
        do i=1, numel
          read(nout,'(3f20.8)') (w0(j,i),j = 1,ndm)
        enddo
      else
        write(nout,'(a)')'Velocidade (n-1)' 
        do i=1, numel
          write(nout,'(3f20.8)') (w0(j,i),j = 1,ndm)
        enddo
      endif      
c .....................................................................
c
c ... pressao
      if(load) then
        read(nout,'(80a)') string
        do i=1, numel
          read(nout,'(f20.8)') p(i)
        enddo
      else    
        write(nout,'(a)')'pressao (n)'       
        do i=1, numel
          write(nout,'(f20.8)') p(i)
        enddo
      endif
c .....................................................................
c
c ... Temperatura
      if(sEnergy) then         
        if(load) then
          read(nout,'(80a)') string
          do i=1, numel
            read(nout,'(f20.8)') en(i)
          enddo
        else    
          write(nout,'(a)')'en (n)'       
          do i=1, numel
            write(nout,'(f20.8)') en(i)
          enddo
        endif
c .....................................................................
c
c ... Temperatura    
        if(load) then
          read(nout,'(80a)') string
          do i=1, numel
            read(nout,'(f20.8)') en0(i)
          enddo
        else 
          write(nout,'(a)')'en (n-1)'        
          do i=1, numel
            write(nout,'(f20.8)') en0(i)
          enddo
        endif
      endif 
c .....................................................................
c
c ... massa especefica    
      if(load) then
        read(nout,'(80a)') string
        do i=1, numel
          read(nout,'(3f20.8)') ro(1,i),ro(2,i),ro(3,i)
        enddo
      else 
        write(nout,'(a)')'ro      '        
        do i=1, numel
          write(nout,'(3f20.8)') ro(1,i),ro(2,i),ro(3,i)
        enddo
      endif
c .....................................................................
c
c ...
      close(nout)
      return
      end
c ********************************************************************* 
