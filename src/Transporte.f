c *********************************************************************
c * TRANSPORTE: Monta a matriz de coeficente e calcula do residuo R   *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine transporte(x      ,ix     ,e
     .                     ,ie     ,nelcon ,pedgeT1
     .                     ,sedgeT1,w      ,num
     .                     ,ls     ,sx
     .                     ,gradT1 ,fluxlT1,rCt1
     .                     ,adT1   ,auT1   ,alT1
     .                     ,bT1    ,bT10
     .                     ,iaT1   ,jaT1   ,t1
     .                     ,t10    ,un     ,rO
     .                     ,numel  ,nnode  ,ndm
     .                     ,nen    ,nshared,ndfT1
     .                     ,dt     ,t      ,matrizT1
     .                     ,neqT1  ,nadT1  ,solverT1
     .                     ,solvTol,maxIt  ,istep
     .                     ,unsymT1,solvT1 ,maxItT1
     .                     ,bs)
      implicit none
      include 'time.fi'
      integer nnode,numel,ndm,nen,nshared,ndfT1
      real*8 x(ndm,*),e(*),sedgeT1(*),rO(3,*),w(ndm,*),ls(*),un(*)
      integer ix(nen+1,*),ie(*),pedgeT1(nshared+1,*)
      integer num(*),nelcon(nshared,*)
c ... gradientes
      real*8 gradT1(ndm,*),fluxlT1(*)
c ... sistema
      real*8  adT1(*),auT1(*),alT1(*),bT1(*),bT10(*)
      integer iaT1(*),jaT1(*)
      real*8  t1(*),t10(*)
      real*8  rCt1(*),sx(*)
      integer neqT1,nadT1,maxIt,matrizT1,solverT1
      real*8 solvTol,dot,rT1,rT10
c ...
      real*8  t,dt,solvT1
      integer istep,iloop,maxItT1
      logical bs,unsymT1
c ...
      elmTime   = get_time() - elmTime
      call pform(adT1,auT1,alT1,iaT1,jaT1
     .          ,num,un,t1,t10,t10
     .          ,gradT1,gradT1,gradT1
     .          ,fluxlT1,fluxlT1,fluxlT1,rO
     .          ,w,bT10,rCt1,x
     .          ,sedgeT1,e,fluxlT1,fluxlT1,ie,nelcon
     .          ,pedgeT1,ix
     .          ,numel,ndm,nen,nen,ndfT1,ndfT1,dt
     .          ,0.0d0,matrizT1,unsymT1,4,1,4,.true.,.false.
     .          ,.false.,bs)
      elmTime   = get_time() - elmTime
c ...       
      call aequalb(t10,t1,neqT1)
c ......................................................................
c
c ... loop nao linear
      do iloop = 1, maxItT1
        print*,'Transporte da propriedade T1'
c ... recontrucao do gradiente de en
        grTime = get_time() - grTime 
        call gform(t1,gradT1,fluxlT1,x,sedgeT1,e,ls,ie,nelcon
     .            ,pedgeT1,ix,numel,ndm,nen,nen,ndfT1,ndfT1,2,1,4)
        grTime = get_time() - grTime
c ......................................................................
c
c ... momtagem do sistema At=F
        elmTime   = get_time() - elmTime   
        call pform(adT1,auT1,alT1,iaT1,jaT1
     .            ,num,un,t1,t1,t1
     .            ,gradT1,gradT1,gradT1
     .            ,fluxlT1,fluxlT1,fluxlT1,rO
     .            ,w,bT1,rCt1,x
     .            ,sedgeT1,e,fluxlT1,fluxlT1,ie,nelcon
     .            ,pedgeT1,ix
     .            ,numel,ndm,nen,nen,ndfT1,ndfT1,dt
     .           ,0.0d0,matrizT1,unsymT1,1,1,4,.true.,.true.
     .           ,.false.,bs)
        elmTime   = get_time() - elmTime
c ...          
        call asumb(bT1,bT10,neqT1,bT1)
        call getRes(sx,bT10,num,neqT1,numel,.false.)
        call asumb(rCt1,bT10,neqT1,rCt1)
c .....................................................................
c
c ...
        write(*,'(1x,a,i5)')    'It T1 transport:                     ' 
     .                     ,iloop
c .....................................................................
c
c ...
        rT1  = dsqrt(dot(rCt1,rCt1,neqT1)) 
        if(iloop .eq. 1) rT10 = rT1
        if(rT1/rT10 .lt. solvT1) return
c .....................................................................
c
c ...
        write(*,'(1x,a,es20.8)')'residuo:                             '
     .                          ,rT1/rT10    
c .....................................................................
c
c ... Solver At1 = bT1
        solvTime = get_time() - solvTime 
        call getRes(t1,sx,num,neqT1,numel,.true.)
        call solver(adT1,auT1,alT1,sx
     .             ,bT1,iaT1,jaT1,neqT1,nadT1
     .             ,.true.,solvTol,maxIt,solverT1,matrizT1,2)
        call getRes(t1,sx,num,neqT1,numel,.false.)
        solvTime = get_time() - solvTime
c .....................................................................
      enddo
c ...
      return
      end
c *********************************************************************       
