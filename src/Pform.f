c *********************************************************************
c * PFORM: Monta a matriz de coeficente e calcula do residuo R        *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ad     - nao definido                                             *
c * au     - nao definido                                             *
c * al     - nao definido                                             *
c * ia     - ponteiro do para o csr                                   *
c * ja     - estrutura do csr                                         *
c * num    - renumeracaos do elementos para as equacoes               *
c *  u     - variavel na celula                                       *
c *  u0    - variavel na celula                                       *
c *  u1    - variavel na celula                                       *
c *  grad  - gradiente da celula                                      *
c *  grad1 - gradiente da celula                                      *
c *  grad2 - gradiente da celula                                      *
c *  div   - diveregente da velocidade                                *
c *  fluxl - limitador de fluxo da celula                             *
c *     iM - nao definido                                             *
c *     ro - massa especifica da celula(pi+1,pt)                      *   
c *  w     - campo de velocidade por celula                           *
c *  f     - nao definido                                             *
c * resCell- nao definido                                             *
c *  x     - coordenadas nodais                                       *
c * sedge  - condicoes de contorno nas arestas                        *
c *  e     - propriedade                                              *
c * ddU    - campo d o algoritimo simple                              *
c * eddyVisc - viscosidade turbulenta                                 *
c *  ie    - tipo do elemento do material                             *
c * nelcon - vizinho da celula                                        *
c * pedge  - tipo de condicao de contorno                             *
c * ix     - conectividade dos elementos                              *
c * numel  - numero de celula                                         *
c * ndm    - numero de dimensoes                                      *
c * nen    - numero de nos por celula                                 *
c * nshared- numero de celula compartilhados pela celula central      *
c * ndf    - graus de liberdade da propriedade transportada           *
c * ndfF   - graus de liberdade do fluido                             *
c * dt     - passo de tempo                                           *
c * alpha  -                                                          *
c * struc  - estrutura da matriz A                                    *
c *        1 -    CSR ( ad - diag prin; a - superio e inferior)       *
c *        2 -    CSR (  a - diag prin, superior e inferior           *
c *        3 -    CSRC( ad - diag prin, au - superior, al - inferior) *
c *        4 - ELLPACK(ad - diag prin, al - inferior e superior)      *
c * unsym  - matriz nao simetrica                                     *
c * iws    - instrucao de montagem                                    *
c * iws1   - 1 - x, 2 - y, 3 - z                                      *
c * forces - mantagem da vetor de forcas                              *
c * matrix - montagem da matrix                                       *
c * momentum - motagem da matriz de momentunm de navier stockes       *
c * calD     - calculo do campo D                                     *
c * resid    - motagem do residuo                                     *
c * bs       - Euller Backward de segunda ordem                       *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * ad     - valores da diagonal da matriz de coeficientes            *
c * al,au  - valores fora da diagonal de coeficientes                 *
c *  f     - vetor de forcao ou residuo da equacao                    *
c * resCell- residuo da celula                                        *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine pform(ad,au,al,ia,ja,num,un,u,u0,u1,grad,grad1,grad2
     .                ,div,fluxl,iM,rO,w,f,resCell,x,sedge,e,ddU
     .                ,eddyVisc 
     .                ,ie,nelcon
     .                ,pedge,ix
     .                ,numel,ndm,nen,nshared,ndf,nst,dt
     .                ,alpha,struc,unsym,iws,iws1,lib,forces,matriz
     .                ,momentum,calD,resid,bs)
      implicit none
      include 'openmp.fi'
      integer numel,nen,ndm,ndf,nst,nshared,dum
      integer struc,iws,iws1,lib
      logical unsym,forces,matriz,momentum,resid,calD,bs
      integer ia(*),ja(*),num(*),ix(nshared+1,numel),ie(*)
      integer nelcon(nshared,*),pedge(nshared+1,numel)
      real*8 ad(*),au(*),al(*),u(ndf,numel),u0(ndf,numel),u1(ndf,numel)
      real*8 grad(ndm,numel),grad1(ndm,numel),grad2(ndm,numel)
      real*8 div(numel),eddyVisc(numel)
      real*8 ddU(nst-1,numel),iM(2*(nst-1),numel),uN(ndf,*)
      real*8 f(ndf,*),sedge(nst,nshared+1,numel),w(ndm,*),fluxl(numel)
      real*8 x(ndm,*),e(10,*),dt,alpha,resCell(*),rO(3,*)
c ... variaveis locais por elemento      
      real*8  sp,la(7),lsedge(28)
      real*8  lx(168),le(10,7)
      real*8  p(3),lUn(24)
      real*8  ldd(21)
      real*8  lfluxl(7)
      real*8  liM(42)
      real*8  lRo(3,7),ldiv(7),lEddyVisc(7)
      real*8  lu(28),lu0(28),lu1(28)
      real*8  lw(21),lgrad(21),lgrad1(21),lgrad2(21)
      integer ld(7),lpedge(7),lviz(6)
c ... variaveis locais
      integer nel,i,ii,j,jj,k,kk
      integer viznel,no,ma,type,idCell       
c .....................................................................
c
      idCell = nshared + 1
c ... openmp
      if(openmpCell) then
c ... loop nas celulas
c$omp parallel private(nel,i,ii,j,jj,k,vizNel,no,ma,type)
c$omp.private(sp,la,lsedge,lx,le,p,lUn,ldd,lfluxl,liM,lRo,ldiv)
c$omp.private(lEddyVisc,lu,lu0,lu1,lw,lgrad,lgrad1,lgrad2)
c$omp.private(ld,lpedge,lviz)
c$omp.shared(numel,nen,ndm,ndf,nst,nshared,dum,struc,iws,iws1,lib)
c$omp.shared(unsym,forces,matriz,momentum,resid,calD,bs)
c$omp.shared(ia,ja,num,ix,ie,nelcon,pedge)
c$omp.shared(ad,au,al,u,u0,u1,grad,grad1,grad2,div,ddU,eddyVisc)
c$omp.shared(iM,f,sedge,w)
c$omp.shared(fluxl,x,e,dt,alpha,resCell,rO,un) num_threads(nThreadsCell)
c$omp do
        do nel = 1, numel
          sp   = 0.0d0
          ld(1:idCell) = 0
c .....................................................................        
c
c ... 
          do i = 1, nshared + 1 
            do k = 1, nst
              ii      = (i-1)*nst + k  
              lsedge(ii) = sedge(k,i,nel)
            enddo
            lpedge(i)     = pedge(i,nel)
          enddo
c .....................................................................
c
c ... loop da celula central
          lEddyVisc(idCell) = eddyVisc(nel) 
          ld(idCell)        = num(nel)
          lfluxl(idCell)    = fluxl(nel)
          ldiv(idCell)      =   div(nel)
          lrO(1,idCell)     = rO(1,nel)
          lrO(2,idCell)     = rO(2,nel)
          lrO(3,idCell)     = rO(3,nel)
c ...
          kk = nst -1
          do k = 1, kk 
            ii = nshared*kk + k  
            ldd(ii)  = ddU(k,nel)
          enddo
c ..
          kk = 2*(nst -1)
          do k =1 , kk
            ii       = nshared*kk + k  
            liM(ii)  = iM(k,nel)
          enddo
c ...
          do k = 1, ndf
            ii      = nshared*ndf + k 
            lu(ii)  =  u(k,nel)
            lu0(ii) = u0(k,nel)
            lu1(ii) = u1(k,nel)  
          enddo
c ...
          do k = 1, ndm
            ii         = nshared*ndm + k
            lw(ii)     =     w(k,nel)   
            lgrad(ii)  =  grad(k,nel)
            lgrad1(ii) = grad1(k,nel)
            lgrad2(ii) = grad2(k,nel)
          enddo
c ...
 
          do i = 1, nen
            no        = ix(i,nel)
            ii = nshared*nen*ndm 
            do k = 1, ndm
              jj     = ii + (i-1)*ndm + k
              lx(jj) = x(k,no)
            enddo
          enddo
c .....................................................................
c
c ... loop na aresta        
          do j = 1 , nshared
            viznel            = nelcon(j,nel)
            lviz(j)           = viznel 
            if( viznel .gt. 0) then
              lEddyVisc(j)      = eddyVisc(vizNel) 
              ldiv(j)           =  div(vizNel)
              lrO(1,j)          = rO(1,vizNel)
              lrO(2,j)          = rO(2,vizNel)
              lrO(3,j)          = rO(3,vizNel)
              kk = nst -1
              do k = 1, kk 
                ii = (j-1)*kk + k  
                ldd(ii)  = ddU(k,vizNel)
              enddo
              kk = 2*(nst -1)
              do k =1 , kk
                ii      = (j-1)*kk + k 
                liM(ii) = iM(k,vizNel)
              enddo
              ld(j)     = num(viznel)
              lfluxl(j) = fluxl(viznel)
              ma = ix(nen+1,viznel)
              do k = 1, ndf
                ii      = (j-1)*ndf + k 
                lu(ii)  =  u(k,viznel)
                lu0(ii) = u0(k,viznel)
                lu1(ii) = u1(k,vizNel)
              enddo
              do kk = 1, 10
                le(kk,j) = e(kk,ma)
              enddo  
              do k = 1, ndm
                ii          = (j-1)*ndm + k
                lw(ii)      =     w(k,vizNel)
                lgrad(ii)   =  grad(k,vizNel)
                lgrad1(ii)  = grad1(k,vizNel)
                lgrad2(ii)  = grad2(k,vizNel)
              enddo
c ... loop nos vertices do elemento vizinho          
              do i = 1, nen
                no = ix(i,vizNel)
                ii = (j-1)*nen*ndm
                do k=1, ndf
                  jj  = (i-1)*ndf + k 
                  lUn(jj) = un(k,no)
                enddo
                do k = 1, ndm
                  jj     = ii + (i-1)*ndm + k
                  lx(jj) = x(k,no)
                enddo
              enddo  
            endif  
          enddo
          ma            = ix(nen+1,nel)
          type = ie(ma)
          do i = 1, 10
            le(i,idCell) = e(i,ma)
          enddo
c ...
          call celllib(la    ,lx     ,lun   ,lu         ,lu0 
     .                ,lu1   ,lRo    ,lgrad ,lgrad1     ,lgrad2
     .                ,ldiv  ,liM    ,lfluxl,lw         ,ldd  
     .                ,le    ,sp     ,p     ,lsedge     ,lpedge
     .                ,lviz  ,dum    ,nen   ,nshared    ,ndm  
     .                ,type  ,iws    ,iws1  ,lib        ,nel 
     .                ,dt    ,alpha  ,dum   ,lEddyVisc  ,bs)
c .....................................................................
c
c ... calculo do residuo
          if(resid) resCell(nel)   = sp
c ... interpolacao pelo coeficiente dos momentos          
          if(momentum) then
            kk             = 2*(nst -1)
            ii             = nshared*kk  
            iM(iws1,nel)   = liM(ii+iws1)
            iM(iws1+2,nel) = liM(ii+iws1+2)
          endif  
c ... campo D            
          if(calD) then  
            kk             = nst -1
            ii             = nshared*kk + iws1  
            ddU(iws1,nel)  = ldd(ii)
          endif          
c .....................................................................
c
c ...
          if(matriz .or. forces) then
            call assbly(ad,au,al,f,ia,ja,la,p,ld,nshared,struc,unsym
     .                 ,forces,matriz)
          endif
c .....................................................................
      enddo
c .....................................................................
c$omp end parallel
c .....................................................................
c
c ... sequencial      
      else
        do nel = 1, numel
          sp   = 0.0d0
          ld(1:idCell) = 0
c .....................................................................        
c
c ... 
          do i = 1, nshared + 1 
            do k = 1, nst
              ii      = (i-1)*nst + k  
              lsedge(ii) = sedge(k,i,nel)
            enddo
            lpedge(i)     = pedge(i,nel)
          enddo
c .....................................................................
c
c ... loop da celula central
          lEddyVisc(idCell) = eddyVisc(nel) 
          ld(idCell)        = num(nel)
          lfluxl(idCell)    = fluxl(nel)
          ldiv(idCell)      =   div(nel)
          lrO(1,idCell)     = rO(1,nel)
          lrO(2,idCell)     = rO(2,nel)
          lrO(3,idCell)     = rO(3,nel)
c ...
          kk = nst -1
          do k = 1, kk 
            ii = nshared*kk + k  
            ldd(ii)  = ddU(k,nel)
          enddo
c ..
          kk = 2*(nst -1)
          do k =1 , kk
            ii       = nshared*kk + k  
            liM(ii)  = iM(k,nel)
          enddo
c ...
          do k = 1, ndf
            ii      = nshared*ndf + k 
            lu(ii)  =  u(k,nel)
            lu0(ii) = u0(k,nel)
            lu1(ii) = u1(k,nel)  
          enddo
c ...
          do k = 1, ndm
            ii         = nshared*ndm + k
            lw(ii)     =     w(k,nel)   
            lgrad(ii)  =  grad(k,nel)
            lgrad1(ii) = grad1(k,nel)
            lgrad2(ii) = grad2(k,nel)
          enddo
c ...
 
          do i = 1, nen
            no        = ix(i,nel)
            ii = nshared*nen*ndm 
            do k = 1, ndm
              jj     = ii + (i-1)*ndm + k
              lx(jj) = x(k,no)
            enddo
          enddo
c .....................................................................
c
c ... loop na aresta        
          do j = 1 , nshared
            viznel            = nelcon(j,nel)
            lviz(j)           = viznel 
            if( viznel .gt. 0) then
              lEddyVisc(j)      = eddyVisc(vizNel) 
              ldiv(j)           =  div(vizNel)
              lrO(1,j)          = rO(1,vizNel)
              lrO(2,j)          = rO(2,vizNel)
              lrO(3,j)          = rO(3,vizNel)
              kk = nst -1
              do k = 1, kk 
                ii = (j-1)*kk + k  
                ldd(ii)  = ddU(k,vizNel)
              enddo
              kk = 2*(nst -1)
              do k =1 , kk
                ii      = (j-1)*kk + k 
                liM(ii) = iM(k,vizNel)
              enddo
              ld(j)     = num(viznel)
              lfluxl(j) = fluxl(viznel)
              ma = ix(nen+1,viznel)
              do k = 1, ndf
                ii      = (j-1)*ndf + k 
                lu(ii)  =  u(k,viznel)
                lu0(ii) = u0(k,viznel)
                lu1(ii) = u1(k,vizNel)
              enddo
              do kk = 1, 10
                le(kk,j) = e(kk,ma)
              enddo  
              do k = 1, ndm
                ii          = (j-1)*ndm + k
                lw(ii)      =     w(k,vizNel)
                lgrad(ii)   =  grad(k,vizNel)
                lgrad1(ii)  = grad1(k,vizNel)
                lgrad2(ii)  = grad2(k,vizNel)
              enddo
c ... loop nos vertices do elemento vizinho          
              do i = 1, nen
                no = ix(i,vizNel)
                ii = (j-1)*nen*ndm
                do k=1, ndf
                  jj  = (i-1)*ndf + k 
                  lUn(jj) = un(k,no)
                enddo
                do k = 1, ndm
                  jj     = ii + (i-1)*ndm + k
                  lx(jj) = x(k,no)
                enddo
              enddo  
            endif  
          enddo
          ma            = ix(nen+1,nel)
          type = ie(ma)
          do i = 1, 10
            le(i,idCell) = e(i,ma)
          enddo
c ...
          call celllib(la    ,lx     ,lun   ,lu         ,lu0 
     .                ,lu1   ,lRo    ,lgrad ,lgrad1     ,lgrad2
     .                ,ldiv  ,liM    ,lfluxl,lw         ,ldd  
     .                ,le    ,sp     ,p     ,lsedge     ,lpedge
     .                ,lviz  ,dum    ,nen   ,nshared    ,ndm  
     .                ,type  ,iws    ,iws1  ,lib        ,nel 
     .                ,dt    ,alpha  ,dum   ,lEddyVisc  ,bs)
c .....................................................................
c
c ... calculo do residuo
          if(resid) resCell(nel)   = sp
c ... interpolacao pelo coeficiente dos momentos          
          if(momentum) then
            kk             = 2*(nst -1)
            ii             = nshared*kk  
            iM(iws1,nel)   = liM(ii+iws1)
            iM(iws1+2,nel) = liM(ii+iws1+2)
          endif  
c ... campo D            
          if(calD) then  
            kk             = nst -1
            ii             = nshared*kk + iws1  
            ddU(iws1,nel)  = ldd(ii)
          endif          
c .....................................................................
c
c ...
          if(matriz .or. forces) then
            call assbly(ad,au,al,f,ia,ja,la,p,ld,nshared,struc,unsym
     .                 ,forces,matriz)
          endif
c .....................................................................
      enddo
c .....................................................................
      endif
      return
      end
c *********************************************************************
c
c *********************************************************************
c * GFORM: Reconstruacao de gradient                                  *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c *  u     - valor da solucao a ser reconstruida por celula           *
c *  grad  - nao definido                                             *
c *  fluxl - limitador de fluxo da celula                             *
c *  x     - coordenadas nodais                                       *
c * sedge  - condicoes de contorno nas arestas                        *
c *  e     - propriedade                                              *
c *  ls    - least square                                             *
c *  ie    - tipo do elemento do material                             *
c * nelcon - vizinho da celula                                        *
c * pedge  - tipo de condicao de contorno                             *
c * ix     - conectividade dos elementos                              *
c * numel  - numero de celula                                         *
c * ndm    - numero de dimensoes                                      *
c * nen    - numero de nos por celula                                 *
c * nshared- numero de celula compartilhados pela celula central      *
c * ndf    - graus de liberdade                                       *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c *  grad  - gradiente nas celulas                                    *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine gform(u,grad,fluxl,x,sedge,e,ls,ie,nelcon,pedge,ix
     .                ,numel,ndm,nen,nshared,ndf,nst,iws,iws1,lib)  
      implicit none
      include 'openmp.fi'
      integer numel,nen,ndm,ndf,nst,nshared
      real*8 u(ndf,*),grad(ndm,*),dum,fluxl(numel)
      real*8 sedge(nst,nshared+1,*),e(10,*)
      real*8 x(ndm,*),ls(ndm,nshared,*)
      integer nelcon(nshared,*),pedge(nshared+1,*),ix(nshared+1,*),ie(*)
      integer type,iws,iws1,lib
c ... variaveis locais por elemento      
      integer lviz(6),lpedge(7)
      real*8  lgrad(3),lfluxl(2),le(10,7)
      real*8  lsedge(28) !(max 4 * 7)
      real*8  lx(168) !(max 3 * 8 * 7)
      real*8  lu(21),lls(18)
c ... variaveis locais      
      integer nel,i,ii,j,jj,k
      integer viznel,no,ma
c ........................................................................
c
c ... openmp
      if(openmpCell) then
c ... loop nas celulas
c$omp parallel private(nel,i,ii,j,jj,k,vizNel,no,ma,type)
c$omp.private(lviz,lpedge,lgrad,lfluxl,le,lsedge,lx,lu,lls)
c$omp.shared(numel,nen,ndm,ndf,nst,nshared,lib,iws,iws1)
c$omp.shared(u,grad,dum,fluxl,sedge,e,x,ls,nelcon,pedge,ix,ie) 
c$omp.num_threads(nThreadsCell)
c$omp do
        do nel = 1, numel
c ... loop da celula central
          do k = 1, ndf
            ii = nshared*ndf + k 
            lu(ii) = u(k,nel)
          enddo
c ...
          do k = 1, ndm
            lgrad(k) = grad(k,nel)
          enddo   
c ... 
          do i = 1, nen
            do k = 1, nst
              ii = (i-1)*nst + k 
              lsedge(ii) = sedge(k,i,nel)
            enddo
            lpedge(i) = pedge(i,nel)
            no        = ix(i,nel)
            ii = nshared*nen*ndm 
            do k = 1 , ndm
              jj = ii + (i-1)*ndm + k
              lx(jj) = x(k,no)
            enddo
          enddo
c ... least square          
c          do j = 1, nshared
c            do k = 1, ndm
c              ii = (j-1)*ndm + k 
c              lls(ii)=ls(k,j,nel) 
c            enddo  
c          enddo  
c ... loop na aresta        
          do j = 1 , nshared
            viznel = nelcon(j,nel)
            lviz(j)= viznel
            if( viznel .gt. 0) then
              do k = 1, ndf
                ii     = (j-1)*ndf + k 
                lu(ii) = u(k,viznel)
              enddo
c ... loop nos vertices do elemento vizinho          
              do i = 1, nen
                no = ix(i,viznel)
                ii = (j-1)*nen*ndm
                do k = 1, ndm
                  jj = ii + (i-1)*ndm + k
                  lx(jj)=x(k,no)
                enddo
              enddo  
            endif  
          enddo
          ma            = ix(nen+1,nel)
          do i = 1, 10
            le(i,nshared+1) = e(i,ma)
          enddo  
          type = ie(ma)
c ...
          call celllib(dum    ,lx     ,dum    ,lu     ,dum
     .                ,dum    ,dum    ,lgrad  ,dum    ,dum  
     .                ,dum    ,dum    ,lfluxl ,dum    ,dum  
     .                ,le     ,dum    ,dum    ,lsedge ,lpedge
     .                ,lviz   ,lls    ,nen    ,nshared,ndm 
     .                ,type   ,iws    ,iws1   ,lib    ,nel
     .                ,dum    ,dum    ,dum    ,dum    ,.false.)
c .....................................................................
          do k =1 , ndm
            grad(k,nel) = lgrad(k)
          enddo
          fluxl(nel) = lfluxl(1)
        enddo
c$omp enddo
c$omp end parallel
c .....................................................................
c
c ... sequencial
      else
c ... loop nas celulas
        do nel = 1, numel
c ... loop da celula central
          do k = 1, ndf
            ii = nshared*ndf + k 
            lu(ii) = u(k,nel)
          enddo
c ...
          do k = 1, ndm
            lgrad(k) = grad(k,nel)
          enddo   
c ... 
          do i = 1, nen
            do k = 1, nst
              ii = (i-1)*nst + k 
              lsedge(ii) = sedge(k,i,nel)
            enddo
            lpedge(i) = pedge(i,nel)
            no        = ix(i,nel)
            ii = nshared*nen*ndm 
            do k = 1 , ndm
              jj = ii + (i-1)*ndm + k
              lx(jj) = x(k,no)
            enddo
          enddo
c ... least square      
c          do j = 1, nshared
c            do k = 1, ndm
c              ii = (j-1)*ndm + k 
c              lls(ii)=ls(k,j,nel) 
c            enddo  
c          enddo  
c ... loop na aresta        
          do j = 1 , nshared
            viznel = nelcon(j,nel)
            lviz(j)= viznel
            if( viznel .gt. 0) then
              do k = 1, ndf
                ii     = (j-1)*ndf + k 
                lu(ii) = u(k,viznel)
              enddo
c ... loop nos vertices do elemento vizinho          
              do i = 1, nen
                no = ix(i,viznel)
                ii = (j-1)*nen*ndm
                do k = 1, ndm
                  jj = ii + (i-1)*ndm + k
                  lx(jj)=x(k,no)
                enddo
              enddo  
            endif  
          enddo
          ma            = ix(nen+1,nel)
          do i = 1, 10
            le(i,nshared+1) = e(i,ma)
          enddo  
          type = ie(ma)
c ...
         call celllib(dum    ,lx     ,dum    ,lu     ,dum 
     .               ,dum    ,dum    ,lgrad  ,dum    ,dum  
     .               ,dum    ,dum    ,lfluxl ,dum    ,dum
     .               ,le     ,dum    ,dum    ,lsedge ,lpedge
     .               ,lviz   ,lls    ,nen    ,nshared,ndm   
     .               ,type   ,iws    ,iws1   ,lib    ,nel 
     .               ,dum    ,dum    ,dum    ,dum    ,.false.)
c .....................................................................
          do k =1 , ndm
            grad(k,nel) = lgrad(k)
          enddo
          fluxl(nel) = lfluxl(1)
        enddo
c .....................................................................
      endif
      return
      end
c *********************************************************************
c
c *********************************************************************
c * CELLPARAMETER: calulo de parametros por celula                    *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c *  mP    - paramentros por celula (cfl, Reynalds, peclet)           *
c *  rO    - massa especifica da celula(pi+1,pt)                      *   
c *  x     - coordenadas nodais                                       *
c * eddyVis- viscosidade trubulenta                                   *
c * gradU1 - gradiente de u1                                          *
c * gradU2 - gradiente de u2                                          *
c *  e     - propriedade                                              *
c * pedge  - arestas com condicao de contorno                         *
c * sedge -  condicao de contorno nas arestas                         *
c *  w     - campo de velocidade por celula                           *
c *  ie    - tipo do elemento do material                             *
c * nelcon - vizinho da celula                                        *
c * ix     - conectividade dos elementos                              *
c * ndm    - numero de dimensoes                                      *
c * nen    - numero de nos por celula                                 *
c * nshared- numero de celula compartilhados pela celula central      *
c * nst    -                                                          *
c * dt     - passo de tempo                                           *
c * alpha  -                                                          *
c * iws    - instrucao de montagem                                    *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * mP     - paramentros por celula (cfl, Reynalds, peclet)           *
c * ----------------------------------------------------------------- *
c * OBS: (cfl,Reynalds,volume,vm) lib 1                               *
c * OBS:                 para lib 2                                   *
c *********************************************************************
      subroutine cellParameter(mP      ,rO     ,w      ,x  ,eddyVisc
     .                        ,gradU1  ,gradU2
     .                        ,e       ,sedge  ,pedge  ,ie ,nelcon
     .                        ,ix      ,numel  ,ndm    ,nen,nshared 
     .                        ,nst     ,dt     ,iws    ,lib)
      implicit none
      include 'openmp.fi'
      integer numel,nen,ndm,nshared
      integer nel,i,ii,j,jj,k,kk
      real*8 w(ndm,*),mP(20,*),sedge(nst,nshared+1,numel)
      real*8 x(ndm,*),e(10,*),dt,rO(3,*),eddyVisc(*)
      real*8  gradU1(ndm,*),gradU2(ndm,*)
      integer nelcon(nshared,*),pedge(nshared+1,numel)
      integer ix(nshared+1,numel),ie(*)
c ... variaveis locais por elemento      
      real*8 lx(168),le(10,9)
      real*8 lw(21),lsedge(28)
c ... variaveis locais
      real*8 lRo(3,7),lmP(20),lEddyVisc(7)
      real*8 lgrad1(21),lgrad2(21)
      integer lviz(7),lpedge(7)
c ...
      integer idCell
c ...
      real*8  ddum
      integer idum
      logical ldum      
c .....................................................................
      integer viznel,no,ma,nst
      integer type,iws,lib
c ... openmp
      idCell = nshared + 1
      if(openmpCell)Then
c$omp parallel private(nel,i,ii,j,jj,k,kk,vizNel,no,ma,type)
c$omp.private(lx,le,lw,lsedge,lRo,lmP,lEddyVisc,lgrad1,lgrad2)
c$omp.private(lviz,lpedge)
c$omp.private(ddum,idum,ldum)
c$omp.shared(numel,nen,ndm,nst,nshared,idCell)
c$omp.shared(w,mP,x,eddyVisc,e,dt,Ro,sedge,nelcon,pedge)
c$omp.shared(gradU1,gradU2,ix,ie) 
c$omp.num_threads(nThreadsCell)
c$omp do
c ... zerando os arranjos locais
        do nel = 1, numel
          lmP(1:20)= 0.0d0
          do i = 1, nshared+1
            lrO(1,i)  = 0.0d0
            lrO(2,i)  = 0.0d0
            lrO(3,i)  = 0.0d0
            do k = 1,ndm
              ii         = (i-1)*ndm + k
              lw(ii)     = 0.d0 
            enddo
            do j = 1, nen
              ii = (i-1)*nen*ndm
              do k = 1, ndm
                jj = ii + (j-1)*ndm + k
                lx(jj) = 0.d0
              enddo
            enddo 
          enddo
c .....................................................................        
c
c ... loop da celula central
          lEddyVisc(idCell) = eddyVisc(nel) 
          do i = 1, nshared + 1 
            do k = 1, nst
              ii      = (i-1)*nst + k  
              lsedge(ii) = sedge(k,i,nel)
            enddo
            lpedge(i)    = pedge(i,nel)
          enddo
c ......................................................................
c
c ...
          lrO(1,nshared+1)  = rO(1,nel)
          lrO(2,nshared+1)  = rO(2,nel)
          lrO(3,nshared+1)  = rO(3,nel)
          do k = 1, ndm
            ii         = nshared*ndm + k
            lw(ii)     =     w(k,nel)
            lgrad1(ii) = gradU1(k,nel)
            lgrad2(ii) = gradU2(k,nel)
          enddo       
c .....................................................................
c
c ...
          do i = 1, nen
            no        = ix(i,nel)
            ii = nshared*nen*ndm 
            do k = 1, ndm
              jj     = ii + (i-1)*ndm + k
              lx(jj) = x(k,no)
            enddo
          enddo
c .....................................................................
c
c ... loop na aresta        
          do j = 1 , nshared
            viznel    = nelcon(j,nel)
            lviz(j)   = viznel 
            if( viznel .gt. 0) then
              lEddyVisc(j)  = eddyVisc(vizNel) 
              lrO(1,j)     = rO(1,viznel)
              lrO(2,j)     = rO(2,viznel)
                     ma    = ix(nen+1,viznel)
              do kk = 1, 10
                le(kk,j) = e(kk,ma)
              enddo  
              do k = 1, ndm
                ii          = (j-1)*ndm + k
                lw(ii)      =    w(k,vizNel)
              enddo
c ... loop nos vertices do elemento vizinho          
              do i = 1, nen
                no = ix(i,vizNel)
                ii = (j-1)*nen*ndm
                do k = 1, ndm
                  jj     = ii + (i-1)*ndm + k
                  lx(jj) = x(k,no)
                enddo
              enddo  
            endif  
          enddo
          ma            = ix(nen+1,nel)
          type = ie(ma)
          do i = 1, 10
            le(i,nshared+1) = e(i,ma)
          enddo
c ...
          call celllib(ddum   ,lx     ,ddum   ,ddum     ,ddum
     .                ,ddum   ,lRo    ,lgrad1 ,lgrad2   ,ddum  
     .                ,ddum   ,ddum   ,ddum   ,lw       ,ddum  
     .                ,le     ,ddum   ,ddum   ,lsedge   ,lpedge
     .                ,lviz   ,ddum   ,nen    ,nshared  ,ndm   
     .                ,type   ,iws    ,idum   ,lib      ,nel   
     .                ,dt     ,ddum   ,lmP    ,lEddyVisc,ldum)
c .....................................................................
c
c ...
          do i = 1, 20
            Mp(i,nel) = lmP(i)
          enddo
c .................................................................
        enddo
c$omp end parallel
c .....................................................................
c
c ... sequencial
      else
c ... zerando os arranjos locais
        do nel = 1, numel
          lmP(1:20)= 0.0d0
          do i = 1, nshared+1
            lrO(1,i)  = 0.0d0
            lrO(2,i)  = 0.0d0
            lrO(3,i)  = 0.0d0
            do k = 1,ndm
              ii         = (i-1)*ndm + k
              lw(ii)     = 0.d0 
            enddo
            do j = 1, nen
              ii = (i-1)*nen*ndm
              do k = 1, ndm
                jj = ii + (j-1)*ndm + k
                lx(jj) = 0.d0
              enddo
            enddo 
          enddo
c .....................................................................        
c
c ... loop da celula central
          lEddyVisc(idCell) = eddyVisc(nel) 
          do i = 1, nshared + 1 
            do k = 1, nst
              ii      = (i-1)*nst + k  
              lsedge(ii) = sedge(k,i,nel)
            enddo
            lpedge(i)    = pedge(i,nel)
          enddo
c ......................................................................
c
c ...
          lrO(1,nshared+1)  = rO(1,nel)
          lrO(2,nshared+1)  = rO(2,nel)
          lrO(3,nshared+1)  = rO(3,nel)
          do k = 1, ndm
            ii         = nshared*ndm + k
            lw(ii)     =     w(k,nel)
            lgrad1(ii) = gradU1(k,nel)
            lgrad2(ii) = gradU2(k,nel)
          enddo       
c .....................................................................
c
c ...
          do i = 1, nen
            no        = ix(i,nel)
            ii = nshared*nen*ndm 
            do k = 1, ndm
              jj     = ii + (i-1)*ndm + k
              lx(jj) = x(k,no)
            enddo
          enddo
c .....................................................................
c
c ... loop na aresta        
          do j = 1 , nshared
            viznel    = nelcon(j,nel)
            lviz(j)   = viznel 
            if( viznel .gt. 0) then
              lEddyVisc(j)  = eddyVisc(vizNel) 
              lrO(1,j)     = rO(1,viznel)
              lrO(2,j)     = rO(2,viznel)
                     ma    = ix(nen+1,viznel)
              do kk = 1, 10
                le(kk,j) = e(kk,ma)
              enddo  
              do k = 1, ndm
                ii          = (j-1)*ndm + k
                lw(ii)      =    w(k,vizNel)
              enddo
c ... loop nos vertices do elemento vizinho          
              do i = 1, nen
                no = ix(i,vizNel)
                ii = (j-1)*nen*ndm
                do k = 1, ndm
                  jj     = ii + (i-1)*ndm + k
                  lx(jj) = x(k,no)
                enddo
              enddo  
            endif  
          enddo
          ma            = ix(nen+1,nel)
          type = ie(ma)
          do i = 1, 10
            le(i,nshared+1) = e(i,ma)
          enddo
c ...
          call celllib(ddum   ,lx     ,ddum   ,ddum     ,ddum
     .                ,ddum   ,lRo    ,lgrad1 ,lgrad2   ,ddum  
     .                ,ddum   ,ddum   ,ddum   ,lw       ,ddum  
     .                ,le     ,ddum   ,ddum   ,lsedge   ,lpedge
     .                ,lviz   ,ddum   ,nen    ,nshared  ,ndm   
     .                ,type   ,iws    ,idum   ,lib      ,nel   
     .                ,dt     ,ddum   ,lmP    ,lEddyVisc,ldum)
c .....................................................................
c
c ...
          do i = 1, 20
            Mp(i,nel) = lmP(i)
          enddo
c .................................................................
        enddo
      endif
c .....................................................................   
      return
      end
c *********************************************************************
c
c *********************************************************************
c * POSPROVELOCITYFIELD: pos processador do campo de veocidades       *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c *  x     - coordenadas nodais                                       *
c *  e     - propriedade                                              *
c *  w     - campo de velocidade por celula                           *
c *  ie    - tipo do elemento do material                             *
c * nelcon - vizinho da celula                                        *
c * ix     - conectividade dos elementos                              *
c * wP     - nao definido                                             *
c * numel  - numero de celula                                         *
c * ndm    - numero de dimensoes                                      *
c * nen    - numero de nos por celula                                 *
c * nshared- numero de celula compartilhados pela celula central      *
c * iws    - instrucao de montagem                                    *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * wP     - campo de velocidade de processado                        *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine posProVelocityField(w,x,e,ie,nelcon,ix,wP,pedge,numel
     .                              ,ndm,nen,nshared,iws,lib)
      implicit none
      include 'openmp.fi'
      integer numel,nen,ndm,nshared
      real*8 w(ndm,*),wP(ndm,*)
      real*8 x(ndm,*),e(10,*)
      integer nelcon(nshared,*),pedge(nshared+1,numel)
      integer ix(nshared+1,numel),ie(*)
c ... variaveis locais por elemento      
      real*8 lx(168),le(10,7)
      real*8 lw(21)
c ... variaveis locai
      integer lviz(7),lpedge(7)
c ...
      integer nel,i,j,k,ii,jj
c .....................................................................
      real*8  ddum
      integer idum
      logical ldum      
c .....................................................................
      integer viznel,no,ma
      integer type,iws,lib
c ... openmp
      if(openmpCell) then
c$omp parallel private(nel,i,ii,j,jj,k,vizNel,no,ma,type)
c$omp.private(lx,le,lw,lviz,lpedge,ddum,idum,ldum)
c$omp.shared(w,wp,x,e,nelcon,pedge,ix,ie)
c$omp.shared(numel,nen,ndm,nshared) num_threads(nThreadsCell)
c ... zerando os arranjos locais
c$omp do
        do nel = 1, numel
c          do i = 1, nshared+1
c            do k = 1,ndm
c              ii         = (i-1)*ndm + k
c              lw(ii)     = 0.d0 
c            enddo
c            do j = 1, nen
c              ii = (i-1)*nen*ndm
c              do k = 1, ndm
c                jj = ii + (j-1)*ndm + k
c                lx(jj) = 0.d0
c              enddo
c            enddo  
c          enddo
c .....................................................................        
c
c ... loop da celula central
          do k = 1, ndm
            ii         = nshared*ndm + k
            lw(ii)     =     w(k,nel)  
          enddo
c .....................................................................
c
c ...
          do i = 1, nshared + 1 
            lpedge(i)     = pedge(i,nel)
          enddo        
          do i = 1, nen
            no        = ix(i,nel)
            ii = nshared*nen*ndm 
            do k = 1, ndm
              jj     = ii + (i-1)*ndm + k
              lx(jj) = x(k,no)
            enddo
          enddo
c .....................................................................
c
c ... loop na aresta        
          do j = 1 , nshared
            viznel    = nelcon(j,nel)
            lviz(j)   = viznel 
            if( viznel .gt. 0) then
              ma = ix(nen+1,viznel)
              do k = 1, 10
                le(k,j) = e(k,ma)
              enddo 
              do k = 1, ndm
                ii          = (j-1)*ndm + k
                lw(ii)      =    w(k,vizNel)
              enddo
c ... loop nos vertices do elemento vizinho          
              do i = 1, nen
                no = ix(i,vizNel)
                ii = (j-1)*nen*ndm
                do k = 1, ndm
                  jj     = ii + (i-1)*ndm + k
                  lx(jj) = x(k,no)
                enddo
              enddo  
            endif  
          enddo
          ma            = ix(nen+1,nel)
          type = ie(ma)
          do i = 1, 10
            le(i,nshared+1) = e(i,ma)
          enddo
c ...
          call celllib(ddum   ,lx     ,ddum   ,ddum   ,ddum 
     .                ,ddum   ,ddum   ,ddum   ,ddum   ,ddum
     .                ,ddum   ,ddum   ,ddum   ,lw     ,ddum 
     .                ,le     ,ddum   ,ddum   ,ddum   ,lpedge 
     .                ,lviz   ,ddum   ,nen    ,nshared,ndm  
     .                ,type   ,iws    ,idum   ,lib    ,nel   
     .                ,ddum   ,ddum   ,ddum   ,ddum   ,ldum)
c .....................................................................
c
c ...
         do k = 1, ndm
           ii        = nshared*ndm + k
           wP(k,nel) = lw(ii)
         enddo   
c .....................................................................
        enddo
c$omp end parallel
c .....................................................................
c
c ... sequencial      
      else
        do nel = 1, numel
c         do i = 1, nshared+1
c            do k = 1,ndm
c              ii         = (i-1)*ndm + k
c              lw(ii)     = 0.d0 
c            enddo
c            do j = 1, nen
c              ii = (i-1)*nen*ndm
c              do k = 1, ndm
c                jj = ii + (j-1)*ndm + k
c                lx(jj) = 0.d0
c              enddo
c            enddo  
c          enddo
c .....................................................................
c
c ... loop da celula central
          do k = 1, ndm
            ii         = nshared*ndm + k
            lw(ii)     =     w(k,nel)  
          enddo
c .....................................................................
c
c ...
          do i = 1, nshared + 1 
            lpedge(i)     = pedge(i,nel)
          enddo        
          do i = 1, nen
            no        = ix(i,nel)
            ii = nshared*nen*ndm 
            do k = 1, ndm
              jj     = ii + (i-1)*ndm + k
              lx(jj) = x(k,no)
            enddo
          enddo
c .....................................................................
c
c ... loop na aresta        
          do j = 1 , nshared
            viznel    = nelcon(j,nel)
            lviz(j)   = viznel 
            if( viznel .gt. 0) then
              ma = ix(nen+1,viznel)
              do k = 1, 10
                le(k,j) = e(k,ma)
              enddo 
              do k = 1, ndm
                ii          = (j-1)*ndm + k
                lw(ii)      =    w(k,vizNel)
              enddo
c ... loop nos vertices do elemento vizinho          
              do i = 1, nen
                no = ix(i,vizNel)
                ii = (j-1)*nen*ndm
                do k = 1, ndm
                  jj     = ii + (i-1)*ndm + k
                  lx(jj) = x(k,no)
                enddo
              enddo  
            endif  
          enddo
          ma            = ix(nen+1,nel)
          type = ie(ma)
          do i = 1, 10
            le(i,nshared+1) = e(i,ma)
          enddo
c ...
          call celllib(ddum   ,lx     ,ddum   ,ddum   ,ddum 
     .                ,ddum   ,ddum   ,ddum   ,ddum   ,ddum
     .                ,ddum   ,ddum   ,ddum   ,lw     ,ddum 
     .                ,le     ,ddum   ,ddum   ,ddum   ,lpedge 
     .                ,lviz   ,ddum   ,nen    ,nshared,ndm  
     .                ,type   ,iws    ,idum   ,lib    ,nel   
     .                ,ddum   ,ddum   ,ddum   ,ddum   ,ldum)
c .....................................................................
c
c ...
         do k = 1, ndm
           ii        = nshared*ndm + k
           wP(k,nel) = lw(ii)
         enddo   
c .................................................................
        enddo
      endif
c .....................................................................
      return
      end
c *********************************************************************
c
c *********************************************************************
c * CLOSEDPRESURE: calculo o fator para pressap                       *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * temp   - temperatura                                              *
c *  e     - propriedade                                              *
c *  w     - campo de velocidade por celula                           *
c *  ie    - tipo do elemento do material                             *
c * nelcon - vizinho da celula                                        *
c * ix     - conectividade dos elementos                              *
c * cPc    - nao definido                                             *
c * numel  - numero de celula                                         *
c * ndm    - numero de dimensoes                                      *
c * nen    - numero de nos por celula                                 *
c * nshared- numero de celula compartilhados pela celula central      *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * cPc    - somatorio de   V/temp                                    *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine closedPresure(temp,x,e,ie,nelcon,ix,cPc,numel
     .                        ,ndm,nen,nshared,iws,lib)
      implicit none
      include 'openmp.fi'
      integer numel,nen,ndm,nshared
      real*8 temp(*)
      real*8 x(ndm,*),e(10,*),cPc
      integer nelcon(nshared,*)
      integer ix(nshared+1,numel),ie(*)
c ... variaveis locais por elemento      
      real*8 lx(168),le(10,7)
      real*8 lTemp
c ... variaveis locai
      integer lviz(7)
c ...
      integer nel,i,j,k,ii,jj
c .....................................................................
      real*8  ddum
      integer idum
      logical ldum      
c .....................................................................
      integer viznel,no,ma
      integer type,iws,lib
c ... openmp
      cPc = 0.d0
      if(openmpCell) then
c$omp parallel private(nel,i,ii,j,jj,k,vizNel,no,ma,type)
c$omp.private(lx,le,lviz,ddum,idum,ldum,ltemp)
c$omp.shared(x,e,nelcon,ix,ie,cPc)
c$omp.shared(numel,nen,ndm,nshared,iws,lib) num_threads(nThreadsCell)

c ... zerando os arranjos locais
c$omp do
        do nel = 1, numel
c          do i = 1, nshared+1
c            do j = 1, nen
c              ii = (i-1)*nen*ndm
c              do k = 1, ndm
c                jj = ii + (j-1)*ndm + k
c                lx(jj) = 0.d0
c              enddo
c            enddo  
c          enddo
c .....................................................................        
c
c ... loop da celula central
          do k = 1, ndm
            lTemp  = temp(nel) 
          enddo
c .....................................................................
c
c ...
          do i = 1, nen
            no        = ix(i,nel)
            ii = nshared*nen*ndm 
            do k = 1, ndm
              jj     = ii + (i-1)*ndm + k
              lx(jj) = x(k,no)
            enddo
          enddo
c .....................................................................
c
c ... loop na aresta        
          do j = 1 , nshared
            viznel    = nelcon(j,nel)
            lviz(j)   = viznel 
            if( viznel .gt. 0) then
              ma = ix(nen+1,viznel)
              do k = 1, 10
                le(k,j) = e(k,ma)
              enddo 
c ... loop nos vertices do elemento vizinho          
              do i = 1, nen
                no = ix(i,vizNel)
                ii = (j-1)*nen*ndm
                do k = 1, ndm
                  jj     = ii + (i-1)*ndm + k
                  lx(jj) = x(k,no)
                enddo
              enddo  
            endif  
          enddo
          ma            = ix(nen+1,nel)
          type = ie(ma)
          do i = 1, 10
            le(i,nshared+1) = e(i,ma)
          enddo
c ...
          call celllib(ddum   ,lx     ,ddum   ,ddum   ,ddum 
     .                ,ddum   ,ddum   ,ddum   ,ddum   ,ddum  
     .                ,ddum   ,ddum   ,ddum   ,ddum   ,ddum
     .                ,le     ,ltemp  ,ddum   ,ddum   ,idum   
     .                ,lviz   ,ddum   ,nen    ,nshared,ndm
     .                ,type   ,iws    ,idum   ,lib    ,nel  
     .                ,ddum   ,ddum   ,ddum   ,ddum   ,ldum)
c .....................................................................
c
c ...
c$omp critical
           cPc      =  cPc + ltemp
c$omp end critical
c 
c .....................................................................
        enddo
c$omp end parallel
c .....................................................................
c
c ... sequencial      
      else
        do nel = 1, numel
c          do i = 1, nshared+1
c            do j = 1, nen
c              ii = (i-1)*nen*ndm
c              do k = 1, ndm
c                jj = ii + (j-1)*ndm + k
c                lx(jj) = 0.d0
c              enddo
c            enddo  
c          enddo
c .....................................................................        
c
c ... loop da celula central
          do k = 1, ndm
            lTemp  = temp(nel) 
          enddo
c .....................................................................
c
c ...
          do i = 1, nen
            no        = ix(i,nel)
            ii = nshared*nen*ndm 
            do k = 1, ndm
              jj     = ii + (i-1)*ndm + k
              lx(jj) = x(k,no)
            enddo
          enddo
c .....................................................................
c
c ... loop na aresta        
          do j = 1 , nshared
            viznel    = nelcon(j,nel)
            lviz(j)   = viznel 
            if( viznel .gt. 0) then
              ma = ix(nen+1,viznel)
              do k = 1, 10
                le(k,j) = e(k,ma)
              enddo 
c ... loop nos vertices do elemento vizinho          
              do i = 1, nen
                no = ix(i,vizNel)
                ii = (j-1)*nen*ndm
                do k = 1, ndm
                  jj     = ii + (i-1)*ndm + k
                  lx(jj) = x(k,no)
                enddo
              enddo  
            endif  
          enddo
          ma            = ix(nen+1,nel)
          type = ie(ma)
          do i = 1, 10
            le(i,nshared+1) = e(i,ma)
          enddo
c ...
          call celllib(ddum   ,lx     ,ddum   ,ddum   ,ddum 
     .                ,ddum   ,ddum   ,ddum   ,ddum   ,ddum  
     .                ,ddum   ,ddum   ,ddum   ,ddum   ,ddum
     .                ,le     ,ltemp  ,ddum   ,ddum   ,idum   
     .                ,lviz   ,ddum   ,nen    ,nshared,ndm
     .                ,type   ,iws    ,idum   ,lib    ,nel  
     .                ,ddum   ,ddum   ,ddum   ,ddum   ,ldum)
c .....................................................................
c
c ...
           cPc      = cPc + ltemp
c .....................................................................
        enddo
       endif
c .....................................................................
      return
      end
c *********************************************************************
c
c *********************************************************************
c * LSFORM: Calculo da matriz Least square                            *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c *  u     - temperatura na celula                                    *
c *  un    - temperatura notal                                        *
c *  grad  - nao definido                                             *
c *  x     - coordenadas nodais                                       *
c * sedge  - condicoes de contorno nas arestas                        *
c *  ie    - tipo do elemento do material                             *
c * nelcon - vizinho da celula                                        *
c * pedge  - tipo de condicao de contorno                             *
c * ix     - conectividade dos elementos                              *
c * numel  - numero de celula                                         *
c * ndm    - numero de dimensoes                                      *
c * nen    - numero de nos por celula                                 *
c * nshared- numero de celula compartilhados pela celula central      *
c * ndf    - graus de liberdade                                       *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c *  grad  - gradiente nas celulas                                    *
c * ----------------------------------------------------------------- *
c *********************************************************************
c      subroutine lsform(ls,x,sedge,ie,nelcon,pedge,ix,lsedge,lx
c     .                 ,lpedge,lls,lviz,numel,ndm,nen,nshared)
c      implicit none
c      integer numel,nen,ndm,nshared
c      integer nel,i,j,k
c      real*8 dum
c      real*8 sedge(nshared+1,*),ls(ndm,nshared,*)
c      real*8 x(ndm,*)
c      integer nelcon(nshared,*),pedge(nshared+1,*),ix(nshared+1,*)
c ... variaveis locais por elemento      
c      real*8 lsedge(nshared+1),lls(ndm,nshared)
c      real*8 lx(ndm,nen,nshared+1),le(10,nshared+1)
c      integer viznel,no,ma,lpedge(nshared+1),ie(*)
c      integer lviz(nshared),type
c      call azero(le,10*(nshared+1))
c ... loop nas celulas
c      do nel = 1, numel
c ... zerando os arranjos locais
c        do i = 1, nshared+1
c          do j = 1, nen
c            do k = 1, ndm
c              lx(k,j,i) = 0.d0
c            enddo
c          enddo  
c        enddo
c        do i = 1, nshared
c          do j = 1, ndm
c            lls(j,i) = 0.d0
c          enddo  
c        enddo
c .....................................................................
c
c ... loop da celula central
c        do i = 1, nen
c          lsedge(i) = sedge(i,nel)
c          lpedge(i) = pedge(i,nel)
c          no        = ix(i,nel)
c          do k = 1 , ndm
c            lx(k,i,nshared+1) = x(k,no)
c         enddo
c       enddo
c ... loop na aresta        
c        do j = 1 , nshared
c          viznel = nelcon(j,nel)
c          lviz(j)= viznel 
c          if( viznel .gt. 0) then
c            ma = ix(nen+1,viznel)
c ... loop nos vertices do elemento vizinho          
c            do i = 1, nen
c              no = ix(i,viznel)
c              do k = 1, ndm
c                lx(k,i,j)=x(k,no)
c              enddo
c            enddo  
c          endif  
c        enddo
c        ma            = ix(nen+1,nel)
c        type = ie(ma)
c        call celllib(dum    ,lx    ,dum ,dum,dum
c     .              ,dum    ,dum   ,dum ,dum,dum
c     .              ,lsedge ,lpedge,lviz,lls,nen
c     .              ,nshared,ndm   ,type,  3,dum
c     .              ,dum    ,nel   ,dum ,dum,.false.)
c .....................................................................
c        do j = 1 , nshared
c          do k =1 , ndm
c            ls(k,j,nel) = lls(k,j)
c          enddo 
c        enddo
c      enddo
c      return
c      end
c *********************************************************************
c
c *********************************************************************
c * UFORMNODE: calcula a media nodal das temperaturas                 *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * un     - nao definido                                             *
c *  u     - variavel na celula                                       *
c * grad   - gradiente                                                *
c * fluxl  - limitador de fluxo                                       *
c *  x     - coordenadas                                              *
c * mdf    - nao definido                                             *
c * ix     - conectividade                                            *
c * ie     - tipo do elemento do material                             *
c * pedge  - tipo de condicao de contorno                             *
c * sedge  - condicoes de contorno nas arestas                        *
c * viz    - vizinhos                                                 *
c * md     - nao definido                                             *
c * nnode  - numero de nos                                            *
c * numel  - numero de elementos                                      *
c * ndm    - numero de dimensoes                                      *
c * nst    - graus de liperdade agrupados                             *
c * nen    - numero de nos por elementos                              *
c * ndf    - graus de liberdade                                       *
c * nshared- numero de celula compartilhados pela celula central      *
c * icod   - 1 - media aritimetica simples                            *
c *        - 2 - media aritimetica ponderada pela distancia           *
c *        - 3 - media calculada pelo gradiente                       *
c * icod2  - consitencia com a condicao de controno                   * 
c *        - 1 - velocidade                                           *
c *        - 2 - pressao                                              *
c *        - 3 - temperatrua                                          *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * un     - valores medios nodais                                    *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine uformnode(un    ,u       ,grad ,fluxl
     .                    ,x     ,mdf     ,ix   ,ie
     .                    ,md    ,pedge   ,sedge,viz  
     .                    ,nnode ,numel    
     .                    ,ndm   ,nen     ,nst 
     .                    ,ndf   ,nshared ,icod,icod2)
      implicit none
      real*8 un(ndf,*),u(ndf,*),xc(3),dx,d,mdf(*),x(ndm,*)
      real*8 fluxl(*),grad(ndf,ndm,*)
      integer nnode,numel,ndf,md(*),ix(nen+1,*),nen,ndm,nshared,nst
      integer ie(*),pedge(nshared+1,*),viz(nshared,*)
      real*8  sedge(nst,nshared+1,numel)
      integer acod,sn(2,4),type,ma
      integer i,k,j,nel,icod,icod2,no1,no2
      call mzero(md,nnode)
      call azero(mdf,nnode)
      call azero(un,ndf*nnode)
      xc(1:3) = 0.0d0
c ... media simples     
      if( icod .eq. 1) then
        do nel = 1, numel
          do j = 1, nen
            no1 = ix(j,nel)
            do i = 1, ndf
              un(i,no1) = un(i,no1) + u(i,nel)
            enddo
            md(no1) = md(no1) + 1
          enddo  
        enddo
        do j = 1, nnode  
          do i = 1, ndf
            un(i,j) = un(i,j)/md(j)
          enddo  
        enddo
c .....................................................................
c
c ... media ponderada pela distancia
      elseif(icod.eq.2) then
        do nel = 1, numel
          do k = 1, ndm 
            xc(k) = 0.0d0
          enddo
c ... calculo do centroide            
          do j = 1, nen
            no1 = ix(j,nel)
            do k = 1, ndm 
              xc(k) = xc(k) + x(k,no1)
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
c .....................................................................      
          do j = 1, nen
            no1 = ix(j,nel)
            d  = 0.0d0
            do k = 1, ndm
              dx    = x(k,no1) - xc(k)
              d     = d + dx*dx
            enddo
            d = dsqrt(d)
            do i = 1, ndf
              un(i,no1) = un(i,no1) + d*u(i,nel)
            enddo  
            mdf(no1) = mdf(no1) + d
          enddo  
        enddo
        do j = 1, nnode  
          do i = 1, ndf
            un(i,j) = un(i,j)/mdf(j)
          enddo  
        enddo
c .....................................................................
c
c ... obtida por variacao linear na celula apartir do gradiente
      elseif(icod .eq. 3) then  
        do nel = 1, numel
          do k = 1, ndm 
            xc(k) = 0.0d0
          enddo
c ... calculo do centroide            
          do j = 1, nen
            no1 = ix(j,nel)
            do k = 1, ndm 
              xc(k) = xc(k) + x(k,no1)
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
c .....................................................................      
          do j = 1, nen
            no1 = ix(j,nel)
            do i = 1, ndf
              d  = 0.0d0
              do k = 1, ndm
                dx = x(k,no1) - xc(k)
c                 d  = d + fluxl(nel)*grad(i,k,nel)*dx
                 d = d + grad(i,k,nel)*dx
              enddo
              un(i,no1) = un(i,no1) + u(i,nel) + d
            enddo
            md(no1) = md(no1) + 1
          enddo  
        enddo
        do j = 1, nnode  
          do i = 1, ndf
            un(i,j) = un(i,j)/md(j)
          enddo  
        enddo
c .....................................................................
c
      endif
c .....................................................................
c
c ... campo de de velocidade
      if(icod2 .eq. 1) then
        do nel = 1, numel
          ma   = ix(nen+1,nel)
          type = ie(ma)
          call setcell(type,sn,acod)
          do j = 1, nshared
            if( viz(j,nel) .lt. 0 ) then
              if(pedge(j,nel) .eq. 0 ) then
                no1      = ix(sn(1,j),nel)
                no2      = ix(sn(2,j),nel)
                do i = 1, ndf
                  un(i,no1) = 0.0d0
                  un(i,no2) = 0.0d0
                enddo  
              else if(pedge(j,nel) .eq. 1 ) then
                no1      = ix(sn(1,j),nel)
                no2      = ix(sn(2,j),nel)
                do i = 1, ndf
                  un(i,no1) = sedge(i,j,nel)
                  un(i,no2) = sedge(i,j,nel)
                enddo  
              else if(pedge(j,nel) .eq. 2 ) then
                no1      = ix(sn(1,j),nel)
                no2      = ix(sn(2,j),nel)
                do i = 1, ndf
                  un(i,no1) = sedge(i,j,nel)
                  un(i,no2) = sedge(i,j,nel)
                enddo  
              endif
            endif
          enddo
        enddo
c .....................................................................
c
c ... campo de pressao
      else if(icod2 .eq. 2) then
        do nel = 1, numel
          ma   = ix(nen+1,nel)
          type = ie(ma)
          call setcell(type,sn,acod)
          do j = 1, nshared
            if( viz(j,nel) .lt. 0 ) then
              if(pedge(j,nel) .eq. 1 ) then
                no1      = ix(sn(1,j),nel)
                no2      = ix(sn(2,j),nel)
                do i = 1, ndf
                  un(i,no1) = sedge(nst,j,nel)
                  un(i,no2) = sedge(nst,j,nel)
                enddo  
              endif
            endif
          enddo
        enddo
c .....................................................................
c
c ... campo de temperatura
      else if(icod2 .eq. 3) then
        do nel = 1, numel
          ma   = ix(nen+1,nel)
          type = ie(ma)
          call setcell(type,sn,acod)
          do j = 1, nshared
            if( viz(j,nel) .lt. 0 ) then
              if(pedge(j,nel) .eq. 1 ) then
                no1      = ix(sn(1,j),nel)
                no2      = ix(sn(2,j),nel)
                do i = 1, ndf
                  un(i,no1) = sedge(nst,j,nel)
                  un(i,no2) = sedge(nst,j,nel)
                enddo  
              endif
            endif
          enddo
        enddo
      endif
c .....................................................................
c
c ...
      return
      end
c **********************************************************************
c
c **********************************************************************
c * CBOUND:  : calcula a media nodal das temperaturas                 *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * snode  - condicao de contorno de temp por no                      *
c * fnode  - condicao de contorno de flux por no                      *
c * sedge  - nao definido                                             *
c * ix     - conectividade                                            *
c * pnode  - restricoes por no                                        *
c * pedge  - nao definido                                             *
c * numel  - numero de elementos                                      *
c * nen    - numero de nos por elementos                              *
c * nshared- faces(arestas) compartilhadas                            *
c * ndf    - graus de liberdade                                       *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * pedge  - restricoes por elmento                                   *
c * sedge  - valores das condicao de contorno por elemento            *
c * ----------------------------------------------------------------- *
c **********************************************************************
      subroutine cbound(snode,fnode,sedge,ix,pnode,pedge,numel,nen
     .                 ,nshared,ndf)
      implicit none
      integer nshared,nen,ndf,numel
      integer ix(nen+1,*),pnode(ndf,*),pedge(nshared+1,*)
      real*8  snode(ndf,*),fnode(ndf,*),sedge(nshared+1,*)
      integer nel,j
      integer lno,lpn(6),lpe(6)
      real*8  ltn(6),lfn(6),lse(6)
c ...
      do nel = 1, numel
        do j = 1, nen
          lno    = ix(j,nel)
          lpn(j) = pnode(1,lno)
          ltn(j) = snode(1,lno)
          lfn(j) = fnode(1,lno)
          lpe(j) = 0
          lse(j) = 0.d0
        enddo
        if( nshared .eq. 3) call boundtria(lpn,lpe,ltn,lfn,lse)  
        if( nshared .eq. 4) call boundquad(lpn,lpe,ltn,lfn,lse)  
        do j = 1, nshared
          pedge(j,nel) = pedge(j,nel) + lpe(j)
          sedge(j,nel) = sedge(j,nel) + lse(j)
        enddo  
      enddo
c .....................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
c * TFORM : calcula do flux apartir do gradiente                       *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * grad   - gradiente                                                *
c * flux   - nap definido                                             *
c * e      - propriedade                                              *
c * ie     - tipo do elemento do material                             *
c * ix     - conectividade                                            *
c * numel  - numero de elementos                                      *
c * ndm    - dimensao do problema                                     *
c * nen    - numero de nos por elementos                              *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * flux   - fluxo                                                    *
c * ----------------------------------------------------------------- *
c **********************************************************************
      subroutine tform(grad,flux,e,ix,numel,ndm,nen)
      implicit none
      real*8 grad(ndm,*),flux(ndm,*)
      real*8 e(10,*),k
      integer ix(nen+1,*),numel,ndm,nen
      integer i,j,ma
      do i = 1, numel
        ma = ix(nen+1,i)
        k  = e(1,ma)
        do j = 1, ndm
          flux(j,i) = -k*grad(j,i)
        enddo  
      enddo
      return
      end
c **********************************************************************
c
c *********************************************************************
c * gradEtoGradT : calcula do gradiente da temperatura apartir do     * 
c * gradietne da entalpia                                             * 
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * gradE  - gradiente da entalpia                                    *
c * flux   - nap definido                                             *
c * e      - propriedade                                              *
c * ie     - tipo do elemento do material                             *
c * ix     - conectividade                                            *
c * numel  - numero de elementos                                      *
c * ndm    - dimensao do problema                                     *
c * nen    - numero de nos por elementos                              *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * gradT  - gradiente da temperatura                                 *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine gradEtoGradT(gradE,gradT,e,ix,numel,ndm,nen)
      implicit none
      real*8 gradE(ndm,*),gradT(ndm,*)
      real*8 e(10,*),cp
      integer ix(nen+1,*),numel,ndm,nen
      integer i,j,ma
      do i = 1, numel
        ma = ix(nen+1,i)
        cp =  e(4,ma)
        do j = 1, ndm
          gradT(j,i) = gradE(j,i)/cp
        enddo  
      enddo
      return
      end
c *********************************************************************
c
c *********************************************************************
c * JacobionMatrixUform : matriz jacobiana da velocidade              *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * jMu    - nao definido                                             *
c * gradU1 - gradiente u1                                             *
c * gradU2 - gradiente u2                                             *
c * numel  - numero de elementos                                      *
c * ndm    - dimensao do problema                                     *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * jMu    - matrix jacobiana da velocidade                           *
c * ----------------------------------------------------------------- *
c **********************************************************************
      subroutine JacobionMatrixUform(jMu,gradU1,gradU2,numel,ndm)
      implicit none
      real*8 jMu(ndm,ndm,*),gradU1(ndm,*),gradU2(ndm,*)
      integer numel,ndm
      integer i,j
      do i = 1, numel
        do j = 1, ndm
          jMu(1,j,i) = gradU1(j,i)
          jMu(2,j,i) = gradU2(j,i)
        enddo
      enddo
      return
      end
c **********************************************************************
c
c **********************************************************************
c * tRO   : inicializa o vetor como a massa especificas               *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ro     - nao definido                                             *
c * e      - propriedade                                              *
c * ie     - tipo do elemento do material                             *
c * ix     - conectividade                                            *
c * numel  - numero de elementos                                      *
c * nen    - numero de nos por elemento                               *
c * flag   - .true. massa especifica ambienta calculada pela lei dos  *
c *          gases ideias                                             *
c *         .false. massa especifica nos tres niveis de tempo igual a *
c *         ambiental                                                 *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * ro     - massa especifica                                         *
c * ----------------------------------------------------------------- *
c **********************************************************************
      subroutine tRo(ro,e,ie,ix,numel,nen,flag)
      implicit none
      real*8 ro(3,*)
      real*8 e(10,*),k
      integer ix(nen+1,*),ie(*),numel,nen
      integer i,ma
      logical flag
      if( flag) then
        do i = 1, numel
          ma      = ix(nen+1,i)
          e(2,ma) = ro(1,i)
        enddo
      else
        do i = 1, numel
          ma      = ix(nen+1,i)
          k       =  e(2,ma)
          ro(1,i) =  k
          ro(2,i) =  k 
          ro(3,i) =  k
         enddo
      endif
      end
c **********************************************************************
c
c *********************************************************************
c * TURBULENCELES :modelo de turbulence LES                           *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c *  w      - campo de velociade                                      *
c *  p      - campo de pressao                                        *
c *  gradU1 - gradiente do campo de velocidade u1                     *
c *  gradU2 - gradiente do campo de velocidade u2                     *
c *  x      - coordenadas nodais                                      *
c * eddyVisc- viscoside turbulenta                                    *
c * cs      - coeficiente smagorinsk dinamico                         *
c * yPlus    - distancia adiminensional a parede                      *
c * sedge   - condicoes de contorno nas arestas                       *
c *  e      - propriedade                                             *
c *  ie     - tipo do elemento do material                            *
c * nelcon  - vizinho da celula                                       *
c * pedge   - tipo de condicao de contorno                            *
c * ix      - conectividade dos elementos                             *
c * numel   - numero de celula                                        *
c * ndm     - numero de dimensoes                                     *
c * nen     - numero de nos por celula                                *
c * nshared - numero de celula compartilhados pela celula central     *
c * ndf     - graus de liberdade                                      *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * eddyVisc- viscoside turbulenta                                    *
c * cs      - coeficiente smagorinsk dinamico                         *
c * p       - campo de pressao modificado                             *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine turbulenceLes(w,gradU1,gradU2,x,eddyVisc
     .                        ,cs,yPlus                          
     .                        ,sedge,e,ie,nelcon,pedge,ix
     .                        ,numel,ndm,nen,nshared,ndf
     .                        ,nst,iws,iws1,lib)  
      include 'openmp.fi'
      integer numel,nen,ndm,ndf,nst,nshared
      real*8 gradU1(ndm,*),gradU2(ndm,*)
      real*8 ddum
      real*8 sedge(nst,nshared+1,*),e(10,*)
      real*8 x(ndm,*),w(ndm,*)
      real*8 eddyVisc(*),cs(*),yPlus(*)
      integer nelcon(nshared,*),pedge(nshared+1,*),ix(nshared+1,*),ie(*)
      integer type,iws,iws1,lib,idCell
c ... variaveis locais por elemento      
      integer lviz(6),lpedge(7)
      real*8  lgradU1(21),lgradU2(21),le(10,7)
      real*8  lsedge(28) !(max 4 * 7)
      real*8  lx(168)    !(max 3 * 8 * 7)
      real*8  lw(21)     !(max 3 * 7)
      real*8  leddyVisc(7),lCs(7)
      real*8  lyPlus
c ... variaveis locais      
      integer nel,i,ii,j,jj,k,kk
      integer viznel,no,ma
c ........................................................................
c
c ... openmp
      if(openmpCell) then
        idCell = nshared + 1
c ... loop nas celulas
c$omp parallel default(none) private(nel,i,ii,j,jj,k,vizNel,no,ma,type)
c$omp.private(lviz,lpedge)
c$omp.private(lgradU1,lgradU2,le,lsedge,lx,lw)
c$omp.private(leddyVisc,lCs,lyPlus)
c$omp.shared(numel,nen,ndm,ndf,nst,idCell,nshared,lib,iws,iws1)
c$omp.shared(w,p,gradU1,gradU2,x,eddyVisc,cs,yPlus,ddum)
c$omp.shared(sedge,e,nelcon,pedge,ix,ie) 
c$omp.num_threads(nThreadsCell)
c$omp do
        do nel = 1, numel
c ... loop da celula central
          lCs(idCell)   = cs(nel)
          do k = 1, ndm
            ii         = nshared*ndm + k
            lw(ii)     =     w(k,nel)   
            lgradU1(ii) = gradU1(k,nel)
            lgradU2(ii) = gradU2(k,nel)
          enddo
c ... 
          do i = 1, nen
            do k = 1, nst
              ii = (i-1)*nst + k 
              lsedge(ii) = sedge(k,i,nel)
            enddo
            lpedge(i) = pedge(i,nel)
            no        = ix(i,nel)
            ii = nshared*nen*ndm 
            do k = 1 , ndm
              jj = ii + (i-1)*ndm + k
              lx(jj) = x(k,no)
            enddo
          enddo
c ... loop na aresta        
          do j = 1 , nshared
            viznel   = nelcon(j,nel)
            lviz(j)  = viznel
            if( viznel .gt. 0) then
              lCs(j)  = cs(viznel)
              ma = ix(nen+1,viznel)
              do kk = 1, 10
                le(kk,j) = e(kk,ma)
              enddo  
              do k = 1, ndm
                ii          = (j-1)*ndm + k
                lw(ii)      =     w(k,vizNel)
                lgradU1(ii)  = gradU1(k,vizNel)
                lgradU2(ii)  = gradU2(k,vizNel)
              enddo
c ... loop nos vertices do elemento vizinho          
              do i = 1, nen
                no = ix(i,vizNel)
                ii = (j-1)*nen*ndm
                do k = 1, ndm
                  jj     = ii + (i-1)*ndm + k
                  lx(jj) = x(k,no)
                enddo
              enddo  
            endif  
          enddo
          ma            = ix(nen+1,nel)
          do i = 1, 10
            le(i,idCell) = e(i,ma)
          enddo  
          type = ie(ma)
c ...
          call celllib(ddum   ,lx     ,ddum   ,ddum     ,lyPlus
     .                ,lCs    ,ddum   ,lgradU1,lgradU2  ,ddum   
     .                ,ddum   ,ddum   ,ddum   ,lw       ,ddum
     .                ,le     ,ddum   ,ddum   ,lsedge   ,lpedge
     .                ,lviz   ,ddum   ,nen    ,nshared  ,ndm   
     .                ,type   ,iws    ,iws1   ,lib      ,nel 
     .                ,ddum   ,ddum   ,ddum   ,leddyVisc,.false.)
c .....................................................................
           eddyVisc(nel) = leddyVisc(idCell)
                cs(nel)  = lCs(idCell)
             yPlus(nel)  = lyPlus
        enddo
c$omp enddo
c$omp end parallel
c .....................................................................
c
c ... sequencial
      else
        idCell = nshared + 1
        do nel = 1, numel
c ... loop da celula central
          lCs(idCell)   = cs(nel)
          do k = 1, ndm
            ii         = nshared*ndm + k
            lw(ii)     =     w(k,nel)   
            lgradU1(ii) = gradU1(k,nel)
            lgradU2(ii) = gradU2(k,nel)
          enddo
c ... 
          do i = 1, nen
            do k = 1, nst
              ii = (i-1)*nst + k 
              lsedge(ii) = sedge(k,i,nel)
            enddo
            lpedge(i) = pedge(i,nel)
            no        = ix(i,nel)
            ii = nshared*nen*ndm 
            do k = 1 , ndm
              jj = ii + (i-1)*ndm + k
              lx(jj) = x(k,no)
            enddo
          enddo
c ... loop na aresta        
          do j = 1 , nshared
            viznel = nelcon(j,nel)
            lviz(j)= viznel
            if( viznel .gt. 0) then
              lCs(j)   = cs(viznel)
              ma = ix(nen+1,viznel)
              do kk = 1, 10
                le(kk,j) = e(kk,ma)
              enddo  
              do k = 1, ndm
                ii          = (j-1)*ndm + k
                lw(ii)      =     w(k,vizNel)
                lgradU1(ii)  = gradU1(k,vizNel)
                lgradU2(ii)  = gradU2(k,vizNel)
              enddo
c ... loop nos vertices do elemento vizinho          
              do i = 1, nen
                no = ix(i,vizNel)
                ii = (j-1)*nen*ndm
                do k = 1, ndm
                  jj     = ii + (i-1)*ndm + k
                  lx(jj) = x(k,no)
                enddo
              enddo  
            endif  
          enddo
          ma            = ix(nen+1,nel)
          do i = 1, 10
            le(i,idCell) = e(i,ma)
          enddo  
          type = ie(ma)
c ...
          call celllib(ddum   ,lx     ,ddum   ,ddum     ,lyPlus
     .                ,lCs    ,ddum   ,lgradU1,lgradU2  ,ddum   
     .                ,ddum   ,ddum   ,ddum   ,lw       ,ddum
     .                ,le     ,ddum   ,lpres  ,lsedge   ,lpedge
     .                ,lviz   ,ddum   ,nen    ,nshared  ,ndm   
     .                ,type   ,iws    ,iws1   ,lib      ,nel 
     .                ,ddum   ,ddum   ,ddum   ,leddyVisc,.false.)
c .....................................................................
           eddyVisc(nel) = leddyVisc(idCell)
                 cs(nel) = lCs(idCell)  
             yPlus(nel)  = lyPlus        
        enddo
c .....................................................................
      endif
      return
      end
c *********************************************************************
c
c *********************************************************************
c * MEANVEL : Media temporal das velocidades                          *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * wMean - nao definido                                              *
c * w     - velocidade instantanea                                    *
c * numel - numero de elementos                                       *
c * ndm   - numero de dimensoes                                       *
c * t     - tempo total                                               *
c * dt    - delta temporal                                            *
c * fend  - true - finaliza a media                                   *
c *         false- calcula a integral temporal                        *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * wMean - media temporal das velocidade     (fend = .true.)         *
c * wMean - integral temporal das velocidades (fend = .false.)        *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine meanVel(wMean,w,numel,ndm,t,dt,fend)
      implicit none
      real(8) w(ndm,*),wMean(ndm,*),dt,t
      integer ndm,numel
      integer i,j
      logical fend
c ...
      if(fend) then
        do i = 1, numel
          do j = 1, ndm
            wMean(j,i) = wMean(j,i)/t
          enddo
        enddo
c .....................................................................
c
c ...
      else
        do i = 1, numel
          do j = 1, ndm
            wMean(j,i) = wMean(j,i) + w(j,i)*dt
          enddo
        enddo
      endif
c .....................................................................
c
c ...
      return
      end
c *********************************************************************
c
c *********************************************************************
c * VELRSM : Media temporal rsm das velocidades                       *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * wRsm  - nao definido                                              *
c * w     - velocidade instantanea                                    *
c * numel - numero de elementos                                       *
c * ndm   - numero de dimensoes                                       *
c * t     - tempo total                                               *
c * dt    - delta temporal                                            *
c * fend  - true - finaliza a media                                   *
c *         false- calcula a integral temporal                        *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * wRsm  - media rsm temporal das velocidade     (fend = .true.)     *
c * wRsm  - integral rsm temporal das velocidades (fend = .false.)    *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine velRsm(wRsm,w,numel,ndm,t,dt,fend)
      implicit none
      real(8) w(ndm,*),wRsm(ndm,*),dt,t
      integer ndm,numel
      integer i,j
      logical fend
c ...
      if(fend) then
        do i = 1, numel
          do j = 1, ndm
            wRsm(j,i) = dsqrt(wRsm(j,i)/t)
          enddo
        enddo
c .....................................................................
c
c ...
      else
        do i = 1, numel
          do j = 1, ndm
            wRsm(j,i) = wRsm(j,i) + w(j,i)*w(j,i)*dt
          enddo
        enddo
      endif
c .....................................................................
c
c ...
      return
      end
c *********************************************************************
c
c *********************************************************************
c * VARIANCIA : calculo da variancia padrão                           *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * var   - nao definico                                              *
c * wRsm  - media rsm                                                 *
c * wMeam - media simples                                             *
c * numel - numero de elementos                                       *
c * ndm   - numero de dimensoes                                       *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c * var   - variancia                                                 *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine velVariancie(var,wRsm,wMean,numel,ndm)
      implicit none
      real(8) var(ndm,*),wRsm(ndm,*),wMean(ndm,*)
      integer ndm,numel
      integer i,j
c ...
      do i = 1, numel
        do j = 1, ndm
          var(j,i) = wRsm(j,i)*wRsm(j,i) - wMean(j,i)*wMean(j,i)
        enddo
      enddo
c .....................................................................
c
c ...
      return
      end
c *********************************************************************
