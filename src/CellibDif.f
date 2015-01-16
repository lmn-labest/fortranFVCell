c **********************************************************************
c * CELL_D  : Celula 2D de com o termo de correção para o fluxo        *
c * difusivo                                                           *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * a      - nao definido                                              *
c * x      - coordenadas dos vertices da celula central e do seus      *
c * vizinhos                                                           *
c * u      - temperatura por celula                                    *
c * u0     - temperatura da celula do passo anterior                   *
c * w      - campo de velociade conhecido                              *
c * k      - propriedades da celula e do seus vizinhos                 *
c * grad   - gradiente reconstruido na celula                          *
c * fluxl  - limitador de fluxo na celulas                             *
c * sp     - nao definido                                              *
c *  p     - nao definido                                              *
c * sedge  - valores da condicoes de contorno por aresta               *
c * dt     - passo de tempo                                            *  
c * pedge  - tipo de condicao de contorno por aresta                   *
c * viz    - vizinhos da celula                                        *
c * nen    - numeros de nos por celula                                 *
c * nshared- numeros de faces por celulas                              *
c * iws    - 1 - sistema de equacoes                                   *
c *        - 2 - reconstruacao de gradiente                            *
c * nel    - numero da celula                                          *
c * sn(2,i)- nos da aresta i                                           *
c * acod   - codigo para o calculo da area                             *
c *         ( 3 - triangulo; 4 - quadrilatero)                         *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * a      - coeficientes                                              *
c * p      - residuo                                                   *
c * sp     - condicao de contorno essencial                            *
c * grad   - gradiente (GREEN GAUSS LINEAR)                            *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine cell_d_ggl(a,x,u,u0,grad,fluxl,k,sp,p,sedge,dt,pedge
     .                     ,viz,nshared,ndm,iws,nel,teta,sn,acod,bs)
      implicit none
c ... variavel externas
      real*8 a(*),sp,p,sedge(*),ap,grad(ndm,*),fluxl(nshared+1)
      real*8 teta,dt
      integer pedge(*),viz(*),iws
      integer nel,acod,nshared,ndm,sn(2,*)
c ... variavel internas
      real*8 xc(3,5),x(ndm,nshared,nshared+1)
      real*8 dfd,nk,ke,kf,cd,kcuf,gfn
      real*8 gfk,umin,umax,r,gpn,kc,m,cc,cvc,du
      real*8 um(4),kis(3,4),ksi(3,4),eta(3,4),n(3,4),meta(4),mksi(4)
      real*8 u(5),u0(5),ca(4),k(10,*),xm(3,4),gf(2),area(5),uf
      real*8 areacell,ug(4),xcg(3,4),alpha
      real*8 limitv,eps,df(3,4),par(2),aux,pl
      integer i,j,l,idcell,viznel,icod
      logical bs
      parameter (eps= 1.0d-8)
      external limitv,areacell
c ...
      idcell = nshared + 1
c ......................................................................
c
c ...
       call cellGeom2D(x,mksi,meta,ksi,eta,xc,area,n,xm,df,ca
     .                ,viz,nshared,ndm,sn,acod)
c .....................................................................
c
c ... propriedeade da celula
      kc      = k(1,nshared+1)
      cc      = k(2,nshared+1)
c .....................................................................
c
c ...
      goto(100,200,300,400) iws
  100 continue
c ...
      p  = 0.d0
      sp = 0.d0
      cd = 0.d0
      do i = 1, nshared
        viznel = viz(i)
c ... condicao de contorno      
        if( viznel .lt. 0) then
          a(i) = 0.d0
c ... fluxo prescrito 
          if(pedge(i) .eq. 0 ) then
            p  = p + meta(i)*sedge(i)
c ... temperatura prescrita            
          else 
            ap = kc*meta(i)/ca(i)
            sp = sp + ap
            p  = p  + ap*sedge(i)
          endif
        else
c ... interpolacao da propriedades
          alpha = ca(i)/mksi(i)
c           kf    = (1.0d0-alpha)*kc             + alpha*k(1,i)
c ... media harmonica  
          kf   = alpha/kc + (1.0d0-alpha)/k(1,i) 
          kf   = 1.0d0/kf 
          gf(1) = (1-alpha)*grad(1,idcell) + alpha*grad(1,i)
          gf(2) = (1-alpha)*grad(2,idcell) + alpha*grad(2,i)
          gfk   =    gf(1)*ksi(1,i)   +    gf(2)*ksi(2,i)
c ... compact stencil for face gradiente ( Darwish - 2003)
          gf(1) = gf(1) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(1,i)
          gf(2) = gf(2) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(2,i)
c ... produtos interno         
          nk    =   n(1,i)*ksi(1,i) +   n(2,i)*ksi(2,i)
          gfk   =    gf(1)*ksi(1,i) +    gf(2)*ksi(2,i)
          gfn   =      gf(1)*n(1,i) +    gf(2)*n(2,i)
c ... difusao direta
          dfd   = (kf*meta(i))/(nk*mksi(i))        
c ... fator de corecao
          cd    = kf*meta(i)*(gfn-gfk/nk)
c          write(3,'(i5,3es16.6)'),nel,dfd,cd,nk
c ...
          a(i)  = dfd
            p   = p + cd
        endif
      enddo
c ...
      aux  = teta*dt
      sp   = sp*aux
      p    =  p*aux
      a(idcell) = sp
      do i = 1, nshared
        a(i)     = a(i)      * aux
        a(idcell)= a(idcell) + a(i)
      enddo  
c ..................................................................... 
c
c
c ... matriz de massa 
      m    =  cc*area(idcell)
c ... Backward de segunda ordem 
      if(bs) m = 1.5d0*m
      a(idcell) = a(idcell) + m
c ... calculo do residuo R(i) = F - Ku(i) 
      p = p - a(idcell)*u(idcell)
      do i = 1, nshared
        p = p + a(i)*u(i)
      enddo
c .....................................................................     
      return
c ......................................................................
c
c ...
  200 continue
c ... criando celulas fantasma para o contorno  
      do i = 1, nshared
        viznel = viz(i)
        ug(i)  = u(i)
c ... contorno      
        if( viznel .lt. 0) then
c ... gerando o centroide da celula fantasmas        
          do j = 1, ndm
            xcg(j,i) = xc(j,idcell) + 2.d0*ca(i)*n(j,i)
          enddo
c ... gerando a condicao de contorno da celula fantasma
c
c ... considerando o gradiente normal para estimar o potencial no centro
c     da celula fantasmas)
          if(pedge(i) .eq. 0 ) then
            ug(i) = u(idcell) + 2.d0*sedge(i)/kc*ca(i)
c ... temperatura prescrita (extrapolacao linear para a centro do
c     potencial da celula fantasmas)
          else
            ug(i) = 2.d0*sedge(i)-u(idcell)
          endif
        else
          do j = 1, ndm
            xcg(j,i) = xc(j,i)
          enddo
        endif
      enddo  
c ... vetor que une os centroides dos vizinho 
      do i = 1, ndm
        do l = 1, nshared
c ... centroide l  
          eta(i,l) = xcg(i,sn(2,l)) - xcg(i,sn(1,l))    
        enddo
      enddo
c ... vetor normal com o modulo igual a arestas
      do l = 1, nshared
        n(1,l) = eta(2,l)
        n(2,l) =-eta(1,l)
      enddo
c .....................................................................
c ... area                 
      area(1) = areacell(eta,acod)
c .....................................................................
c
c ...
      do l = 1, nshared
        um(l) = 0.5d0*(ug(sn(2,l))+ug(sn(1,l))) 
      enddo
c .....................................................................
c
c ... Reconstrucao linear Green-Gauss 
      do i = 1, nshared
        uf = um(i)
        grad(1,1) = grad(1,1) + uf*n(1,i)
        grad(2,1) = grad(2,1) + uf*n(2,i)
      enddo
      grad(1,1) = grad(1,1) / area(1)
      grad(2,1) = grad(2,1) / area(1)
c ... funcao de limatacao do fluxo
      icod = 3
      fluxl(1) = 1.0d0
c ... obtendo o valor máximo da solucao
      umax = u(idcell)
      umin = umax
      do i = 1, nshared
        if(umax .lt. ug(i)) umax = ug(i)
        if(umin .gt. ug(i)) umin = ug(i)  
      enddo
c ... 
      par(1) = area(idcell)
      par(2) = 1.0d0
      do i = 1, nshared
        gpn = grad(1,1)*df(1,i) + grad(2,1)*df(2,i) + eps
        uf = u(idcell) + gpn 
        fluxl(2) = 1.0d0
        if( uf .gt. u(idcell) ) then
          r         = ( umax - u(idcell) ) / gpn
          fluxl(2)  = limitv(r,par,icod)
        elseif(uf.lt.u(idcell)) then
          r = ( umin - u(idcell) ) / gpn
          fluxl(2)  = limitv(r,par,icod)
        endif
c ... menor parametro
        fluxl(1) = min(fluxl(1),fluxl(2))
      enddo
c ... fluxl >= 0
      fluxl(1) = max(0.d0,fluxl(1))
      return
c .....................................................................
c
c ...
  300 continue
      return
c ... F = fontes no volume + Vpu0 - dt*(1-teta)*F 
  400 continue
c ...
      p   = 0.d0
      pl  = 0.d0
      sp  = 0.d0
      cd  = 0.d0
      cvc = 0.d0
      dfd = 0.d0
c .....................................................................
c
c ... Vp(n)u(n)
      m    = cc*area(idcell)
      p    = m*u(idcell)
c ... Backward de segunda ordem 
      if(bs) p = 2*p - 0.5d0*m*u0(idcell)
c .....................................................................
c
c ... VdtQ
      if(pedge(idcell) .eq. 1)  p = p + sedge(idcell)*area(idcell)*dt
c .....................................................................
c
c ... Euler Bakward
      if(teta .eq. 1.0d0) return
c .....................................................................
c
c ...
      do i = 1, nshared
        viznel = viz(i)
        if(viznel .lt. 0) then
c ... fluxo convectivo na face de contorno          
          if(pedge(i) .eq. 0) then
            pl = pl - meta(i)*sedge(i)
          else
            ap   = kc*meta(i)/ca(i)
            pl   = pl  - ap*(sedge(i)-u(idcell))
          endif
        else 
c ... interpolacao da propriedades
          alpha = ca(i)/mksi(i)
c           kf    = (1.0d0-alpha)*kc             + alpha*k(1,i)
c ... media harmonica  
          kf   = alpha/kc + (1.0d0-alpha)/k(1,i) 
          kf   = 1.0d0/kf 
          gf(1) = (1-alpha)*grad(1,idcell) + alpha*grad(1,i)
          gf(2) = (1-alpha)*grad(2,idcell) + alpha*grad(2,i)
c ... produtos interno         
          nk    =   n(1,i)*ksi(1,i) + n(2,i)*ksi(2,i)
          gfn   =      gf(1)*n(1,i) +  gf(2)*n(2,i)
          gfk   =    gf(1)*ksi(1,i) +  gf(2)*ksi(2,i)
c ... interpolacao de gradiente (Darwish, F. Maukalled - 2003)
          gf(1) = gf(1) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(1,i)
          gf(2) = gf(2) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(2,i)
c ... 
          gf(1) = (1-alpha)*grad(1,idcell) + alpha*grad(1,i)
          gf(2) = (1-alpha)*grad(2,idcell) + alpha*grad(2,i)
c ... difusao direta        
          dfd   = (kf*meta(i))/(nk*mksi(i))
c ... fator de corecao
          du    = u(i) - u(idcell)
          cd    = kf*meta(i)*(gfn-gfk/nk)
          pl    = pl - dfd*du + cd   
        endif
      enddo
c ......................................................................
      p = p - dt*(1.0d0-teta)*pl
c ......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
c * CELL_D_WLS   : Celula quadrilatera de com o termpo de correção     *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * a      - nao definido                                              *
c * x      - coordenadas dos vertices da celula central e do seus      *
c * vizinhos                                                           *
c * u      - temperatura por celula                                    *
c * u0     - temperatura da celula do passo anterior                   *
c * w      - campo de velociade conhecido                              *
c * k      - propriedades da celula e do seus vizinhos                 *
c * grad   - gradiente reconstruido na celula                          *
c * fluxl  - limitador de fluxo na celulas                             *
c * sp     - nao definido                                              *
c *  p     - nao definido                                              *
c * sedge  - valores da condicoes de contorno por aresta               *
c * dt     - passo de tempo                                            *  
c * pedge  - tipo de condicao de contorno por aresta                   *
c * viz    - vizinhos da celula                                        *
c * nen    - numeros de nos por celula                                 *
c * nshared- numeros de faces por celulas                              *
c * iws    - 1 - sistema de equacoes                                   *
c *        - 2 - reconstruacao de gradiente                            *
c *        - 3 - matrix least squere                                   *
c *        - 4 - matrix least square                                   *
c * nel    - numero da celula                                          *
c * sn(2,i)- nos da aresta i                                           *
c * acod   - codigo para o calculo da area                             *
c *         ( 3 - triangulo; 4 - quadrilatero)                         *                                         *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * a      - coeficientes                                              *
c * p      - residuo                                                   *
c * sp     - condicao de contorno essencial                            *
c * grad   - gradiente (WEIGHTED LEAST SQUARE)                         *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine cell_d_Wls(a,x,u,u0,grad,fluxl,k,sp,p,sedge,dt,pedge
     .                     ,viz,ls,nshared,ndm,iws,nel,teta,sn,acod,bs)
      implicit none
      real*8 areacell,gf(2),gfn,gfk,dx(2,4),ata(2,2),cc,m,dt
      real*8 a(*),sp,p(1),xc(3,5),x(ndm,nshared,nshared+1),dfd,nk,um(4)
      real*8 cd,kf
      real*8 sedge(*),ksi(3,4),eta(3,4),n(3,4),meta(4),mksi(4),area(5)
      real*8 kc,un(4),u(5),ca(4),alpha,k(10,*),xm(3,4),ap,grad(ndm,*)
      real*8 det,ls(ndm,nshared),aux,w(4),acod,fluxl(nshared+1)
      real*8 tdt,teta,u0(5),uf,df(3,4),ug(4),xcg(3,4),du,pl
      integer pedge(*),viz(*),viznel,iws,nel,sn(2,*)
      integer ndm,nen,nshared,i,j,l,idcell
      logical bs
c ...
      idcell = nshared + 1
c .....................................................................
c
c ...
       call cellGeom2D(x,mksi,meta,ksi,eta,xc,area,n,xm,df,ca
     .                ,viz,nshared,ndm,sn,acod)
c .....................................................................
c ... propriedeade da celula
      kc      = k(1,nshared+1)
      cc      = k(2,nshared+1)
c .....................................................................
c
c ... criando celulas fantasma para o contorno  
      do i = 1, nshared
        viznel = viz(i)
        ug(i) = u(i)
c ... contorno      
        if( viznel .lt. 0) then
c ... gerando o centroide da celula fantasmas        
          do j = 1, ndm
            xcg(j,i) = xc(j,idcell) + 2.d0*ca(i)*n(j,i)
          enddo
c ... gerando a condicao de contorno da celula fantasma
c ... considerando o gradiente normal para estimar o potencial no centro
c     da celula fantasmas)
          if(pedge(i) .eq. 0 ) then
            ug(i) = u(idcell) + 2.d0*sedge(i)/kc*ca(i)
c ... temperatura prescrita (extrapolacao linear para a centro o
c     potencial da celula fantasmas)
          else
            ug(i) = 2.d0*sedge(i)-u(idcell)
          endif
        else
          do j = 1, ndm
            xcg(j,i) = xc(j,i)
          enddo
        endif     
      enddo  
c .....................................................................
c
c ... 
      goto(100,200,300,400) iws
c .....................................................................
c
  100 continue
c ...
      p  = 0.d0
      sp = 0.d0
      do i = 1, nshared
        viznel = viz(i)
c ... condicao de contorno      
        if( viznel .lt. 0) then
          a(i) = 0.d0
c ... fluxo prescrito 
          if(pedge(i) .eq. 0 ) then
c        print*,'flux',nel
            p  = p + meta(i)*sedge(i)
c ... temperatura prescrita            
          else
            ap = kc*meta(i)/ca(i)
            sp = sp + ap
            p  = p  + ap*sedge(i)
          endif
        else
c ... interpolacao da propriedaes
          alpha = ca(i)/mksi(i)
          kf    = (1-alpha)*kc + alpha*k(1,i)
          gf(1) = (1-alpha)*grad(1,idcell) + alpha*grad(1,i)
          gf(2) = (1-alpha)*grad(2,idcell) + alpha*grad(2,i)
c ... produtos interno         
          nk    =   n(1,i)*ksi(1,i) +   n(2,i)*ksi(2,i)
          gfk   =    gf(1)*ksi(1,i) +    gf(2)*ksi(2,i)
          gfn   =      gf(1)*n(1,i) +    gf(2)*n(2,i)
c ... interpolacao de gradiente (Darwish, F. Maukalled - 2003)
          gf(1) = gf(1) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(1,i)
          gf(2) = gf(2) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(2,i)
c ... produtos internos
          gfn   =      gf(1)*n(1,i) +  gf(2)*n(2,i)
          gfk   =    gf(1)*ksi(1,i) +    gf(2)*ksi(2,i)
c ... fator de corecao
           cd    = kf*meta(i)*(gfn-gfk/nk)
c          write(3,'(i5,3es16.6)'),nel,dfd,cd,nk
c ... difusao direta        
          dfd   = (kf*meta(i))/(nk*mksi(i))
c ...
          a(i)  = dfd
            p   = p + cd
        endif
      enddo
c .....................................................................
c
c ...
      tdt  = teta*dt
      sp   = sp*tdt
      p    =  p*tdt
      a(idcell) = sp
      do i = 1, nshared
        a(i)     = a(i)      * tdt
        a(idcell)= a(idcell) + a(i)
      enddo  
c .....................................................................
c
c ... matriz de massa 
      m    = cc*area(idcell)
c ... Backward de segunda ordem 
      if(bs) m = 1.5d0*m
      a(idcell) = a(idcell) + m
c ... calculo do residuo R(i) = F - Ku(i) 
      p = p - a(idcell)*u(idcell)
      do i = 1, nshared
        p = p + a(i)*u(i)
      enddo
c .....................................................................     
      return
c .....................................................................
c
c ...
  200 continue
c ...
      do j = 1, nshared
        um(j) =ug(j)-u(idcell) 
      enddo
c .....................................................................
c
c ... Reconstrucao Least -Square
      grad(1,1) = 0.0d0
      grad(2,1) = 0.0d0
      do j = 1, nshared
        grad(1,1) = grad(1,1) + ls(1,j) * um(j)
        grad(2,1) = grad(2,1) + ls(2,j) * um(j)
      enddo
c ... funcao de limatacao do fluxo para o termo dvectivo 
      fluxl(1) = 1.0d0
      return
c .....................................................................
c
c ...
  300 continue
c ...
      do j = 1, nshared
         do i = 1, ndm
          dx(i,j) = xcg(i,j) - xc(i,idcell)
        enddo   
      enddo
c .....................................................................
c
c ... pesos 1.0/|r-rp|
      do j = 1, nshared
        w(j) = 1.0d0/dsqrt(dx(1,j)*dx(1,j)+dx(2,j)*dx(2,j))
      enddo
c .....................................................................
c
c ... Matriz AtWA
c     | [wj(dxj)^2  + ...]    [wj(dxjdyj) + ...] |
c     |                                          | 
c     | [wj(dxjdyj) + ...]    [wj(dyj)^2  + ...] |
      ata(1,1) = 0.0d0
      ata(2,1) = 0.0d0 
      ata(1,2) = 0.0d0 
      ata(2,2) = 0.0d0
      do j = 1, nshared
       ata(1,1) = ata(1,1) + w(j)*dx(1,j)*dx(1,j) 
       ata(1,2) = ata(1,2) + w(j)*dx(1,j)*dx(2,j)
       ata(2,2) = ata(2,2) + w(j)*dx(2,j)*dx(2,j) 
      enddo
      ata(2,1) = ata(1,2) 

c ... determinante de AtWA     
      det      = 1.0d0/(ata(1,1)*ata(2,2) - ata(1,2)*ata(2,1))
c ... Invertendo AtA
      ata(1,2) = -ata(1,2)*det
      ata(2,1) = -ata(2,1)*det
      aux      =  ata(2,2)*det
      ata(2,2) =  ata(1,1)*det
      ata(1,1) =  aux
c ...(AtWA)-1AWt
      do j = 1, nshared
        ls(1,j)  = (ata(1,1)*dx(1,j) + ata(2,1)*dx(2,j))*w(j)
        ls(2,j)  = (ata(1,2)*dx(1,j) + ata(2,2)*dx(2,j))*w(j)
      enddo
      return
c ... F = fontes no volume + Vpu0 - dt*(1-teta)*F 
  400 continue
c ...
      p   = 0.d0
      pl  = 0.d0
      sp  = 0.d0
      dfd = 0.d0
c .....................................................................
c
c ... Vp(n)u(n)
      m    = cc*area(idcell)
      p    = m*u(idcell)
c ... Backward de segunda ordem 
      if(bs) p = 2*p - 0.5d0*m*u0(idcell)
c .....................................................................
c
c ... VdtQ
      if(pedge(idcell) .eq. 1)  p = p + sedge(idcell)*area(idcell)*dt
c .....................................................................
c
c ... Euler Bakward
      if(teta .eq. 1.0d0) return
c .....................................................................
c
c ...
      do i = 1, nshared
        viznel = viz(i)
        if(viznel .lt. 0) then
c ... fluxo convectivo na face de contorno          
          if(pedge(i) .eq. 0) then
            pl = pl - meta(i)*sedge(i)
          else
            ap  = kc*meta(i)/ca(i)
            pl  = pl  - ap*(sedge(i)-u(idcell))
          endif
        else 
c ... interpolacao da propriedades
          alpha = ca(i)/mksi(i)
c           kf    = (1.0d0-alpha)*kc             + alpha*k(1,i)
c ... media harmonica  
          kf    = alpha/kc + (1.0d0-alpha)/k(1,i) 
          kf    = 1.0d0/kf 
          gf(1) = (1-alpha)*grad(1,idcell) + alpha*grad(1,i)
          gf(2) = (1-alpha)*grad(2,idcell) + alpha*grad(2,i)
c ... produtos interno         
          nk    =   n(1,i)*ksi(1,i) + n(2,i)*ksi(2,i)
          gfn   =      gf(1)*n(1,i) +  gf(2)*n(2,i)
          gfk   =    gf(1)*ksi(1,i) +  gf(2)*ksi(2,i)
c ... interpolacao de gradiente (Darwish, F. Maukalled - 2003)
          gf(1) = gf(1) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(1,i)
          gf(2) = gf(2) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(2,i)
c ...          
          gfn   =      gf(1)*n(1,i) +  gf(2)*n(2,i)
          gfk   =    gf(1)*ksi(1,i) +  gf(2)*ksi(2,i)
c ... difusao direta         
          dfd   = (kf*meta(i))/(nk*mksi(i))
c ... fator de corecao
          du    = u(i) - u(idcell)
          cd    = kf*meta(i)*(gfn-gfk/nk)
          pl    = pl - dfd*du + cd   
        endif
      enddo
c ......................................................................
      p = p - dt*(1.0d0-teta)*pl     
c ......................................................................
      return
      end
c *********************************************************************
