c **********************************************************************
c * CELL_A_EB_GGL:Celula quadrilatera com o termpo de                  *
c * correção com base na aresta (advecao pura)                         *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * a      - nao definido                                              *
c * x      - coordenadas dos vertices da celula central e do seus      *
c * vizinhos                                                           *
c * u      - temperatura por celula                                    *
c * u0     - temperatura por celula do passo anterior                  *
c * w      - campo de velociade conhecido                              *
c * grad   - gradiente reconstruido na celula                          *
c * fluxl  - limitador de fluxo na celulas                             *
c * k      - propriedades da celula e do seus vizinhos                 *
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
c * nel    - numero da celula                                          *                            *
c * sn(2,i)- nos da aresta i                                           *
c * acod   - codigo para o calculo da area                             *
c *         ( 3 - triangulo; 4 - quadrilatero)                         *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * a      - coeficientes                                              *
c * p      - residuo                                                   *
c * grad   - gradiente (GREEN GAUSS LINEAR)                            *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine cell_a_eb_ggl(a,x,u,u0,w,grad,fluxl,k
     .          ,sp,p,sedge,dt,pedge,viz,nshared,ndm,iws,nel,teta,sn
     .          ,acod,bs)
      implicit none
      real*8 areacell,uf,w(ndm,*),wfn,wf(2),fluxl(nshared+1)
      real*8 r,limit,du,rof,gpk,roc,cc,dt,df(3,4),gfk,gf(2)             
      real*8 a(*),sp,p(1),xc(3,5),x(ndm,nshared,nshared+1),um(4),cv,cvc
      real*8 sedge(*),ksi(3,4),eta(3,4),n(3,4),meta(4),mksi(4),area(5)
      real*8 u(5),ca(4),alpha,k(10,*),xm(3,4),km(3,4),grad(ndm,*)
      real*8 eps,ug(4),xcg(3,4),m,u0(5),tdt,pl,teta
      integer pedge(*),viz(*),viznel,iws,nel,l
      integer ndm,nshared,i,j,icod,acod,idcell,sn(2,*)
      logical bs
      parameter (eps= 1.0d-8)
c ... 
c ... 1 - UD - upwind de primeira ordem
c ... 2 - LUD - upwind de segunda ordem 
c ... 4 - Van Albada - TVD
c ... 5 - Mid - Mod - TVD
c ... 6 - MUSCL - TVD
c ... 7 - OSHER - TVD
      icod = 1   
      idcell = nshared + 1   
c .....................................................................
c
c ...      
      call cellGeom2D(x,mksi,meta,ksi,eta,xc,area,n,xm,df,ca
     .               ,viz,nshared,ndm,sn,acod)
c .....................................................................
c
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
c ... temperatura prescrita (extrapolacao linear para a centro o
c     potencial da celula fantasmas)
          if(pedge(i) .eq. 0 ) then
            ug(i)  = u(idcell)
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
c ... vetor que une a intersecão da reta entre os dois centroides 
c     ao ponto médio da aresta
      call vectorKm2d(xcg,xc,x,xm,km,sn,nshared,ndm)
c .....................................................................      
c
c ... propriedeade da celula
      roc     = k(1,nshared+1)
      cc      = k(2,nshared+1)
      icod    = k(3,nshared+1)
c .....................................................................
c
c ...
      goto(100,200,300,400) iws
c .....................................................................
  100 continue
c ...
      cvc = 0.d0
      do i = 1, nshared
        viznel = viz(i)
c ... condicao de contorno      
        if( viznel .lt. 0) then
          a(i) = 0.d0
c ... fluxo convectivo na face de contorno 
          wfn = w(1,idcell)*n(1,i) +  w(2,idcell)*n(2,i)         
          cv  = roc*wfn*meta(i)
c ... fluxo convectivo primeira ordem            
c ... correcao do fluxo convectivo
          if(cv.gt.0.0d0)then
            sp  = sp + cv
          else  
             p  = p  - cv*sedge(i)
          endif  
        else
c ... interpolacao da propriedaes
          alpha = area(i)/(area(i)+area(idcell)) 
          rof   = (1-alpha)*roc    + alpha*k(1,i)
          wf(1) = (1-alpha)*w(1,idcell) + alpha*w(1,i)
          wf(2) = (1-alpha)*w(2,idcell) + alpha*w(2,i)
c ... produtos interno         
          wfn   =      wf(1)*n(1,i) +  wf(2)*n(2,i)
c ... fluxo convectivo de primeira ordem
          cv = rof*wfn*meta(i)
c ... fluxo convectivo de ordem superior
          if(cv.lt.0.0d0) then
            du  = u(idcell) - u(i) + eps
            gpk = grad(1,i)*(xc(1,idcell)-xc(1,i))
     .          + grad(2,i)*(xc(2,idcell)-xc(2,i))
          else  
            du  = u(i) - u(idcell) + eps 
            gpk = (grad(1,idcell)*ksi(1,i)
     .          + grad(2,idcell)*ksi(2,i))*mksi(i)
          endif
c ... interpolacao unidirecional
          alpha = ca(i)/mksi(i)    
          r     = 2.0d0*gpk/du - 1.0d0
          cvc   = alpha*limit(r,icod)*du
c ....................................................................
c
c ... correcao nao limitada para malha na estruturadas
          gf(1) = (1-alpha)*grad(1,idcell) + alpha*grad(1,i)
          gf(2) = (1-alpha)*grad(2,idcell) + alpha*grad(2,i)
          gfk   =    gf(1)*ksi(1,i)   +    gf(2)*ksi(2,i)
c ... compact stencil for face gradiente ( Darwish - 2003)
          gf(1) = gf(1) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(1,i)
          gf(2) = gf(2) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(2,i)
c .....
          cvc   = cvc +  gf(1)*km(1,i) + gf(2)*km(2,i)
c .....................................................................
          a(i) = max(-cv,0.0d0)
          sp   = sp   + cv
          p    =  p   - cv*cvc
c          write(3,'(i6,9es16.6)')nel,cv,cvc
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
      m    = roc*cc*area(idcell)
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
c ... recontrucao de gradiente
  200 continue
c .....................................................................      
c
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
c ... funcao de limatacao do fluxo para o termo dvectivo 
      fluxl(1) = 1.0d0
      return
c .....................................................................
c
c ...
  300 continue
      return
c .....................................................................
c
c ... F = fontes no volume + Vpu0 - dt*(1-teta)*F
  400 continue
c ...
      p   = 0.d0
      pl  = 0.d0
      sp  = 0.d0
      cv  = 0.d0
      cvc = 0.d0
c .....................................................................
c
c ... Vp(n)u(n)
      m    = roc*cc*area(idcell)
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
          wfn = w(1,idcell)*n(1,i) +  w(2,idcell)*n(2,i)
          cv  = roc*wfn*meta(i)
          if(pedge(i) .eq. 0) then
            if(cv .gt. 0.0d0) pl = pl + cv*u(idcell) 
          else
c ... fluxo advectivo de primeira ordem            
            if(cv.lt.0.0d0) then
              pl = pl  + cv*sedge(i)
            endif
          endif
        else 
c ... interpolacao da propriedaes
          alpha = area(i)/(area(i)+area(idcell))
          rof   = (1-alpha)*roc       + alpha*k(1,i)
          wf(1) = (1-alpha)*w(1,idcell)    + alpha*w(1,i)
          wf(2) = (1-alpha)*w(2,idcell)    + alpha*w(2,i)
c ... produtos interno         
          wfn   =      wf(1)*n(1,i) +  wf(2)*n(2,i)
c ... fluxo convectivo de primeira ordem
          cv    = rof*wfn*meta(i)
c ... fluxo convectivo de ordem superior
          if(cv.lt.0.0d0) then
           du   = u(idcell) - u(i) + eps            
           gpk  = grad(1,i)*(xc(1,idcell)-xc(1,i))
     .          + grad(2,i)*(xc(2,idcell)-xc(2,i))
          else  
           du  = u(i) - u(idcell) + eps    
           gpk = (grad(1,idcell)*ksi(1,i)
     .         +  grad(2,idcell)*ksi(2,i))*mksi(i)
          endif
c ... interpolacao unidirecional
          alpha = ca(i)/mksi(i)              
          r   = 2.0d0*gpk/du - 1.0d0
          cvc = alpha*limit(r,icod)*du
c ... correcao nao limitada para malha na estruturadas
          gf(1) = (1-alpha)*grad(1,idcell) + alpha*grad(1,i)
          gf(2) = (1-alpha)*grad(2,idcell) + alpha*grad(2,i)
          gfk   =    gf(1)*ksi(1,i)   +    gf(2)*ksi(2,i)
c ... compact stencil for face gradiente ( Darwish - 2003)
          gf(1) = gf(1) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(1,i)
          gf(2) = gf(2) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(2,i)
c ...
          cvc   = cvc +  gf(1)*km(1,i) + gf(2)*km(2,i)
c .....................................................................
          pl = pl+max(cv,0.0d0)*u(idcell)+min(cv,0.d0)*u(i)+cv*cvc
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
c * CELL_A_VB_GGL: Celula 2D correção da baseada no volume             *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * a      - nao definido                                              *
c * x      - coordenadas dos vertices da celula central e do seus      *
c * vizinhos                                                           *
c * u      - temperatura por celula                                    *
c * u0     - temperatura por celula do passo anterior                  *
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
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * a      - coeficientes                                              *
c * p      - residuo                                                   *
c * sp     - condicao de contorno essencial                            *
c * grad   - gradiente (GREEN GAUSS LINEAR)                            *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine cell_a_vb_ggl(a,x,u,u0,w,grad,fluxl,k,sp,p,sedge
     .                        ,dt,pedge,viz,nshared,ndm,iws,nel
     .                        ,teta,sn,acod,bs)
      implicit none
      real*8 areacell,uf,w(ndm,*),wf(ndm),wfn,par(2)
      real*8 r,limitv,du,rof,gpn,roc,cc,dt,teta,df(3,4)           
      real*8 a(*),sp,p,xc(3,5),x(ndm,nshared,nshared+1),um(4),cv,cvc
      real*8 sedge(*),ksi(3,4),eta(3,4),n(3,4),meta(4),mksi(4),area(5)
      real*8 u(5),ca(4),alpha,k(10,*),xm(3,4),grad(ndm,*),pl,tdt
      real*8 eps,ug(4),xcg(3,4),m,u0(5),fluxl(nshared+1),umax,umin
      parameter ( eps= 1.0d-8)
      integer pedge(*),viz(*),viznel,iws,nel,icod
      integer ndm,nshared,i,j,l,acod,sn(2,*),idcell
      logical bs
c ...
      icod   = 1
      idcell = nshared + 1
c .....................................................................
c
c ...
      call cellGeom2D(x,mksi,meta,ksi,eta,xc,area,n,xm,df,ca
     .               ,viz,nshared,ndm,sn,acod)
c .....................................................................
c
c ... propriedeade da celula
      roc     = k(1,nshared+1)
      cc      = k(2,nshared+1)
      icod    = k(3,nshared+1)
c .....................................................................

c ... 
      goto(100,200,300,400) iws
c .....................................................................
c
  100 continue
c ...
      cvc = 0.d0
      do i = 1, nshared
        viznel = viz(i)
c ... condicao de contorno      
        if( viznel .lt. 0) then
          a(i) = 0.d0
c ... fluxo convectivo na face de contorno          
          wfn = w(1,idcell)*n(1,i) +  w(2,idcell)*n(2,i)
          cv  = roc*wfn*meta(i)
c ... fluxo convectivo primeira ordem            
          if(cv.gt.0.0d0)then
            sp  = sp + cv
          else  
            p   = p  - cv*sedge(i)
          endif  
        else
c ... interpolacao da propriedaes
          alpha = area(i)/(area(i)+area(idcell)) 
          rof   = (1-alpha)*roc    + alpha*k(1,i)
          wf(1) = (1-alpha)*w(1,idcell) + alpha*w(1,i)
          wf(2) = (1-alpha)*w(2,idcell) + alpha*w(2,i)
c ... produtos interno         
          wfn   =      wf(1)*n(1,i) +  wf(2)*n(2,i)
c ... fluxo convectivo de primeira ordem
          cv = rof*wfn*meta(i)
c ... fluxo convectivo de ordem superior
          if(cv.lt.0.0d0) then
c ... distancia normal a aresta    
            gpn  = grad(1,i)*(xm(1,i)-xc(1,i))
     .           + grad(2,i)*(xm(2,i)-xc(2,i))
            gpn  = fluxl(i) *gpn
          else  
            gpn = grad(1,idcell)*df(1,i)+grad(2,idcell)*df(2,i)
            gpn = fluxl(idcell) * gpn
          endif  
          cvc = gpn
c .....................................................................
          a(i) = max(-cv,0.0d0)
          sp   = sp   + cv
          p    =  p   - cv*cvc
c          write(3,'(i6,9es16.6)'),nel,cv,cvc,ca(i),mksi(i)
        endif
      enddo
c ....................................................................
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
      m    = roc*cc*area(idcell)
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
c ... temperatura prescrita (extrapolacao linear para a centro o
c     potencial da celula fantasmas)
          if(pedge(i) .eq. 0 ) then
            ug(i)  = u(idcell)
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
c ... funcao de limatacao do fluxo para o termo dvectivo 
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
        elseif(uf.eq.u(idcell))then
          fluxl(2)  = 1.0d0
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
c .....................................................................
c
c ...
  400 continue
c ...
      p   = 0.d0
      pl  = 0.d0
      sp  = 0.d0
      cvc = 0.d0
c .....................................................................
c
c ... Vp(n)u(n)
      m    = roc*cc*area(idcell)
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
          wfn = w(1,idcell)*n(1,i) +  w(2,idcell)*n(2,i)
          cv  = roc*wfn*meta(i)
          if(pedge(i) .eq. 0) then
            if(cv.gt.0.0d0) pl = pl + cv*u(idcell) 
          else
c ... fluxo advectivo de primeira ordem            
            if(cv.lt.0.0d0) then
              pl = pl  + cv*sedge(i)
            endif
          endif
        else 
c ... interpolacao da propriedaes
          alpha = area(i)/(area(i)+area(idcell))
          rof   = (1-alpha)*roc       + alpha*k(2,i)
          wf(1) = (1-alpha)*w(1,idcell)    + alpha*w(1,i)
          wf(2) = (1-alpha)*w(2,idcell)    + alpha*w(2,i)
c ... produtos interno         
          wfn   =      wf(1)*n(1,i) +  wf(2)*n(2,i)
c ... fluxo convectivo de primeira ordem
          cv    = rof*wfn*meta(i)
c ... fluxo convectivo de ordem superior
          if(cv.lt.0.0d0) then
c ... distancia normal a aresta    
           gpn  = grad(1,i)*(xm(1,i)-xc(1,i))
     .          + grad(2,i)*(xm(2,i)-xc(2,i))
           gpn  = fluxl(i) *gpn
         else  
           gpn = grad(1,idcell)*df(1,i)+grad(2,idcell)*df(2,i)
           gpn = fluxl(idcell) *gpn
          endif  
          cvc = gpn
c .....................................................................
          pl = pl+max(cv,0.0d0)*u(idcell)+min(cv,0.d0)*u(i)+cv*cvc
        endif
      enddo
c ......................................................................
      p = p - dt*(1.0d0-teta)*pl     
c ......................................................................

      return 
      end
c **********************************************************************
