c **********************************************************************
c * CELL_E_EB_GGL : Celula 2D com o termpo de correção para o fluxo    *
c * divusivo e convectivo para equacao da energia formulacao gas ideal *
c * (EDGE-BASE-TVD)                                                    *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * a      - nao definido                                              *
c * x      - coordenadas dos vertices da celula central e do seus      *
c * vizinhos                                                           *
c * u      - temperatura por celula                                    *
c * u0     - temperatura da celula do passo anterior                   *
c * grad   - gradiente reconstruido na celula                          *
c * fluxl  - limitador de fluxo na celulas                             *
c * k      - propriedades da celula e do seus vizinhos                 *
c * rm     - nao definido                                              *
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
c * p      - vetor de forcas                                           *
c * rm     - residou da equacao da quantidade de movmentos             *
c * grad   - gradiente (GREEN GAUSS LINEAR)                            *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine cell_e_eb_idg_ggl(a,x,u,u0,rO,w,grad,fluxl,k,rm,p,sedge
     .                            ,dt,pedge,viz,nshared,ndm,iws,nel,sn
     .                            ,acod,bs)
      implicit none
      include 'idealGas.fi'
      real*8 areacell,uf,gf(2),gfn,gfk,w(ndm,*),wf(ndm),wfn,cv,cvc
      real*8 beta,limit,du,rof,gpk,roc,vm,u0(5),m,cc,dt,df(3,4),rm
      real*8 fluxl(2),km(3,4),rO(3,nshared+1),gama            
      real*8 a(*),sp,p,xc(3,5),x(ndm,nshared,nshared+1),dfd,nk,um(4),cd
      real*8 sedge(*),ksi(3,4),eta(3,4),n(3,4),meta(4),mksi(4),area(5)
      real*8 kc,u(5),ca(4),alpha,k(10,*),xm(3,4),ap,grad(ndm,*)
      real*8 eps,peclet,ug(4),xcg(3,4),kf,cEsp,mKsiF(4),tRef,ZERO
      integer pedge(*),viz(*),viznel,iws,nel,tRo2,enthalpia
      integer ndm,nshared,i,j,icod,acod,idcell,sn(2,*),l
      logical bs
      parameter (ZERO  = 1.0d-14)
      parameter (eps = 1.0d-14)
c ... 
c ... 1 - UD - upwind de primeira ordem
c ... 2 - LUD - upwind de segunda ordem
c ... 3 - Van Leer - TVD
c ... 4 - Van Albada - TVD
c ... 5 - Mid - Mod - TVD
c ... 6 - MUSCL - TVD
c ... 7 - OSHER - TVD
      icod   = 1 
      idcell = nshared + 1
      tRo2    = 3
c .....................................................................
c
c ...      
      call cellGeom2D(x,mksi,meta,ksi,eta,xc,area,n,xm,df,ca
     .               ,viz,nshared,ndm,sn,acod)
c .....................................................................
c
c ...
      if(iws  .eq. 5) goto 500 
      if(iws  .eq. 6) goto 600 
c .....................................................................
c
c ... propriedeade da celula
      roc     = rO(tRo2,idCell)
      kc      =  k(3,idCell)
      cEsp    =  k(4,idCell)
      tRef    =  k(5,idCell)
      cc      =  k(6,idCell)
      icod    =  k(7,idCell)
c ... numero d peclet
      vm     = dsqrt(w(1,idcell)*w(1,idcell) + w(2,idcell)*w(2,idcell))
      peclet = cEsp*roc*vm*dsqrt(area(idcell))/kc
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
c
c ... considerando
          if(pedge(i) .eq. 2) then
            enthalpia = cEsp*(sedge(i)-tRef) 
c            enthalpia = cEsp*(sedge(i)+tConv)           
            ug(i) = enthalpia
          else 
            gama = 2.d0*ca(i)
            ug(i) = u(idCell) 
     .            + gama*(grad(1,1)*n(1,i)+grad(2,1)*n(2,i))
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
      call vectorKm2d(xcg,xc,x,xm,km,mKsiF,sn,nshared,ndm)
c .....................................................................      
c        
c ... 
      goto(100,200,300,400) iws
c .....................................................................
c
  100 continue
c ...
      p   = 0.d0
      sp  = 0.d0
      cd  = 0.d0
      cvc = 0.d0
      do i = 1, nshared
        viznel = viz(i)
c ... condicao de contorno      
        if( viznel .lt. 0) then
          a(i) = 0.d0
c ... fluxo convectivo na face de contorno          
          wfn = w(1,idcell)*n(1,i) +  w(2,idcell)*n(2,i)
          cv  = roc*wfn*meta(i)
c ... condicao de parede (fluxo prescrito)
          if(pedge(i) .eq. 0 ) then
            p   = p  + sedge(i)*meta(i)
c ... condicao de parede (temperatura prescrita)
          else if(pedge(i) .eq. 1 ) then
            ap = (kc/cEsp)*meta(i)/ca(i)
            sp = sp + ap
            enthalpia = cEsp*(sedge(i)-tRef) 
c            enthalpia = cEsp*(sedge(i)+tConv)      
            p  = p  + ap*enthalpia
c ... entrada ou saida de massa                      
          else if(pedge(i) .eq. 2) then
c ... fluxo advectivo de primroeira ordem
            if(dabs(wfn) .gt. ZERO) then  
              if(cv.lt.0.0d0) then
                cv  = k(2,idCell)*wfn*meta(i)
                enthalpia = cEsp*(sedge(i)-tRef) 
                p   = p  - cv*enthalpia
              else
               sp  = sp + cv
              endif
            endif  
c ... saida localmente parabolica                      
          elseif(pedge(i) .eq. 4) then
c ... fluxo divusivo zero
            sp  = sp + cv            
          endif
        else
c ... interpolacao da propriedaes
          alpha = mKsiF(i)/mksi(i)
c           kf    = (1.0d0-alpha)*(kc/cEsp) + alpha*(k(3,i)/k(4,i))
c ... media harmonica  
          kf    = alpha/(kc/cEsp) + (1.0d0-alpha)/(k(3,i)/k(4,i)) 
          kf    = 1.0d0/kf 
          wf(1) = (1.0d0-alpha)*w(1,idcell)    + alpha*w(1,i)
          wf(2) = (1.0d0-alpha)*w(2,idcell)    + alpha*w(2,i)
          gf(1) = (1.0d0-alpha)*grad(1,idcell) + alpha*grad(1,i)
          gf(2) = (1.0d0-alpha)*grad(2,idcell) + alpha*grad(2,i)
c ... produtos interno 
          nk    =   n(1,i)*ksi(1,i) + n(2,i)*ksi(2,i)
          wfn   =      wf(1)*n(1,i) +  wf(2)*n(2,i)
          gfk   =    gf(1)*ksi(1,i) +  gf(2)*ksi(2,i)
c ... interpolacao de gradiente (Darwish, F. Maukalled - 2003)
          gf(1) = gf(1) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(1,i)
          gf(2) = gf(2) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(2,i)
c ... produtos interno         
          gfn   =      gf(1)*n(1,i) +  gf(2)*n(2,i)
          gfk   =    gf(1)*ksi(1,i) +  gf(2)*ksi(2,i)      
c ... difusao direta        
          dfd   = (kf*meta(i))/(nk*mksi(i))
c ... fator de corecao
          cd    = kf*meta(i)*(gfn-gfk/nk)
          a(i)  = dfd
          p     = p   + cd
c ... fluxo convectivo de primeira ordem
c         if( wfn .gt. 0) then 
c           rof = rO(tRo2,idCell)
c         elseif( wfn .lt. 0 ) then
c           rof = rO(tRo2,i)
c         else           
c           rof = (1.0d0-alpha)*rO(tRo2,idCell) +  alpha*rO(tRo2,i)     
c         endif
          rof = (1.0d0-alpha)*rO(tRo2,idCell) +  alpha*rO(tRo2,i) 
          cv = rof*wfn*meta(i)
c ... fluxo convectivo de ordem superior
          if(cv.lt.0.0d0) then
            du  = u(idcell) - u(i) + eps 
            gpk = grad(1,i)*(xc(1,idcell)-xc(1,i))
     .          + grad(2,i)*(xc(2,idcell)-xc(2,i))
            alpha = mKsiF(i)/mksi(i)  
            beta     = 2.0d0*gpk/du - 1.0d0
            cvc   = (1.0d0-alpha)*limit(beta,icod)*du
          else  
            du  = u(i) - u(idcell) + eps 
            gpk = (grad(1,idcell)*ksi(1,i)
     .           + grad(2,idcell)*ksi(2,i))*mksi(i)
            alpha = mKsiF(i)/mksi(i)  
            beta     = 2.0d0*gpk/du - 1.0d0
            cvc   = alpha*limit(beta,icod)*du
          endif
c ... interpolacao unidirecional
          cvc   = cvc +  gf(1)*km(1,i) + gf(2)*km(2,i)
c .....................................................................
          a(i) = a(i) + max(-cv,0.0d0)
          sp   = sp   + cv
          p    =  p   - cv*cvc
c          write(3,'(i5,9es16.6)'),nel,dfd,cd,nk,cv,cvc,peclet
        endif
      enddo
c .....................................................................
c
c ...
      a(idcell) = sp
      do i = 1, nshared
        a(idcell) = a(idcell) + a(i)
      enddo
c .....................................................................
c
c ... matriz de massa 
      m    =  rO(tRo2,idCell)*cc*area(idCell)/dt
c ... Backward de segunda ordem 
      if(bs) m = 1.5d0*m
      a(idcell) = a(idCell) + m
c ... calculo do residuo R(i) = F - Ku(i) 
      rm = p - a(idCell)*u(idCell)
      do i = 1, nshared
        rm = rm + a(i)*u(i)
      enddo
c .....................................................................     
      return
c .....................................................................
c
c ... recontrucao de gradiente
  200 continue
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
c
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
      grad(1,1) = 0.0d0
      grad(2,1) = 0.0d0
      do i = 1, nshared
        uf = um(i)
        grad(1,1) = grad(1,1) + uf*n(1,i)
        grad(2,1) = grad(2,1) + uf*n(2,i)
      enddo
      grad(1,1) = grad(1,1) / area(1)
      grad(2,1) = grad(2,1) / area(1)
c ... funcao de limatacao do fluxo para o termo advectivo 
      fluxl(1) = 1.0d0
      fluxl(2) = 1.0d0
      return
c .....................................................................
c
c ... posprocessamento do campo de velocidades
  300 continue
      return
c .....................................................................
c
c ... F = fontes no volume + Vpu0/dt 
  400 continue
c .....................................................................
c
c ... Vp(n)u(n)
      m    = ro(2,idCell)*cc*area(idCell)/dt
      p    = m*u(idCell)
c ... Backward de segunda ordem
      if(bs) then
        m  = ro(1,idCell)*cc*area(idCell)/dt
        p  = 2*p - 0.5d0*m*u0(idCell) 
      endif
c .....................................................................
c
c ... VdtQ
      if(pedge(idCell) .ne. 0)  p = p + sedge(idCell)*area(idCell)
c .....................................................................
      return
c .....................................................................
c
c ... posprocessamento do campo de velocidades
  500 continue
      do i = 1, nshared
        viznel = viz(i)
c ... condicao de contorno      
        if( viznel .lt. 0) then
c ... condicao de parede          
          if(pedge(i) .eq. 0 
     .      .or. pedge(i) .eq. 1) then
            wfn = w(1,idCell)*eta(1,i) + w(2,idCell)*eta(2,i)
            w(1,idCell) = wfn*eta(1,i)
            w(2,idCell) = wfn*eta(2,i)
          endif
        endif
      enddo
      return
c .....................................................................
c
c ... (V/Temp)
  600 continue
      rm = area(idCell)/(rm+273.15d0)
c .....................................................................
c
c .....................................................................
      return
c .....................................................................
      end
c **********************************************************************
c
c **********************************************************************
c * CELL_T_EB_GGL : Celula 2D com o termpo de correção para o fluxo    *
c * divusivo e convectivo para equacao de transporte                   *
c * (EDGE-BASE-TVD)                                                    *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * a      - nao definido                                              *
c * x      - coordenadas dos vertices da celula central e do seus      *
c * vizinhos                                                           *
c * u      - temperatura por celula                                    *
c * u0     - temperatura da celula do passo anterior                   *
c * grad   - gradiente reconstruido na celula                          *
c * fluxl  - limitador de fluxo na celulas                             *
c * k      - propriedades da celula e do seus vizinhos                 *
c * rm     - nao definido                                              *
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
c * p      - vetor de forcas                                           *
c * rm     - residou da equacao da quantidade de movmentos             *
c * grad   - gradiente (GREEN GAUSS LINEAR)                            *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine cell_t_eb_ggl(a,x,u,u0,rO,w,grad,fluxl,k,rm,p,sedge
     .                        ,dt,pedge,viz,nshared,ndm,iws,nel,sn
     .                        ,acod,bs)
      implicit none
      real*8 areacell,uf,gf(2),gfn,gfk,w(ndm,*),wf(ndm),wfn,cv,cvc
      real*8 r,limit,du,rof,gpk,vm,u0(5),m,cc,dt,df(3,4),rm
      real*8 fluxl(nshared+1),km(3,4),rO(3,nshared+1)            
      real*8 a(*),sp,p,xc(3,5),x(ndm,nshared,nshared+1),dfd,nk,um(4),cd
      real*8 sedge(*),ksi(3,4),eta(3,4),n(3,4),meta(4),mksi(4),area(5)
      real*8 kc,u(5),ca(4),alpha,k(10,*),xm(3,4),ap,grad(ndm,*)
      real*8 eps,peclet,ug(4),xcg(3,4),kf,mKsiF(4),tRef
      integer pedge(*),viz(*),viznel,iws,nel,tRo2,specificMassA 
      integer ndm,nshared,i,j,icod,acod,idcell,sn(2,*),l
      logical bs
      parameter (eps = 1.0d-14)
c ... 
c ... 1 - UD - upwind de primeira ordem
c ... 2 - LUD - upwind de segunda ordem
c ... 3 - Van Leer - TVD
c ... 4 - Van Albada - TVD
c ... 5 - Mid - Mod - TVD
c ... 6 - MUSCL - TVD
c ... 7 - OSHER - TVD
      icod   = 1 
      idcell = nshared + 1
      tRo2    = 3
c .....................................................................
c
c ...      
      call cellGeom2D(x,mksi,meta,ksi,eta,xc,area,n,xm,df,ca
     .               ,viz,nshared,ndm,sn,acod)
c .....................................................................
c
c ... propriedeade da celula
      rof          = rO(tRo2,idCell)
      specificMassA =  k(2,idCell)
      kc            =  k(3,idCell)
c      cEsp          =  k(4,idCell)
      cc            =  k(6,idCell)
      icod          =  k(7,idCell)
c ... numero d peclet
      vm     = dsqrt(w(1,idcell)*w(1,idcell) + w(2,idcell)*w(2,idcell))
      peclet = specificMassA*vm*dsqrt(area(idcell))/kc
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
c
c ... considerando
          if(pedge(i) .eq. 2) then
            ug(i) = 2.d0*sedge(i)-u(idcell)
          else 
            ug(i) = u(idcell)  
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
      call vectorKm2d(xcg,xc,x,xm,km,mKsiF,sn,nshared,ndm)
c .....................................................................      
c        
c ... 
      goto(100,200,300,400) iws
c .....................................................................
c
  100 continue
c ...
      p   = 0.d0
      sp  = 0.d0
      cd  = 0.d0
      cvc = 0.d0
      do i = 1, nshared
        viznel = viz(i)
c ... condicao de contorno      
        if( viznel .lt. 0) then
          a(i) = 0.d0
c ... fluxo convectivo na face de contorno          
          wfn = w(1,idcell)*n(1,i) +  w(2,idcell)*n(2,i)
          cv  = rof*wfn*meta(i)
c ... condicao de parede (fluxo prescrito)
          if(pedge(i) .eq. 0 ) then
            p   = p  + sedge(i)*meta(i)
c ... condicao de parede (temperatura prescrita)
          else if(pedge(i) .eq. 1 ) then
            ap = kc*meta(i)/ca(i)
            sp = sp + ap
            p  = p  + ap*sedge(i)
c ... entrada ou saida de massa                      
          else if(pedge(i) .eq. 2) then
c ... fluxo advectivo de primeira ordem
            if(cv.lt.0.0d0) then
               cv  = specificMassA*wfn*meta(i)
               p   = p  - cv*sedge(i)
            else
               sp  = sp + cv
            endif
c ... saida localmente parabolica                      
          elseif(pedge(i) .eq. 4) then
c ... fluxo divusivo zero
            sp  = sp + cv            
          endif
        else
c ... interpolacao da propriedaes
          alpha = mKsiF(i)/mksi(i)
c ... media harmonica  
          kf    = alpha/kc + (1.0d0-alpha)/k(3,i) 
          kf    = 1.0d0/kf
          wf(1) = (1.0d0-alpha)*w(1,idcell)    + alpha*w(1,i)
          wf(2) = (1.0d0-alpha)*w(2,idcell)    + alpha*w(2,i)
          gf(1) = (1.0d0-alpha)*grad(1,idcell) + alpha*grad(1,i)
          gf(2) = (1.0d0-alpha)*grad(2,idcell) + alpha*grad(2,i)
c ... produtos interno 
          nk    =   n(1,i)*ksi(1,i) + n(2,i)*ksi(2,i)
          wfn   =      wf(1)*n(1,i) +  wf(2)*n(2,i)
          gfk   =    gf(1)*ksi(1,i) +  gf(2)*ksi(2,i)
c ... interpolacao de gradiente (Darwish, F. Maukalled - 2003)
          gf(1) = gf(1) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(1,i)
          gf(2) = gf(2) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(2,i)
c ... produtos interno         
          gfn   =      gf(1)*n(1,i) +  gf(2)*n(2,i)
          gfk   =    gf(1)*ksi(1,i) +  gf(2)*ksi(2,i)      
c ... difusao direta        
          dfd   = (kf*meta(i))/(nk*mksi(i))
c ... fator de corecao
          cd    = kf*meta(i)*(gfn-gfk/nk)
          a(i)  = dfd
          p     = p   + cd
c ... fluxo convectivo de primeira ordem
          if( wfn .gt. 0) then 
            rof = rO(tRo2,idCell)
          elseif( wfn .lt. 0 ) then
            rof = rO(tRo2,i)
          else           
            rof = (1.0d0-alpha)*rO(tRo2,idCell) +  alpha*rO(tRo2,i)     
          endif
c           rof = (1.0d0-alpha)*rO(tRo2,idCell) +  alpha*rO(tRo2,i)   
           cv  = rof*wfn*meta(i)
c ... fluxo convectivo de ordem superior
          if(cv.lt.0.0d0) then
            du  = u(idcell) - u(i) + eps 
            gpk = grad(1,i)*(xc(1,idcell)-xc(1,i))
     .          + grad(2,i)*(xc(2,idcell)-xc(2,i))
            alpha = mKsiF(i)/mksi(i)  
            r     = 2.0d0*gpk/du - 1.0d0
            cvc   = (1.0d0-alpha)*limit(r,icod)*du
          else  
            du  = u(i) - u(idcell) + eps 
            gpk = (grad(1,idcell)*ksi(1,i)
     .           + grad(2,idcell)*ksi(2,i))*mksi(i)
            alpha = mKsiF(i)/mksi(i)  
            r     = 2.0d0*gpk/du - 1.0d0
            cvc   = alpha*limit(r,icod)*du
          endif
c ... interpolacao unidirecional
          cvc   = cvc +  gf(1)*km(1,i) + gf(2)*km(2,i)
c .....................................................................
          a(i) = a(i) + max(-cv,0.0d0)
          sp   = sp   + cv
          p    =  p   - cv*cvc
c          write(3,'(i5,9es16.6)'),nel,dfd,cd,nk,cv,cvc,peclet
        endif
      enddo
c .....................................................................
c
c ...
      a(idcell) = sp
      do i = 1, nshared
        a(idcell) = a(idcell) + a(i)
      enddo
c .....................................................................
c
c ... matriz de massa 
      m    =  rO(tRo2,idCell)*cc*area(idCell)/dt
c ... Backward de segunda ordem 
      if(bs) m = 1.5d0*m
      a(idcell) = a(idCell) + m
c ... calculo do residuo R(i) = F - Ku(i) 
      rm = p - a(idCell)*u(idCell)
      do i = 1, nshared
        rm = rm + a(i)*u(i)
      enddo
c .....................................................................     
      return
c .....................................................................
c
c ... recontrucao de gradiente
  200 continue
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
c
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
      grad(1,1) = 0.0d0
      grad(2,1) = 0.0d0
      do i = 1, nshared
        uf = um(i)
        grad(1,1) = grad(1,1) + uf*n(1,i)
        grad(2,1) = grad(2,1) + uf*n(2,i)
      enddo
      grad(1,1) = grad(1,1) / area(1)
      grad(2,1) = grad(2,1) / area(1)
c ... funcao de limatacao do fluxo para o termo advectivo 
      fluxl(1) = 1.0d0
      return
c .....................................................................
c
c ... posprocessamento do campo de velocidades
  300 continue
      return
c .....................................................................
c
c ... F = fontes no volume + Vpu0/dt 
  400 continue
c .....................................................................
c
c ... Vp(n)u(n)
      m    = rO(2,idCell) *cc*area(idCell)/dt
      p    = m*u(idCell)
c ... Backward de segunda ordem
      if(bs) then
        m  = rO(1,idCell)*cc*area(idCell)/dt
        p  = 2*p - 0.5d0*m*u0(idCell) 
      endif
c .....................................................................
c
c ... VdtQ
      if(pedge(idCell) .ne. 0)  p = p + sedge(idCell)*area(idCell)
c .....................................................................
      return
c .....................................................................
      end
c **********************************************************************
c
c **********************************************************************
c * CELL_AD_EB_GGL : Celula quadrilatera de com o termpo de correção   *
c * para o fluxo divusivo e convectivo                                 *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * a      - nao definido                                              *
c * x      - coordenadas dos vertices da celula central e do seus      *
c * vizinhos                                                           *
c * u      - temperatura por celula                                    *
c * u0     - temperatura da celula do passo anterior                   *
c * rO     - massa especifica da celula e sua vizinhas em t e t+1      *
c * grad   - gradiente reconstruido na celula                          *
c * fluxl  - limitador de fluxo na celulas                             *
c * k      - propriedades da celula e do seus vizinhos                 *
c * rm     - nao definido                                              *
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
c * p      - vetor de forcas                                           *
c * rm     - residou da equacao da quantidade de movmentos             *
c * grad   - gradiente (GREEN GAUSS LINEAR)                            *
c * -------------------------------------------------------------------*
c **********************************************************************
c      subroutine cell_ad_eb_ggl(a,x,u,u0,w,grad,fluxl,k,sP,p,sedge,dt
c     .                         ,pedge,viz,nshared,ndm,iws,nel,sn
c     .                         ,teta,acod,bs)
c      implicit none
c      real*8 areacell,uf,gf(2),gfn,gfk,w(ndm,*),wf(ndm),wfn,cv,cvc
c      real*8 r,limit,du,rof,gpk,roc,vm,u0(5),m,cc,dt,df(3,4)
c      real*8 fluxl(nshared+1),km(3,4),rm,sP            
c      real*8 a(*),p,xc(3,5),x(ndm,nshared,nshared+1),dfd,nk,um(4),cd
c      real*8 sedge(*),ksi(3,4),eta(3,4),n(3,4),meta(4),mksi(4),area(5)
c      real*8 kc,u(5),ca(4),alpha,k(10,*),xm(3,4),ap,grad(ndm,*)
c      real*8 eps,peclet,cfl,ug(4),xcg(3,4),tdt,kf,teta,pl
c      integer pedge(*),viz(*),viznel,iws,nel
c      integer ndm,nshared,i,j,icod,acod,idcell,sn(2,*),l
c      logical bs
c      parameter (eps = 1.0d-8)
c ... 
c ... 1 - UD - upwind de primeira ordem
c ... 2 - LUD - upwind de segunda ordem
c ... 3 - Van Leer - TVD
c ... 4 - Van Albada - TVD
c ... 5 - Mid - Mod - TVD
c ... 6 - MUSCL - TVD
c ... 7 - OSHER - TVD
c      icod = 1 
c      idcell = nshared + 1
c .....................................................................
c
c ...      
c       call cellGeom2D(x,mksi,meta,ksi,eta,xc,area,n,xm,df,ca
c     .                ,viz,nshared,ndm,sn,acod)
c .....................................................................
c
c ... propriedeade da celula
c      kc      = k(1,nshared+1)
c      roc     = k(2,nshared+1)
c      cc      = k(3,nshared+1)
c      icod    = k(4,nshared+1)
c ... numero d peclet
c      vm     = dsqrt(w(1,idcell)*w(1,idcell) + w(2,idcell)*w(2,idcell))
c      peclet = roc*vm*dsqrt(area(idcell))/kc
c ... cfl
c      cfl = cc*dt*vm/dsqrt(area(idcell))
c .....................................................................
c
c ... criando celulas fantasma para o contorno  
c      do i = 1, nshared
c        viznel = viz(i)
c        ug(i)  = u(i)
c ... contorno      
c        if( viznel .lt. 0) then
c ... gerando o centroide da celula fantasmas        
c          do j = 1, ndm
c            xcg(j,i) = xc(j,idcell) + 2.d0*ca(i)*n(j,i)
c          enddo
c ... gerando a condicao de contorno da celula fantasma
c
c ... considerando o gradiente normal para estimar o potencial no centro
c     da celula fantasmas)
c          if(pedge(i) .eq. 0 ) then
c            ug(i) = u(idcell) + 2.d0*sedge(i)/kc*ca(i)
c ... temperatura prescrita (extrapolacao linear para a centro o
c     potencial da celula fantasmas)
c          else
c            ug(i) = 2.d0*sedge(i)-u(idcell)
c          endif
c        else
c          do j = 1, ndm
c            xcg(j,i) = xc(j,i)
c          enddo
c        endif
c      enddo
c .....................................................................
c
c ... vetor que une a intersecão da reta entre os dois centroides 
c     ao ponto médio da aresta
c      call vectorKm2d(xcg,xc,x,xm,km,sn,nshared,ndm)
c .....................................................................      
c        
c ... 
c      goto(100,200,300,400) iws
c .....................................................................
c
c  100 continue
c ...
c      p   = 0.d0
c      sp  = 0.d0
c      cd  = 0.d0
c      cvc = 0.d0
c      do i = 1, nshared
c        viznel = viz(i)
c ... condicao de contorno      
c        if( viznel .lt. 0) then
c          a(i) = 0.d0
c ... fluxo convectivo na face de contorno          
c          wfn = w(1,idcell)*n(1,i) +  w(2,idcell)*n(2,i)
c          cv  = roc*wfn*meta(i)
c ... fluxo prescrito 
c          if(pedge(i) .eq. 0 ) then
c ... fluxo difusivo          
c            p   =  p + meta(i)*sedge(i)
c ... fluxo advectivo
c            if(cv.gt.0.0d0)  then
c               p  = p  - sedge(i)/kc*ca(i)
c               sp = sp + cv
c            endif               
c          else
c ... fluxo difusivo     
c            ap  = kc*meta(i)/ca(i)
c            sp  = sp + ap
c            p   = p  + ap*sedge(i)
c ... fluxo advectivo de primeira ordem            
c            if(cv.lt.0.0d0) then
c              p   = p  - cv*sedge(i)
c            endif
c          endif
c        else
c ... interpolacao da propriedaes
c          alpha = ca(i)/mksi(i)
c           kf    = (1.0d0-alpha)*kc             + alpha*k(1,i)
c ... media harmonica  
c          kf   = alpha/kc + (1.0d0-alpha)/k(1,i) 
c          kf   = 1.0d0/kf 
c          rof   = (1.0d0-alpha)*roc            + alpha*k(2,i)
c          wf(1) = (1.0d0-alpha)*w(1,idcell)    + alpha*w(1,i)
c          wf(2) = (1.0d0-alpha)*w(2,idcell)    + alpha*w(2,i)
c          gf(1) = (1.0d0-alpha)*grad(1,idcell) + alpha*grad(1,i)
c          gf(2) = (1.0d0-alpha)*grad(2,idcell) + alpha*grad(2,i)
c          gfk   =    gf(1)*ksi(1,i) +  gf(2)*ksi(2,i)
c ... interpolacao de gradiente (Darwish, F. Maukalled - 2003)
c          gf(1) = gf(1) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(1,i)
c          gf(2) = gf(2) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(2,i)
c ... produtos interno         
c          nk    =   n(1,i)*ksi(1,i) + n(2,i)*ksi(2,i)
c          wfn   =      wf(1)*n(1,i) +  wf(2)*n(2,i)
c          gfn   =      gf(1)*n(1,i) +  gf(2)*n(2,i)
c ... difusao direta        
c          dfd   = (kf*meta(i))/(nk*mksi(i))
c ... fator de corecao
c          cd    = kf*meta(i)*(gfn-gfk/nk)
c          a(i)  = dfd   
c          p     = p   + cd
c ... fluxo convectivo de primeira ordem
c          cv = rof*wfn*meta(i)
c ... fluxo convectivo de ordem superior
c          if(cv.lt.0.0d0) then
c            du  = u(idcell) - u(i) + eps 
c            gpk = grad(1,i)*(xc(1,idcell)-xc(1,i))
c     .          + grad(2,i)*(xc(2,idcell)-xc(2,i))
c            alpha = ca(i)/mksi(i)  
c            r     = 2.0d0*gpk/du - 1.0d0
c            cvc   = (1.0d0-alpha)*limit(r,icod)*du
c          else  
c            du  = u(i) - u(idcell) + eps 
c            gpk = (grad(1,idcell)*ksi(1,i)
c     .           + grad(2,idcell)*ksi(2,i))*mksi(i)
c            alpha = ca(i)/mksi(i)  
c            r     = 2.0d0*gpk/du - 1.0d0
c            cvc   = alpha*limit(r,icod)*du
c          endif
c ... interpolacao unidirecional
c          cvc   = cvc +  gf(1)*km(1,i) + gf(2)*km(2,i)
c .....................................................................
c          a(i) = a(i) + max(-cv,0.0d0)
c          sp   = sp   + cv
c          p    =  p   - cv*cvc
c          write(3,'(i5,9es16.6)'),nel,dfd,cd,nk,cv,cvc,peclet
c        endif
c      enddo
c .....................................................................
c
c ...
c      tdt  = teta*dt
c      sp   = sp*tdt
c      p    =  p*tdt
c      a(idcell) = sp
c      do i = 1, nshared
c        a(i)     = a(i)      * tdt
c        a(idcell)= a(idcell) + a(i)
c      enddo  
c .....................................................................
c
c ... matriz de massa 
c      m    = roc*cc*area(idcell)
c ... Backward de segunda ordem 
c      if(bs) m = 1.5d0*m
c      a(idcell) = a(idcell) + m
c ... calculo do residuo R(i) = F - Ku(i) 
c      p = p - a(idcell)*u(idcell)
c      do i = 1, nshared
c        p = p + a(i)*u(i)
c      enddo
c .....................................................................     
c      return
c .....................................................................
c
c ... recontrucao de gradiente
c  200 continue
c ... vetor que une os centroides dos vizinho        
c ... vetor que une os centroides dos vizinho 
c      do i = 1, ndm
c        do l = 1, nshared
c ... centroide l  
c          eta(i,l) = xcg(i,sn(2,l)) - xcg(i,sn(1,l))    
c        enddo
c      enddo
c ... vetor normal com o modulo igual a arestas
c      do l = 1, nshared
c        n(1,l) = eta(2,l)
c        n(2,l) =-eta(1,l)
c      enddo
c .....................................................................
c ... area                 
c      area(1) = areacell(eta,acod)
c .....................................................................
c
c ...
c      do l = 1, nshared
c        um(l) = 0.5d0*(ug(sn(2,l))+ug(sn(1,l))) 
c      enddo    
c .....................................................................
c
c ... Reconstrucao linear Green-Gauss 
c      do i = 1, nshared
c        uf = um(i)
c        grad(1,1) = grad(1,1) + uf*n(1,i)
c        grad(2,1) = grad(2,1) + uf*n(2,i)
c      enddo
c      grad(1,1) = grad(1,1) / area(1)
c      grad(2,1) = grad(2,1) / area(1)
c ... funcao de limatacao do fluxo para o termo advectivo 
c      fluxl(1) = 1.0d0
c      return
c .....................................................................
c
c ...
c  300 continue
c      return
c .....................................................................
c
c ... F = fontes no volume + Vpu0 - dt*(1-teta)*F 
c  400 continue
c ...
c      p   = 0.d0
c      pl  = 0.d0
c      sp  = 0.d0
c      cd  = 0.d0
c      cvc = 0.d0
c      dfd = 0.d0
c .....................................................................
c
c ... Vp(n)u(n)
c      m    = roc*cc*area(idcell)
c      p    = m*u(idcell)
c ... Backward de segunda ordem 
c      if(bs) p = 2*p - 0.5d0*m*u0(idcell)
c .....................................................................
c
c ... VdtQ
c      if(pedge(idcell) .eq. 1)  p = p + sedge(idcell)*area(idcell)*dt
c .....................................................................
c
c ... Euler Bakward
c      if(teta .eq. 1.0d0) return
c .....................................................................
c
c ...
c      do i = 1, nshared
c        viznel = viz(i)
c        if(viznel .lt. 0) then
c ... fluxo convectivo na face de contorno          
c          wfn = w(1,idcell)*n(1,i) +  w(2,idcell)*n(2,i)
c          cv  = roc*wfn*meta(i)
c          if(pedge(i) .eq. 0) then
c           pl = pl - meta(i)*sedge(i)
c            pl = pl + cv*u(idcell) 
c          else
c            ap  = kc*meta(i)/ca(i)
c            pl  = pl  - ap*(sedge(i)-u(idcell))
c ... fluxo advectivo de primeira ordem            
c            if(cv.lt.0.0d0) then
c              pl = pl  + cv*sedge(i)
c            endif
c          endif
c        else 
c ... interpolacao da propriedaes
c          alpha = ca(i)/mksi(i)
c           kf    = (1.0d0-alpha)*kc             + alpha*k(1,i)
c ... media harmonica  
c          kf   = alpha/kc + (1.0d0-alpha)/k(i,1) 
c          kf   = 1.0d0/kf 
c          rof   = (1.0d0-alpha)*roc       + alpha*k(2,i)
c          wf(1) = (1.0d0-alpha)*w(1,idcell)    + alpha*w(1,i)
c          wf(2) = (1.0d0-alpha)*w(2,idcell)    + alpha*w(2,i)
c          gf(1) = (1.0d0-alpha)*grad(1,idcell) + alpha*grad(1,i)
c          gf(2) = (1.0d0-alpha)*grad(2,idcell) + alpha*grad(2,i)
c ... produtos interno         
c          nk    =   n(1,i)*ksi(1,i) + n(2,i)*ksi(2,i)
c          wfn   =      wf(1)*n(1,i) +  wf(2)*n(2,i)
c          gfn   =      gf(1)*n(1,i) +  gf(2)*n(2,i)
c          gfk   =    gf(1)*ksi(1,i) +  gf(2)*ksi(2,i)
c ... interpolacao de gradiente (Darwish, F. Maukalled - 2003)
c          gf(1) = gf(1) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(1,i)
c          gf(2) = gf(2) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(2,i)
c ... 
c          gf(1) = (1.0d0-alpha)*grad(1,idcell) + alpha*grad(1,i)
c          gf(2) = (1.0d0-alpha)*grad(2,idcell) + alpha*grad(2,i)
c ... difusao direta        
c          dfd   = (kf*meta(i))/(nk*mksi(i))
c ... fator de corecao
c          du    = u(i) - u(idcell)
c          cd    = kf*meta(i)*(gfn-gfk/nk)
c          pl    = pl - dfd*du + cd   
c ... fluxo convectivo de primeira ordem
c          cv    = rof*wfn*meta(i)
c ... fluxo convectivo de ordem superior
c          if(cv.lt.0.0d0) then
c            du    = u(idcell) - u(i) + eps            
c            gpk   = grad(1,i)*(xc(1,idcell)-xc(1,i))
c     .            + grad(2,i)*(xc(2,idcell)-xc(2,i))
c            alpha = ca(i)/mksi(i)      
c             r    = 2.0d0*gpk/du - 1.0d0
c             cvc  = (1.0d0-alpha)*limit(r,icod)*du
c          else  
c            du    = u(i) - u(idcell) + eps    
c            gpk   = (grad(1,idcell)*ksi(1,i)
c     .            +  grad(2,idcell)*ksi(2,i))*mksi(i)
c            alpha = ca(i)/mksi(i)      
c             r    = 2.0d0*gpk/du - 1.0d0
c             cvc  = alpha*limit(r,icod)*du
c          endif
c ... interpolacao unidirecional
c          cvc   = cvc +  gf(1)*km(1,i) + gf(2)*km(2,i)
c .....................................................................
c          pl = pl+max(cv,0.0d0)*u(idcell)+min(cv,0.d0)*u(i)+cv*cvc
c        endif
c      enddo
c ......................................................................
c      p = p - dt*(1.0d0-teta)*pl     
c ......................................................................
c      return
c      end
c **********************************************************************
c
c **********************************************************************
c * CELL_AD_VB : Celula 2D de com o termo de correção para o fluxo     *
c * divusivo e adveccao com funcao limitada por area                   *
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
c      subroutine cell_ad_vb_ggl(a,x,u,u0,w,grad,fluxl,k,sp,p,sedge
c     .                         ,dt,pedge,viz,nshared,ndm,iws
c     .                         ,nel,teta,sn,acod,bs)
c      implicit none
c      real*8 areacell,uf,gf(2),gfn,gfk,w(ndm,*),wf(ndm),wfn,cv,cvc
c      real*8 r,limitv,rof,gpk,roc,vm,u0(5),m,cc,dt,df(3,4),par(3)  
c      real*8 a(*),sp,p,xc(3,5),x(ndm,nshared,nshared+1)
c      real*8 cd,kf,dfd,nk,um(4)
c      real*8 sedge(nshared+1),ksi(3,4),eta(3,4),n(3,4),meta(4)
c      real*8 mksi(4),area(5)
c      real*8 kc,u(5),ca(4),alpha,k(10,*),xm(3,4),ap,grad(ndm,*)
c      real*8 eps,peclet,ug(4),xcg(3,4),fluxl(nshared+1),umax,umin
c      real*8 du,gpn,teta,tdt,pl,cfl
c      integer pedge(nshared+1),viz(nshared+1),viznel,iws,nel
c      integer ndm,nshared,i,j,l,icod,acod,idcell,sn(2,*)
c      logical bs
c      parameter (eps = 1.0d-8)
cc ... 
c ... 1 - UD - upwind de primeira ordem
c ... 2 - SOU - upwind de segunda ordem
c ... 3 - Barth         
c ... 4 - Van Albada - TVD
c      icod   = 1
c      idcell = nshared + 1
c .....................................................................
c
c ...      
c      call cellGeom2D(x,mksi,meta,ksi,eta,xc,area,n,xm,df,ca
c     .               ,viz,nshared,ndm,sn,acod)
c .....................................................................
c
c ... propriedeade da celula
c      kc      = k(1,nshared+1)
c      roc     = k(2,nshared+1)
c      cc      = k(3,nshared+1)
c      icod    = k(4,nshared+1)
c ... numero d peclet
c      vm    = dsqrt(w(1,idcell)*w(1,idcell) + w(2,idcell)*w(2,idcell))
c      peclet = roc*vm*dsqrt(area(idcell))/kc
c .....................................................................
c
c ... criando celulas fantasma para o contorno  
c      do i = 1, nshared
c        viznel = viz(i)
c        ug(i)  = u(i)
c ... contorno      
c        if( viznel .lt. 0) then
c ... gerando o centroide da celula fantasmas        
c          do j = 1, ndm
c            xcg(j,i) = xc(j,idcell) + 2.d0*ca(i)*n(j,i)
c          enddo
c ... gerando a condicao de contorno da celula fantasma
c ... temperatura prescrita (extrapolacao linear para a centro o
c     potencial da celula fantasmas)
c          if(pedge(i) .eq. 0 ) then
c            ug(i) = u(idcell) + 2.d0*sedge(i)/kc*ca(i)
c          else
c            ug(i) = 2.d0*sedge(i)-u(idcell)
c          endif  
c        else
c          do j = 1, ndm
c            xcg(j,i) = xc(j,i)
c          enddo
c        endif
c      enddo  
c .....................................................................
c
c ... 
c      goto(100,200,300,400) iws
c .....................................................................
c
c  100 continue
c ...
c      p   = 0.d0
c      sp  = 0.d0
c      cd  = 0.d0
c      cvc = 0.d0
c      do i = 1, nshared
c        viznel = viz(i)
c ... condicao de contorno      
c        if( viznel .lt. 0) then
c          a(i) = 0.d0
c ... fluxo convectivo na face de contorno          
c          wfn = w(1,idcell)*n(1,i) +  w(2,idcell)*n(2,i)
c          cv  = roc*wfn*meta(i)
c ... fluxo prescrito 
c          if(pedge(i) .eq. 0 ) then
c ... fluxo difusivo          
c            p   =  p + meta(i)*sedge(i)
c ... fluxo advectivo
c            if(cv.gt.0.0d0)  then
c               p  = p  - sedge(i)/kc*ca(i)
c               sp = sp + cv
c            endif               
c          else
c ... fluxo difusivo     
c            ap  = kc*meta(i)/ca(i)
c            sp  = sp + ap
c            p   = p  + ap*sedge(i)
c ... fluxo advectivo de primeira ordem            
c            if(cv.lt.0.0d0) then
c              p   = p  - cv*sedge(i)
c            endif
c          endif
c        else
c ... interpolacao da propriedades
c          alpha = ca(i)/mksi(i)
c ... media harmonica  
c          kf   = alpha/kc + (1.0d0-alpha)/k(1,i) 
c          kf   = 1.0d0/kf 
c          rof   = (1-alpha)*roc + alpha*k(2,i)
c          wf(1) = (1-alpha)*w(1,idcell) + alpha*w(1,i)
c          wf(2) = (1-alpha)*w(2,idcell) + alpha*w(2,i)
c          gf(1) = (1-alpha)*grad(1,idcell) + alpha*grad(1,i)
c          gf(2) = (1-alpha)*grad(2,idcell) + alpha*grad(2,i)
c ... produtos interno         
c          nk    =   n(1,i)*ksi(1,i) + n(2,i)*ksi(2,i)
c          wfn   =      wf(1)*n(1,i) +  wf(2)*n(2,i)
c          gfn   =      gf(1)*n(1,i) +  gf(2)*n(2,i)
c          gfk   =    gf(1)*ksi(1,i) +  gf(2)*ksi(2,i)
c ... interpolacao de gradiente (Darwish, F. Maukalled - 2003)
c          gf(1) = gf(1) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(1,i)
c          gf(2) = gf(2) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(2,i)
c ... produtos internos
c          gfn   =      gf(1)*n(1,i) +  gf(2)*n(2,i)
c          gfk   =    gf(1)*ksi(1,i) +  gf(2)*ksi(2,i)
c ... difusao direta        
c          dfd   = (kf*meta(i))/(nk*mksi(i))
c ... fator de corecao
c          cd    = kf*meta(i)*(gfn-gfk/nk)
c          a(i)  = dfd   
c          p     = p   + cd
c ... fluxo convectivo de primeira ordem
c          cv = rof*wfn*meta(i)
c ... fluxo convectivo de ordem superior
c          if(cv.lt.0.0d0) then
c ... distancia normal a aresta    
c            gpn  = grad(1,i)*(xm(1,i)-xc(1,i))
c     .           + grad(2,i)*(xm(2,i)-xc(2,i))
c            gpn  = fluxl(i) *gpn
c          else  
c            gpn = grad(1,idcell)*df(1,i)+grad(2,idcell)*df(2,i)
c            gpn = fluxl(idcell) *gpn
c          endif  
c          cvc = gpn
c .....................................................................
c          a(i) = a(i) + max(-cv,0.0d0)
c          sp   = sp   + cv
c          p    =  p   - cv*cvc
c          write(3,'(i5,9es16.6)'),nel,dfd,cd,nk,cv,cvc,peclet
c        endif
c      enddo
c .....................................................................
c
c ...
c      tdt  = teta*dt
c      sp   = sp*tdt
c      p    =  p*tdt
c      a(idcell) = sp
c      do i = 1, nshared
c        a(i)     = a(i)      * tdt
c        a(idcell)= a(idcell) + a(i)
c      enddo  
c .....................................................................
c
c ... matriz de massa 
c      m    = roc*cc*area(idcell)
c ... Backward de segunda ordem 
c      if(bs) m = 1.5d0*m
c      a(idcell) = a(idcell) + m
c ... calculo do residuo R(i) = F - Ku(i) 
c      p = p - a(idcell)*u(idcell)
c      do i = 1, nshared
c        p = p + a(i)*u(i)
c      enddo
cc .....................................................................     
c      return
c .....................................................................
c
c ... recontrucao de gradiente
c  200 continue
c .....................................................................      
c
c ... vetor que une os centroides dos vizinho 
c      do i = 1, ndm
c        do l = 1, nshared
c ... centroide l  
c          eta(i,l) = xcg(i,sn(2,l)) - xcg(i,sn(1,l))    
c        enddo
c      enddo
c ... vetor normal com o modulo igual a arestas
c      do l = 1, nshared
c        n(1,l) = eta(2,l)
c        n(2,l) =-eta(1,l)
c      enddo
c .....................................................................
c ... area                 
c      area(1) = areacell(eta,acod)
c .....................................................................
c
c ...
c      do l = 1, nshared
c        um(l) = 0.5d0*(ug(sn(2,l))+ug(sn(1,l))) 
c      enddo
c .....................................................................
c
c ... Reconstrucao linear Green-Gauss 
c      do i = 1, nshared
c        uf = um(i)
c        grad(1,1) = grad(1,1) + uf*n(1,i)
c        grad(2,1) = grad(2,1) + uf*n(2,i)
c      enddo
c      grad(1,1) = grad(1,1) / area(1)
c      grad(2,1) = grad(2,1) / area(1)
c ... funcao de limatacao do fluxo para o termo dvectivo 
c      fluxl(1) = 1.0d0
c ... obtendo o valor máximo da solucao
c      umax = u(idcell)
c      umin = umax
c      do i = 1, nshared
c        if(umax .lt. ug(i)) umax = ug(i)
c        if(umin .gt. ug(i)) umin = ug(i)  
c      enddo
c ... 
c      par(1) = area(idcell)
c      par(2) = 1.0d0
c      do i = 1, nshared
c        gpn = grad(1,1)*df(1,i) + grad(2,1)*df(2,i) + eps
c        uf = u(idcell) + gpn 
c        fluxl(2) = 1.0d0
c        if( uf .gt. u(idcell) ) then
c          r         = ( umax - u(idcell) ) / gpn
c          fluxl(2)  = limitv(r,par,icod)
c        elseif(uf.eq.u(idcell))then
c          fluxl(2)  = 1.0d0
c        elseif(uf.lt.u(idcell)) then
c          r = ( umin - u(idcell) ) / gpn
c          fluxl(2)  = limitv(r,par,icod)
c        endif
c ... menor parametro
c        fluxl(1) = min(fluxl(1),fluxl(2))
c      enddo
c ... fluxl >= 0
c      fluxl(1) = max(0.d0,fluxl(1))
c      return
c .....................................................................
c
c ...
c  300 continue
c      return
c ... F = fontes no volume + Vpu0 - dt*(1-teta)*F 
c  400 continue
c ...
c      p   = 0.d0
c      pl  = 0.d0
c      sp  = 0.d0
c      cd  = 0.d0
c      cvc = 0.d0
c      dfd = 0.d0
c .....................................................................
c
c ... Vp(n)u(n)
c      m    = roc*cc*area(idcell)
c      p    = m*u(idcell)
c ... Backward de segunda ordem 
c      if(bs) p = 2*p - 0.5d0*m*u0(idcell)
c .....................................................................
c
c ... VdtQ
c      if(pedge(idcell) .eq. 1)  p = p + sedge(idcell)*area(idcell)*dt
c .....................................................................
c
c ... Euler Bakward
c      if(teta .eq. 1.0d0) return
c .....................................................................
c
c ...
c      do i = 1, nshared
c        viznel = viz(i)
c        if(viznel .lt. 0) then
c ... fluxo convectivo na face de contorno          
c          wfn = w(1,idcell)*n(1,i) +  w(2,idcell)*n(2,i)
c          cv  = roc*wfn*meta(i)
c          if(pedge(i) .eq. 0) then
c            pl = pl - meta(i)*sedge(i)
c            pl = pl + cv*u(idcell) 
c          else
c            ap  = kc*meta(i)/ca(i)
c            pl  = pl  - ap*(sedge(i)-u(idcell))
c ... fluxo advectivo de primeira ordem            
c            if(cv.lt.0.0d0) then
c              pl = pl  + cv*sedge(i)
c            endif
c          endif
c        else 
c ... interpolacao da propriedades
c          alpha = ca(i)/mksi(i)
c           kf    = (1.0d0-alpha)*kc             + alpha*k(1,i)
c ... media harmonica  
c          kf    = alpha/kc + (1.0d0-alpha)/k(1,i) 
c          kf    = 1.0d0/kf 
c          rof   = (1-alpha)*roc            + alpha*k(2,i)
c          wf(1) = (1-alpha)*w(1,idcell)    + alpha*w(1,i)
c          wf(2) = (1-alpha)*w(2,idcell)    + alpha*w(2,i)
c          gf(1) = (1-alpha)*grad(1,idcell) + alpha*grad(1,i)
c          gf(2) = (1-alpha)*grad(2,idcell) + alpha*grad(2,i)
c ... produtos interno         
c          nk    =   n(1,i)*ksi(1,i) + n(2,i)*ksi(2,i)
c          wfn   =      wf(1)*n(1,i) +  wf(2)*n(2,i)
c          gfn   =      gf(1)*n(1,i) +  gf(2)*n(2,i)
c          gfk   =    gf(1)*ksi(1,i) +  gf(2)*ksi(2,i)
c ... interpolacao de gradiente (Darwish, F. Maukalled - 2003)
c          gf(1) = gf(1) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(1,i)
c          gf(2) = gf(2) + ((u(i) - u(idcell))/mksi(i)-gfk)*ksi(2,i)
c ...          
c          gfn   =      gf(1)*n(1,i) +  gf(2)*n(2,i)
c          gfk   =    gf(1)*ksi(1,i) +  gf(2)*ksi(2,i)
c ... difusao direta        
c          dfd   = (kf*meta(i))/(nk*mksi(i))
c ... fator de corecao
c          du    = u(i) - u(idcell)
c          cd    = kf*meta(i)*(gfn-gfk/nk)
c          pl    = pl - dfd*du + cd   
c ... fluxo convectivo de primeira ordem
c          cv    = rof*wfn*meta(i)
c ... fluxo convectivo de ordem superior
c          if(cv.lt.0.0d0) then
c ... distancia normal a aresta    
c           gpn  = grad(1,i)*(xm(1,i)-xc(1,i))
c     .          + grad(2,i)*(xm(2,i)-xc(2,i))
c           gpn  = fluxl(i) *gpn
c         else  
c           gpn = grad(1,idcell)*df(1,i)+grad(2,idcell)*df(2,i)
c           gpn = fluxl(idcell) *gpn
c          endif  
c          cvc = gpn
c .....................................................................
c          pl = pl+max(cv,0.0d0)*u(idcell)+min(cv,0.d0)*u(i)+cv*cvc
c        endif
c      enddo
c ......................................................................
c      p = p - dt*(1.0d0-teta)*pl     
c ......................................................................
c      return
c      end
c **********************************************************************
      
