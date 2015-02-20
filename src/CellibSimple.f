c **********************************************************************
c * CELL_SI_EB_IDG_GGL : Celula 2D com termo de correção               *
c * para o fluxo difusivo e convectivo para equação de momentun        *
c * formulacao gas ideal                                               *
c * (EDGE-BASE-TVD)                                                    *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * a      - nao definido                                              *
c * x      - coordenadas dos vertices da celula central e do seus      *
c * vizinhos                                                           *
c * u      - variavel da iteracao anterior                             *
c * u0     - variavel da iteracao anterior                             *
c * u1     - variavel da iteracao anterior                             *
c * rO     - massa especifica da celula e sua vizinhas em t -1, t e t+1*
c * w      - campo de velocidade estimado                              *
c * d      - nao definidp                                              *
c * grad1  - gradiente u(iws2) reconstruido da celula                  *
c * grad2  - gradiente u(2 ou 1) reconstruido da celula                *
c * gradP  - gradiente p  reconstruido da celula                       *
c * div    - divergente da velocidade                                  *
c * fluxl  - limitador de fluxo na celulas                             *
c * k      - propriedades da celula e do seus vizinhos                 *
c * rm     - nao definido                                              *
c * im     - nao definido                                              *
c *  p     - nao definido                                              *
c * sedge  - valores da condicoes de contorno por aresta               *
c * dt     - passo de tempo                                            *  
c * pedge  - tipo de condicao de contorno por aresta                   *
c * viz    - vizinhos da celula                                        *
c * nen    - numeros de nos por celula                                 *
c * nshared- numeros de faces por celulas                              *
c * iws1   - 1 - sistema de equacoes                                   *
c *        - 2 - reconstruacao de gradiente                            *
c * iws2   - 1 - direcao 1                                             *
c *        - 2 - direcao 2                                             *
c * nel    - numero da celula                                          *
c * sn(2,i)- nos da aresta i                                           *
c * acod   - codigo para o calculo da area                             *
c *         ( 3 - triangulo; 4 - quadrilatero)                         *
c * mP     - parametros da equacao de momentos(cfl,Reynalds)           *
c * bs     - Euller Backward de segunda ordem                          *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * a      - coeficientes                                              *
c * p      - vetor de forcas                                           *
c * d      - campo d                                                   *
c * rm(2)  - residou da equacao da quantidade de movmentos             *
c * im     - interploacao de momentos para velocidade nas faces        *
c * grad1   - gradiente (GREEN GAUSS LINEAR)                            *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine cell_si_eb_idg_ggl(a,x,u,u0,u1,rO,w,d,grad1,grad2,gradP
     .                             ,div,iM,fluxl,k,rm,p,sedge,dt,pedge
     .                             ,viz,nshared,ndm,iws1,iws2,nel,sn
     .                             ,acod,Mp,bs)
      implicit none
      include 'simple.fi'
      real*8 areacell,uf,gf(2),gfn,gfk,w(ndm,*),p,wf(ndm),wfn,cv,cvc
      real*8 r,limit,du,rof,gpk,vm,m,cc,dt,df(3,4),rO(3,*)
      real*8 fluxl(nshared+1),km(3,4),d(2,nshared+1),Mp(*)
      real*8 viscosity,specificMass,specificMassA,rm,iM(4,nshared+1) 
      real*8 a(*),sp,xc(3,5),x(ndm,nshared,nshared+1),dfd,nk,um(4),cd
      real*8 sedge(3,nshared+1),ksi(3,4),eta(3,4),n(3,4),meta(4)
      real*8 mksi(4),area(5),ap,u0(5),u1(5),kc
      real*8 u(5),ca(4),alpha,k(10,*),xm(3,4),grad1(ndm,*)
      real*8 eps,cfl,re,ug(4),xcg(3,4),kf,div(5),cEsp
      real*8 gradP(ndm,*),grad2(ndm,*)
      real*8 pf,p1,p2,pface,mKsiF(4),specificMassRef,tRef
      integer pedge(*),viz(*),viznel,iws1,iws2,nel,tRo1,tRo2
      integer ndm,nshared,i,j,icod,acod,idcell,sn(2,*),l
      real*8 const,aNb,ZERO
      logical bs
      parameter (ZERO  = 1.0d-7)
      parameter (eps   = 1.0d-14)
      parameter (const = 1.0d60)
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
      tRo1   = 2
      tRo2   = 3
c .....................................................................
c
c ...      
       call cellGeom2D(x,mksi,meta,ksi,eta,xc,area,n,xm,df,ca
     .                ,viz,nshared,ndm,sn,acod)
c .....................................................................
c
c ... propriedeade da celula
      viscosity    =  k(1,idCell)
      specificMass = rO(tRo2,idCell)
      specificMassA=  k(2,idCell)
      kc           =  k(3,idCell)
      cEsp         =  k(4,idCell)
      tRef         =  k(5,idCell)
      cc           =  k(6,idCell)
      icod         =  k(7,idCell)
      vm     = dsqrt(w(1,idcell)*w(1,idcell) + w(2,idcell)*w(2,idcell))
      if( iws1 .eq. 7) goto 700
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
c ... considerando o gradiente normal para estimar o potencial no centro
c     da celula fantasmas)
          if(pedge(i) .eq. 0 ) then
            ug(i) = -u(idCell)
c ... variavel prescrita prescrita (extrapolacao de primeira ordem para a centro o
c     da celula fantasmas)
          else if(pedge(i) .eq. 1 .or. pedge(i) .eq. 2 ) then
            ug(i) = sedge(iws2,i)
c ... condicao de contorno para pressao / locamente parabolica            
          else if(pedge(i) .eq. 3 .or. pedge(i) .eq. 4 ) then
c            gama = 2.d0*ca(i)
            ug(i) = u(idCell) 
c    .            + gama*(grad1(1,1)*n(1,i)+grad1(2,1)*n(2,i))
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
      goto(100,200,300,400,500,600) iws1
c .....................................................................
c
  100 continue
c ...
      p                = 0.d0
      pf               = 0.d0
      sp               = 0.d0
      cd               = 0.d0
      cvc              = 0.d0
      specificMassRef  = 0.d0
      do i = 1, nshared
        viznel = viz(i)
c ... condicao de contorno      
        if( viznel .lt. 0) then
          a(i) = 0.d0
          ap   = viscosity*meta(i)
c ... termos viscosos
          if(iws2 .eq. 1) then
            p = p + ap*(grad1(iws2,idCell)*n(1,i) 
     .            + grad2(iws2,idCell)*n(2,i))
          else
            p = p + ap*(grad1(iws2,idCell)*n(2,i) 
     .            + grad2(iws2,idCell)*n(1,i))
          endif
          p  = p -(2.d0/3.d0)*ap*div(idCell)*n(iws2,i)
c .....................................................................
c
c ... gradiente da pressao com resconstrucao de segunda ordem
c          pFace   = u0(idCell)+(gradP(1,idCell) + gradP(2,idCell))*ca(i)
c ... pressao prescrita
          if(pedge(i) .eq. 3) then  
             pFace   = sedge(3,i)
          else
            pFace   = u0(idCell)
            if(sPressure) then
              pFace = pFace + (gradP(1,idCell)*df(1,i) 
     .                      +  gradP(2,idCell)*df(2,i))
            endif
           endif
           pf      = pf + pFace*meta(i)*n(iws2,i)
c .....................................................................
c
c ... condicao de parede
c          re = (vm*specificMass*dsqrt(area(idcell)))/viscosity
c ... Kempf-Karman (1951)
c          cv   = 0.055*re**(-0.182d0)
c          cv   = dsqrt(0.5d0*cv)
c     .         *(w(1,idCell)*w(1,idCell)+ w(2,idCell)*w(2,idCell))
c          yPlus = (specificMass/viscosity)*(ca(i)*cv)
c .....................................................................

c ... parade  
          if(pedge(i) .eq. 0 ) then
            ap = viscosity*meta(i)/ca(i)
            sp = sp + ap
c ... parade movel 
          else if(pedge(i) .eq. 1) then
            ap = viscosity*meta(i)/ca(i)
            sp = sp + ap
            p  = p  + ap*sedge(iws2,i) 
c ... velocidade prescrita
          else if(pedge(i) .eq. 2 ) then
c ... fluxo advectivo de primeira ordem
            wfn = sedge(1,i)*n(1,i) + sedge(2,i)*n(2,i)
            cv  = specificMassA*wfn*meta(i)
            p   =  p  - cv*sedge(iws2,i)
c ... pressao prescrita
          else if(pedge(i) .eq. 3) then
            wfn = w(1,idCell)*n(1,i) + w(2,idCell)*n(2,i) 
            if(dabs(wfn) .gt. ZERO) then
              if( wfn .gt. 0.d0 ) then 
                cv  = specificMass*wfn*meta(i)
                sp  = sp + cv
              else if(wfn .lt. 0.d0) then
c                wfn = -dsqrt(2.d0*(sedge(3,i)-u0(idCell))/specificMassA)
                cv  =  specificMassA*wfn*meta(i)
                p   =  p + cv*wfn
              endif
          endif
c ... localmente parabolica
          else if(pedge(i) .eq. 4) then
            wfn = w(1,idcell)*n(1,i) + w(2,idcell)*n(2,i)
            cv  = specificMass*wfn*meta(i)
            sp  = sp + cv
          endif
c .....................................................................
c
c ...
        else
c ... interpolacao da propriedaes
          alpha = mKsiF(i)/mksi(i)
c ... media harmonica  
          kf    = alpha/viscosity + (1.0d0-alpha)/k(1,i) 
          kf    = 1.0d0/kf
c .....................................................................
          wf(1) = (1.0d0-alpha)*w(1,idcell)    + alpha*w(1,i)
          wf(2) = (1.0d0-alpha)*w(2,idcell)    + alpha*w(2,i)
          gf(1) = (1.0d0-alpha)*grad1(1,idcell) + alpha*grad1(1,i)
          gf(2) = (1.0d0-alpha)*grad1(2,idcell) + alpha*grad1(2,i)
          gfk   = gf(1)*ksi(1,i) +  gf(2)*ksi(2,i)
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
c            rof   = rO(tRo2,idCell)
c          elseif( wfn .lt. 0 ) then
c            rof   = rO(tRo2,i)
c          else           
c            rof   = (1.0d0-alpha)*rO(tRo2,idCell)  + alpha*rO(tRo2,i)
c          endif
          rof   = (1.0d0-alpha)*rO(tRo2,idCell)  + alpha*rO(tRo2,i)
          cv    = rof*wfn*meta(i)
c ... fluxo convectivo de ordem superior
          if(cv.lt.0.0d0) then
            du  = u(idcell) - u(i) + eps 
            gpk = grad1(1,i)*(xc(1,idcell)-xc(1,i))
     .          + grad1(2,i)*(xc(2,idcell)-xc(2,i))
            alpha = mKsiF(i)/mksi(i)
            r     = 2.0d0*gpk/du - 1.0d0
            cvc   = (1.0d0-alpha)*limit(r,icod)*du
          else  
            du  = u(i) - u(idcell) + eps 
            gpk = (grad1(1,idCell)*ksi(1,i)
     .           + grad1(2,idCell)*ksi(2,i))*mksi(i)
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
c .....................................................................
c
c ... gradiente da pressao com resconstrucao de segunda ordem
          alpha = mKsiF(i)/mksi(i)
          if(sPressure) then
            p1   = u0(idCell) + gradP(1,idCell)*df(1,i) 
     .                        + gradP(2,idCell)*df(2,i)
            p2   = u0(i)      + gradP(1,i)*(xm(1,i)-xc(1,i))
     .                        + gradP(2,i)*(xm(2,i)-xc(2,i))
            pFace = 0.5d0*(p1 + p2)
          else
            p1    = u0(idCell) 
            p2    = u0(i)  
            pFace = (1-alpha)*p1 + alpha*p2
          endif
          pf   = pf + pFAce*meta(i)*n(iws2,i)
c .....................................................................
c
c ... termos nao lineares da viscosidade
          
          alpha = mKsiF(i)/mksi(i)
          gf(1) = (1.0d0-alpha)*grad1(iws2,idcell) + alpha*grad1(iws2,i)
          gf(2) = (1.0d0-alpha)*grad2(iws2,idcell) 
     .          +        alpha *grad2(iws2,i)
          gpk   = (1.0d0-alpha)*div(idCell)  + alpha*div(i)
          if( iws2 .eq. 1) then
            p     = p + kf*meta(i)*(gf(1)*n(1,i)+ gf(2)*n(2,i))
          else if(iws2 .eq. 2 ) then
            p     = p + kf*meta(i)*(gf(1)*n(2,i) + gf(2)*n(1,i))
          endif
c ... divergente
          p     = p -(2.d0/3.d0)*kf*gpk*meta(i)*n(iws2,i)
c ... massa esoecifica dos vizinhos
        endif
      enddo
c .....................................................................
c
c ...
c      p         = p*dt
c      sp        = sp*dt 
      a(idcell) = sp
      do i = 1, nshared
c             a(i) = a(i)*dt
         a(idcell)= a(idcell) + a(i)
      enddo
c ......................................................................
c
c ... derivada temporal
      m    = ro(tRo2,idCell)*cc*area(idcell)/dt
c ... Backward de segunda ordem 
      if(bs) m = 1.5d0*m
      a(idCell) = a(idCell) + m
c .....................................................................
c
c ... under-relaxation(simple)
      a(idCell) = a(idCell)/underU 
      p         = p + (1-underU)*a(idCell)*u(idCell)
c .....................................................................
c
c ... Gravidade
       specificMass    = ro(tRo2,idCell)
c ... massa especifica de referencia (Empuxo= specificMass - specificMassRef) 
       specificMassRef = specificMassA
c ...       
       p   = p + (specificMass-specificMassRef)*g(iws2)*area(idCell)
c      p   = p + (specificMass-specificMassA)*g(iws2)*area(idCell)*dt
c      p   = p + specificMass*g(iws2)*area(idCell)
c .....................................................................
c
c ...
      rm  = 0.0d0
      aNb = 0.0d0
      do i = 1, nshared
        rm =  rm + a(i)*u(i)
        aNb= aNb + a(i)
      enddo
c ... interpolacao
      iM(iws2  ,idCell)  = rm + p
      iM(iws2+2,idCell)  = a(idCell)
c ... campo D para velocidade de correcao
      if(simpleC) then               
        d(iws2,idCell)     = area(idCell) / (a(idCell) -aNb) 
      else
        d(iws2,idCell)     = area(idCell)/ a(idCell)
      endif
c .....................................................................
c
c ... gradP
      p   = p - pf 
c      p   = p - gradP(iws2,idCell)*area(idCell)
c .....................................................................
c
c ... residou da celula
      rm  = p + rm - a(idCell)*u(idCell)
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
      grad1(1,1) = 0.0d0
      grad1(2,1) = 0.0d0
      do i = 1, nshared
        uf = um(i)
        grad1(1,1) = grad1(1,1) + uf*n(1,i)
        grad1(2,1) = grad1(2,1) + uf*n(2,i)
      enddo
      grad1(1,1) = grad1(1,1) / area(1)
      grad1(2,1) = grad1(2,1) / area(1)
c ... funcao de limatacao do fluxo para o termo advectivo 
      fluxl(1) = 1.0d0
      return
c .....................................................................
c
c
c .. F = dt*u0*V
  300 continue
      m    = ro(2,idCell)*cc*area(idCell)/dt
      p    = m*u(idCell)
c ... Backward de segunda ordem
      m    = ro(1,idCell)*cc*area(idcell)/dt 
      if(bs) p = 2*p - 0.5d0*m*u0(idcell)
      return  
c .....................................................................
c
c ... VdtQ
c      if(pedge(idcell) .eq. 1)  p = p + sedge(idcell)*area(idcell)
c .....................................................................
c
c ...
  400 continue
      return
  500 continue
      return
  600 continue
      return
c .....................................................................
c .. CFL e Reynalds
  700 continue
c ... cfl number
      cfl = (dt*vm)/dsqrt(area(idcell))
      Mp(1) = cfl 
c .....................................................................
c
c ... Reynolds number
      re = (specificMass*vm*dsqrt(area(idcell)))/viscosity
      Mp(2) = re
c .....................................................................
c
c ...
      Mp(3) = area(idCell)
c .....................................................................
c
c ...
      Mp(4) = vm          
c .....................................................................
c
c ... Prandtl number
      Mp(5) = (cEsp*viscosity)/kc
c .....................................................................
c
c ... Grashof number
      cd = dsqrt(g(1)*g(1) + g(2)*g(2))
      cd = cd*(dabs(ro(3,idCell) - specificMassA))*specificMassA
      cd = cd * area(idCell) * area(idCell) * area(idCell)
      cd = (cd*ro(3,idCell)*ro(3,idCell))/(viscosity*viscosity)
      Mp(6) = cd
c ......................................................................
c
c ... massa 
      Mp(7) = ro(3,idCell)*area(idCell)
c ......................................................................
c 
c ...
      p   = 0.0d0
      dfd = 0.0d0
      do i = 1, nshared
        viznel = viz(i)
c ... contorno      
        if( viznel .lt. 0) then
c ... velocidade prescrita
          if(pedge(i) .eq. 2 ) then
c ... fluxo advectivo de primeira ordem
            wfn = sedge(1,i)*n(1,i) + sedge(2,i)*n(2,i)
            cv  = specificMassA*wfn*meta(i)
            dfd = dfd + cv
c ... pressao prescrita
          else if(pedge(i) .eq. 3) then
            wfn = w(1,idCell)*n(1,i) + w(2,idCell)*n(2,i) 
            if(dabs(wfn) .gt. ZERO) then
              if( wfn .gt. 0.d0 ) then 
                cv  = specificMass*wfn*meta(i)
                 p  = p + cv
              else if(wfn .lt. 0.d0) then
c              wfn = -dsqrt(2.d0*(sedge(3,i) - u0(idCell))/specificMassA)
                cv  =  specificMassA*wfn*meta(i)
                dfd =  dfd + cv
              endif
            endif
          endif 
        endif 
      enddo
c ......................................................................
c
c ...
      Mp(8) = dfd     ! entrada de massa 
      Mp(9) = p       ! saida de massa
c ......................................................................
c
c ... Energia cinetica 
      Mp(10) = 0.5d0*vm*vm*area(idCell)
c ......................................................................
      return  
      end
c **********************************************************************
c
c **********************************************************************
c * CELL_SI_EB_I_GGL : Celula 2D com termo de correção                 *
c * para o fluxo difusivo e convectivo para equação de momentun        *
c * formulacao incomepressivel                                         *
c * (EDGE-BASE-TVD)                                                    *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * a      - nao definido                                              *
c * x      - coordenadas dos vertices da celula central e do seus      *
c * vizinhos                                                           *
c * u      - variavel da iteracao anterior                             *
c * u0     - variavel da iteracao anterior                             *
c * u1     - variavel da iteracao anterior                             *
c * rO     - massa especifica da celula e sua vizinhas em t -1, t e t+1*
c * w      - campo de velocidade estimado                              *
c * d      - nao definidp                                              *
c * grad1  - gradiente u(iws2) reconstruido da celula                  *
c * grad2  - gradiente u(2 ou 1) reconstruido da celula                *
c * gradP  - gradiente p  reconstruido da celula                       *
c * fluxl  - limitador de fluxo na celulas                             *
c * k      - propriedades da celula e do seus vizinhos                 *
c * rm     - nao definido                                              *
c * im     - nao definido                                              *
c *  p     - nao definido                                              *
c * sedge  - valores da condicoes de contorno por aresta               *
c * dt     - passo de tempo                                            *  
c * pedge  - tipo de condicao de contorno por aresta                   *
c * viz    - vizinhos da celula                                        *
c * nen    - numeros de nos por celula                                 *
c * nshared- numeros de faces por celulas                              *
c * iws1   - 1 - sistema de equacoes                                   *
c *        - 2 - reconstruacao de gradiente                            *
c * iws2   - 1 - direcao 1                                             *
c *        - 2 - direcao 2                                             *
c * nel    - numero da celula                                          *
c * sn(2,i)- nos da aresta i                                           *
c * acod   - codigo para o calculo da area                             *
c *         ( 3 - triangulo; 4 - quadrilatero)                         *
c * mP     - parametros da equacao de momentos(cfl,Reynalds)           *
c * bs     - Euller Backward de segunda ordem                          *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * a      - coeficientes                                              *
c * p      - vetor de forcas                                           *
c * d      - campo d                                                   *
c * rm(2)  - residou da equacao da quantidade de movmentos             *
c * im     - interploacao de momentos para velocidade nas faces        *
c * grad1   - gradiente (GREEN GAUSS LINEAR)                            *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine cell_si_eb_i_ggl(a,x,u,u0,u1,w,d,grad1,grad2,gradP
     .                           ,iM,fluxl,k,rm,p,sedge,dt,pedge
     .                           ,viz,nshared,ndm,iws1,iws2,nel,sn
     .                           ,acod,Mp,bs)
      implicit none
      include 'simple.fi'
      real*8 areacell,uf,gf(2),gfn,gfk,w(ndm,*),p,wf(ndm),wfn,cv,cvc
      real*8 r,limit,du,rof,gpk,vm,m,cc,dt,df(3,4)
      real*8 fluxl(nshared+1),km(3,4),d(2,nshared+1),Mp(*)
      real*8 viscosity,specificMassA,rm,iM(4,nshared+1) 
      real*8 a(*),sp,xc(3,5),x(ndm,nshared,nshared+1),dfd,nk,um(4),cd
      real*8 sedge(3,nshared+1),ksi(3,4),eta(3,4),n(3,4),meta(4)
      real*8 mksi(4),area(5),ap,u0(5),u1(5)
      real*8 u(5),ca(4),alpha,k(10,*),xm(3,4),grad1(ndm,*)
      real*8 eps,cfl,re,ug(4),xcg(3,4),kf
      real*8 gradP(ndm,*),grad2(ndm,*)
      real*8 pf,p1,p2,pface,mKsiF(4),aNb
      integer pedge(*),viz(*),viznel,iws1,iws2,nel
      integer ndm,nshared,i,j,icod,acod,idcell,sn(2,*),l
      real*8 const,ZERO
      real*8 s11,s22,s12,modS
      logical bs
      parameter (ZERO =  1.0d-60)
      parameter (eps   = 1.0d-14)
      parameter (const = 1.0d60)
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
c .....................................................................
c
c ...      
       call cellGeom2D(x,mksi,meta,ksi,eta,xc,area,n,xm,df,ca
     .                ,viz,nshared,ndm,sn,acod)
c .....................................................................
c
c ... propriedeade da celula
      viscosity    =  k(1,idCell)
      specificMassA=  k(2,idCell)
      cc           =  k(6,idCell)
      icod         =  k(7,idCell)
      vm     = dsqrt(w(1,idCell)*w(1,idCell) + w(2,idCell)*w(2,idCell))
      if( iws1 .eq. 7) goto 700
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
c ... considerando o gradiente normal para estimar o potencial no centro
c     da celula fantasmas)
          if(pedge(i) .eq. 0 ) then
           ug(i) = -u(idCell)
c ... variavel prescrita prescrita (extrapolacao de primeira ordem para a centro o
c     da celula fantasmas)
          else if(pedge(i) .eq. 1 .or. pedge(i) .eq. 2) then
            ug(i) = sedge(iws2,i)
c ... pressao prescrita            
          else if(pedge(i) .eq. 3 .or. pedge(i) .eq. 4) then
            ug(i) = u(idCell)          
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
      goto(100,200,300,400,500,600) iws1
c .....................................................................
c
  100 continue
c ...
      p            = 0.d0
      pf           = 0.d0
      sp           = 0.d0
      cd           = 0.d0
      cvc          = 0.d0
      do i = 1, nshared
        viznel = viz(i)
c ... condicao de contorno      
        if( viznel .lt. 0) then
          a(i) = 0.d0
          ap   = viscosity*meta(i)
c ... termos viscosos
          if(iws2 .eq. 1) then
            p = p + ap*(grad1(iws2,idCell)*n(1,i) 
     .            + grad2(iws2,idCell)*n(2,i))
          else
            p = p + ap*(grad1(iws2,idCell)*n(2,i) 
     .            + grad2(iws2,idCell)*n(1,i))
          endif
c .....................................................................
c
c ... gradiente da pressao com resconstrucao de segunda ordem
c ... pressao prescrita
         if(pedge(i) .eq. 3) then  
            pFace   = sedge(3,i)
          else
            pFace   = u0(idCell)
c ... calculo da pressao na face atraves de reconstrucao
            if(sPressure) then
              pFace = pFace + (gradP(1,idCell)*df(1,i) 
     .                      +  gradP(2,idCell)*df(2,i))
            endif   
          endif
          pf      = pf + pFace*meta(i)*n(iws2,i)
c .....................................................................
c
c ... condicao de parede
c          re = (vm*specificMass*dsqrt(area(idcell)))/viscosity
c ... Kempf-Karman (1951)
c          cv   = 0.055*re**(-0.182d0)
c          cv   = dsqrt(0.5d0*cv)
c     .         *(w(1,idCell)*w(1,idCell)+ w(2,idCell)*w(2,idCell))
c          yPlus = (specificMass/viscosity)*(ca(i)*cv)
c .....................................................................

c ... parade  
          if(pedge(i) .eq. 0 ) then
            ap = viscosity*meta(i)/ca(i)
            sp = sp + ap
c ... parade movel 
          else if(pedge(i) .eq. 1) then
            ap = viscosity*meta(i)/ca(i)
            sp = sp + ap
            p  = p  + ap*sedge(iws2,i)
c ... velocidade prescrita
          else if(pedge(i) .eq. 2 ) then
c ... fluxo advectivo de primeira ordem
            wfn = sedge(1,i)*n(1,i) + sedge(2,i)*n(2,i)
            cv  = specificMassA*wfn*meta(i)
            p   =  p  - cv*sedge(iws2,i)
c ... pressao prescrita
          else if(pedge(i) .eq. 3) then
            wfn = w(1,idCell)*n(1,i) + w(2,idCell)*n(2,i)
            if(dabs(wfn) .gt. ZERO) then
             if( wfn .gt. 0.d0 ) then 
               cv  = specificMassA*wfn*meta(i)
                sp  = sp + cv
              else if(wfn .lt. 0.0) then
c                wfn = -dsqrt(2.d0*(sedge(3,i)-u0(idCell))/specificMassA)
                cv  =  specificMassA*wfn*meta(i)
                p   =  p + cv*wfn
              endif
            endif
c ... localmente parabolica
          else if(pedge(i) .eq. 4) then
            wfn = w(1,idcell)*n(1,i) + w(2,idcell)*n(2,i)
            cv  = specificMassA*wfn*meta(i)
            sp  = sp + cv
          endif
c .....................................................................
c
c ...
        else
c ... interpolacao da propriedaes
          alpha = mKsiF(i)/mksi(i)
c ... media harmonica  
          kf    = alpha/viscosity + (1.0d0-alpha)/k(1,i) 
          kf    = 1.0d0/kf
c .....................................................................
          wf(1) = (1.0d0-alpha)*w(1,idcell)    + alpha*w(1,i)
          wf(2) = (1.0d0-alpha)*w(2,idcell)    + alpha*w(2,i)
          gf(1) = (1.0d0-alpha)*grad1(1,idcell) + alpha*grad1(1,i)
          gf(2) = (1.0d0-alpha)*grad1(2,idcell) + alpha*grad1(2,i)
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
          rof   = k(2,idCell)
          cv    = rof*wfn*meta(i)
c ... fluxo convectivo de ordem superior
          if(cv.lt.0.0d0) then
            du  = u(idcell) - u(i) + eps 
            gpk = grad1(1,i)*(xc(1,idcell)-xc(1,i))
     .          + grad1(2,i)*(xc(2,idcell)-xc(2,i))
            alpha = mKsiF(i)/mksi(i)
            r     = 2.0d0*gpk/du - 1.0d0
            cvc   = (1.0d0-alpha)*limit(r,icod)*du
          else  
            du  = u(i) - u(idcell) + eps 
            gpk = (grad1(1,idCell)*ksi(1,i)
     .           + grad1(2,idCell)*ksi(2,i))*mksi(i)
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
c .....................................................................
c
c ... gradiente da pressao com resconstrucao de segunda ordem
          alpha = mKsiF(i)/mksi(i)
          if(sPressure) then
            p1   = u0(idCell) + gradP(1,idCell)*df(1,i) 
     .                        + gradP(2,idCell)*df(2,i)
            p2   = u0(i)      + gradP(1,i)*(xm(1,i)-xc(1,i))
     .                        + gradP(2,i)*(xm(2,i)-xc(2,i))
            pFace = 0.5d0*(p1 + p2)
          else
            p1    = u0(idCell) 
            p2    = u0(i)  
            pFace = (1-alpha)*p1 + alpha*p2
          endif
          pf   = pf + pFAce*meta(i)*n(iws2,i)
c .....................................................................
c
c ... termos nao lineares da viscosidade
          alpha = mKsiF(i)/mksi(i)
          gf(1) = (1.0d0-alpha)*grad1(iws2,idcell) + alpha*grad1(iws2,i)
          gf(2) = (1.0d0-alpha)*grad2(iws2,idcell) 
     .          +        alpha *grad2(iws2,i)
          if( iws2 .eq. 1) then
            p     = p + kf*meta(i)*(gf(1)*n(1,i)+ gf(2)*n(2,i))
          else if(iws2 .eq. 2 ) then
            p     = p + kf*meta(i)*(gf(1)*n(2,i) + gf(2)*n(1,i))
          endif
        endif
      enddo
c .....................................................................
c
c ... 
      if( pedge(idCell) .eq. -2) sp = sp + const
c .....................................................................
c
c ...
      a(idcell) = sp
      do i = 1, nshared
         a(idCell)= a(idCell) + a(i)
      enddo
c ......................................................................
c
c ... derivada temporal
      m    = k(2,idCell)*cc*area(idCell)/dt
c ... Backward de segunda ordem 
      if(bs) m = 1.5d0*m
      a(idCell) = a(idCell) + m
c .....................................................................
c
c ... under-relaxation(simple)
      a(idCell) = a(idCell)/underU 
      p         = p + (1-underU)*a(idCell)*u(idCell)
c .....................................................................
c
c ...
      rm  = 0.0d0
      aNb = 0.0d0
      do i = 1, nshared
        rm  = rm  + a(i)*u(i)
        aNb = aNb + a(i)
      enddo
c .....................................................................
c
c ... interpolacao
      iM(iws2  ,idCell)  = rm + p 
      iM(iws2+2,idCell)  = a(idCell)
c ... campo D para velocidade de correcao
      if(simpleC) then               
        d(iws2,idCell)     = area(idCell) / (a(idCell) -aNb) 
      else
        d(iws2,idCell)     = area(idCell) /  a(idCell) 
      endif
c .....................................................................
c
c ... gradP
      p   = p - pf 
c      p   = p - gradP(iws2,idCell)*area(idCell)
c .....................................................................
c
c ... residou da celula
      rm  = p + rm - a(idCell)*u(idCell)
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
      grad1(1,1) = 0.0d0
      grad1(2,1) = 0.0d0
      do i = 1, nshared
        uf = um(i)
        grad1(1,1) = grad1(1,1) + uf*n(1,i)
        grad1(2,1) = grad1(2,1) + uf*n(2,i)
      enddo
      grad1(1,1) = grad1(1,1) / area(1)
      grad1(2,1) = grad1(2,1) / area(1)
c ... funcao de limatacao do fluxo para o termo advectivo 
      fluxl(1) = 1.0d0
      return
c .....................................................................
c
c .. F = dt*u0*V
  300 continue
      m    = k(2,idCell)*cc*area(idCell)/dt
      p    = m*u(idCell)
c ... Backward de segunda ordem
      m    = k(2,idCell)*cc*area(idcell)/dt 
      if(bs) p = 2*p - 0.5d0*m*u0(idcell)
      return  
c .....................................................................
c
c .....................................................................
c
c ...
  400 continue
      return
  500 continue
      return
  600 continue
      return
c .....................................................................
c .. CFL e Reynalds
  700 continue
c ... cfl numbe
      cfl = (dt*vm)/dsqrt(area(idcell))
      Mp(1) = cfl 
c .....................................................................
c
c ... Reynolds number
      re = (specificMassA*vm*dsqrt(area(idcell)))/viscosity
      Mp(2) = re
c .....................................................................
c
c ...
      Mp(3) = area(idCell)
c .....................................................................
c
c ...
      Mp(4) = vm          
c .....................................................................
c
c ... 
      Mp(5) = 0.0d0              
c .....................................................................
c
c ... 
      Mp(6) = 0.0d0
c ......................................................................
c
c ... massa 
      Mp(7) = k(2,idCell)*area(idCell)
c ......................................................................
c 
c ...
      p   = 0.0d0
      dfd = 0.0d0
      do i = 1, nshared
        viznel = viz(i)
c ... contorno      
        if( viznel .lt. 0) then
c ... velocidade prescrita
          if(pedge(i) .eq. 2 ) then
c ... fluxo advectivo de primeira ordem
            wfn = sedge(1,i)*n(1,i) + sedge(2,i)*n(2,i)
            cv  = specificMassA*wfn*meta(i)
            dfd = dfd + cv
c ... pressao prescrita
          else if(pedge(i) .eq. 3) then
            wfn = w(1,idCell)*n(1,i) + w(2,idCell)*n(2,i)
            if(dabs(wfn) .gt. ZERO) then 
              if( wfn .gt. 0.d0 ) then 
                cv  = specificMassA*wfn*meta(i)
                p   = p + cv
              else if(wfn .lt. 0.d0) then
c               wfn = -dsqrt(2.d0*(sedge(3,i)-u0(idCell))/specificMassA)
                cv  =  specificMassA*wfn*meta(i)
                dfd =  dfd + cv
              endif
            endif
          endif 
        endif 
      enddo
c ......................................................................
c
c ...
      Mp(8) = dfd     ! entrada de massa 
      Mp(9) = p       ! saida de massa
c ......................................................................
c
c ... Energia cinetica 
      Mp(10) = 0.5d0*specificMassA*vm*vm*area(idCell)
c ......................................................................
c
c ... dissipacao de energia da escala resolvivel
c ... S:S - contracao do tensor de tensoes
      s11     =   grad1(1,idCell)                 ! du1/dx1                 
      s22     =   grad2(2,idCell)                 ! du2/dx2
      s12     =   grad1(2,idCell)+grad2(1,idCell) ! du1/dx2 + du2/dx1
      modS   =   2.0d0*(s11*s11+s22*s22)+s12*s12       ! S:S
      Mp(11) =   (viscosity/specificMassA)*modS
c ......................................................................
      return  
      end
c **********************************************************************
c
c **********************************************************************
c * CELL_SI_WM_LES_EB_I_GGL : Celula 2D com termo de correção          *
c * para o fluxo difusivo e convectivo para equação de momentun        *
c * formulacao incomepressivel com lei de parade e LES                 *
c * (EDGE-BASE-TVD)                                                    *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * a        - nao definido                                            *
c * x        - coordenadas dos vertices da celula central e do seus    *
c * vizinhos                                                           *
c * u        - variavel da iteracao anterior                           *
c * u0       - variavel da iteracao anterior                           *
c * u1       - variavel da iteracao anterior                           *
c * rO       - massa especifica da celula e sua vizinhas em            *
c * t -1, t e t+1                                                      *
c * w        - campo de velocidade estimado                            *
c * d        - nao definido                                            *
c * grad1    - gradiente u(iws2) reconstruido da celula                *
c * grad2    - gradiente u(2 ou 1) reconstruido da celula              *
c * gradP    - gradiente p  reconstruido da celula                     *
c * fluxl    - limitador de fluxo na celulas                           *
c * k        - propriedades da celula e do seus vizinhos               *
c * rm       - nao definido                                            *
c * im       - nao definido                                            *
c *  p       - nao definido                                            *
c * sedge    - valores da condicoes de contorno por aresta             *
c * dt       - passo de tempo                                          *  
c * pedge    - tipo de condicao de contorno por aresta                 *
c * viz      - vizinhos da celula                                      *
c * nen      - numeros de nos por celula                               *
c * nshared  - numeros de faces por celulas                            *
c * iws1     - 1 - sistema de equacoes                                 *
c *          - 2 - reconstruacao de gradiente                          *
c * iws2     - 1 - direcao 1                                           *
c *          - 2 - direcao 2                                           *
c * nel      - numero da celula                                        *
c * sn(2,i)  - nos da aresta i                                         *
c * acod     - codigo para o calculo da area                           *
c *          ( 3 - triangulo; 4 - quadrilatero)                        *
c * mP       - parametros da equacao de momentos(cfl,Reynalds)         *
c * eddyVisc - parametros da equacao de momentos(cfl,Reynalds)         *
c * bs       - Euller Backward de segunda ordem                        *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * a      - coeficientes                                              *
c * p      - vetor de forcas                                           *
c * d      - campo d                                                   *
c * rm(2)  - residou da equacao da quantidade de movmentos             *
c * im     - interploacao de momentos para velocidade nas faces        *
c * grad1   - gradiente (GREEN GAUSS LINEAR)                            *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine cell_si_wm_les_eb_i_ggl(a,x,u,u0,u1,w,d,grad1,grad2
     .                                  ,gradP,iM,fluxl,k,rm,p,sedge
     .                                  ,dt,pedge,viz,nshared,ndm
     .                                  ,iws1,iws2,nel,sn
     .                                  ,acod,mP,eddyVisc,bs)
      implicit none
      include 'simple.fi'
      include 'les.fi'
      real*8 areacell,uf,gf(2),gfn,gfk,w(ndm,*),p,wf(ndm),wfn,cv,cvc
      real*8 r,limit,du,rof,gpk,vm,m,cc,dt,df(3,4)
      real*8 fluxl(nshared+1),km(3,4),d(2,nshared+1),mP(*)
      real*8 viscosity,specificMassA,rm,iM(4,nshared+1) 
      real*8 a(*),sp,xc(3,5),x(ndm,nshared,nshared+1),dfd,nk,um(4),cd
      real*8 sedge(3,nshared+1),ksi(3,4),eta(3,4),n(3,4),meta(4)
      real*8 mksi(4),area(5),ap
      real*8 u(5),u0(5),u1(5)
      real*8 stressW(4),yPLus(4),yPLusMax
      real*8 ca(4),alpha,k(10,*),xm(3,4)
      real*8 eps,cfl,re,ug(4),xcg(3,4),kf
      real*8 gradP(ndm,*),grad2(ndm,*),grad1(ndm,*)
      real*8 pf,p1,p2,pface,mKsiF(4),aNb
c ... LES
      real*8 uc0,uc1,dMin
      real*8 eddyVisc(5),modS,modSd,filtro,tfiltro,lMin,viscEf1,viscEf2
      real*8 tFilterModSS(3),tFiltermodS,tFilterS(3),tFilterV(2)
      real*8 mLilly(3),leornadS(3)
      real*8 tFilterVv(3),volW,volWt,vel(2)
      real*8,parameter :: Cs  = 0.1d0
      real*8,parameter :: Cw  = 0.17d0
      real*8,parameter :: r23 = 0.666666666666667d0 ! 2.0/3.0
      real*8,parameter :: r13 = 0.333333333333333d0 ! 1.0/3.0
      real*8,parameter :: VanDriest = 26.d0
      real*8,parameter :: vonKarman = 0.4187d0
      real*8 s11,s22,s12,g11,g12,g21,g22
      logical wall
c .......................................................................
      integer pedge(*),viz(*),viznel,iws1,iws2,nel
      integer ndm,nshared,i,j,icod,acod,idcell,sn(2,*),l
      real*8  const,ZERO
      logical bs
      parameter (ZERO =  1.0d-4)
      parameter (eps   = 1.0d-14)
      parameter (const = 1.0d60)
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
      wall   = .false.
c .....................................................................
c
c ...      
       call cellGeom2D(x,mksi,meta,ksi,eta,xc,area,n,xm,df,ca
     .                ,viz,nshared,ndm,sn,acod)
c .....................................................................
c
c ... propriedeade da celula
      viscosity    =  k(1,idCell)
      specificMassA=  k(2,idCell)
      cc           =  k(6,idCell)
      icod         =  k(7,idCell)
      vm           = dsqrt(w(1,idCell)*w(1,idCell) 
     .                   + w(2,idCell)*w(2,idCell))
c .....................................................................
c
c ... calculo de yPLus e Tensao da parede
      yPLus(1:nshared)   = 0.0d0
      stressW(1:nshared) = 0.0d0
      yPLusMax           = 0.0d0
      dMin               = 0.0d0
      if( iws1  .eq. 1 .or. iws1 .eq. 8 ) then 
        do i = 1, nshared
          viznel = viz(i)
c ...
          if( viznel .lt. 0) then
            if(pedge(i) .eq. 0 ) then
              wall = .true.
c ... velocidade tangencial a face
              wfn = dabs(w(1,idCell)*eta(1,i) + w(2,idCell)*eta(2,i))
              if(wfn .ne. 0.0d0) then
                call wallModel(yPlus(i)     ,stressW(i),viscosity
     .                        ,specificMassA,wfn,ca(i),nel,wallModelF)
              endif
            endif
          endif
c .....................................................................
        enddo
c ...
        yPlusMax = yPlus(1)
        dMin     =    ca(1) 
        do i = 2, nshared
          yPlusMax = max(yPlusMax,yPlus(i))
          dMin     = min(dMin    ,   ca(i)) 
        enddo
      endif
c .....................................................................
c
c ...
      if( iws1 .eq. 7) goto 700
      if( iws1 .eq. 8) goto 800
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
c ... considerando o gradiente normal para estimar o potencial no centro
c     da celula fantasmas)
          if(pedge(i) .eq. 0 ) then
           ug(i) = -u(idCell)
c ... variavel prescrita prescrita (extrapolacao de primeira ordem para a centro o
c     da celula fantasmas)
          else if(pedge(i) .eq. 1 .or. pedge(i) .eq. 2) then
            ug(i) = sedge(iws2,i)
c ... pressao prescrita            
          else if(pedge(i) .eq. 3 .or. pedge(i) .eq. 4) then
            ug(i) = u(idCell)
c ... tensao na parede (         
          else if(pedge(i) .eq. 5) then
            ug(i) = u(idCell) - 2.0d0*ca(i)*sedge(iws2,i)/viscosity 
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
      goto(100,200,300,400,500,600) iws1
c .....................................................................
c
  100 continue
c ...
      p            = 0.d0
      pf           = 0.d0
      sp           = 0.d0
      cd           = 0.d0
      cvc          = 0.d0
      loopCell: do i = 1, nshared
        viznel = viz(i)
c ... condicao de contorno      
        if( viznel .lt. 0) then
          a(i) = 0.d0
          ap   = viscosity*meta(i)
c ... termos viscosos
          if(iws2 .eq. 1) then
            p = p + ap*(grad1(iws2,idCell)*n(1,i) 
     .            + grad2(iws2,idCell)*n(2,i))
          else
            p = p + ap*(grad1(iws2,idCell)*n(2,i) 
     .            + grad2(iws2,idCell)*n(1,i))
          endif
c .....................................................................
c
c ... gradiente da pressao com resconstrucao de segunda ordem
c ... pressao prescrita
          if(pedge(i) .eq. 3) then  
            pFace   = sedge(3,i)
          else
            pFace   = u0(idCell)
c ... calculo da pressao na face atraves de reconstrucao
            if(sPressure) then
              pFace = pFace + (gradP(1,idCell)*df(1,i) 
     .                      +  gradP(2,idCell)*df(2,i))
            endif   
          endif
          pf      = pf + pFace*meta(i)*n(iws2,i)
c .....................................................................
c
c ... parade  
          if(pedge(i) .eq. 0 ) then
c ... sinal  
            wfn = w(1,idCell)*eta(1,i) + w(2,idCell)*eta(2,i)
            if(wfn .gt. 0.0d0) then
              stressW(i) = -stressW(i)*eta(iws2,i)*meta(i)
            else
              stressW(i) =  stressW(i)*eta(iws2,i)*meta(i)
            endif
c             print*,yPlus,stressW,w(iws2,idCell),wfn,eta(iws2,i),nel
            p  = p + stressW(i)
c ... parade movel 
          else if(pedge(i) .eq. 1) then
            ap = viscosity*meta(i)/ca(i)
            sp = sp + ap
            p  = p  + ap*sedge(iws2,i)
c ... velocidade prescrita
          else if(pedge(i) .eq. 2 ) then
c ... fluxo advectivo de primeira ordem
            wfn = sedge(1,i)*n(1,i) + sedge(2,i)*n(2,i)
            cv  = specificMassA*wfn*meta(i)
            p   =  p  - cv*sedge(iws2,i)
c ... pressao prescrita
          else if(pedge(i) .eq. 3) then
            wfn = w(1,idCell)*n(1,i) + w(2,idCell)*n(2,i)
            if( wfn .gt. 0.d0 ) then 
              cv  = specificMassA*wfn*meta(i)
              sp  = sp + cv
            else if(wfn .lt. 0.0) then
              cv  =  specificMassA*wfn*meta(i)
              p   =  p + cv*wfn
            endif
c ... localmente parabolica
          else if(pedge(i) .eq. 4) then
            wfn = w(1,idcell)*n(1,i) + w(2,idcell)*n(2,i)
            cv  = specificMassA*wfn*meta(i)
            sp  = sp + cv
c ... tensao tagencial
          else if(pedge(i) .eq. 5) then
            p   =  p + sedge(iws2,i)*meta(i)
c .....................................................................
          endif
c ...................................................................... 
c
c ...
        else
c ... interpolacao da propriedaes
          alpha = mKsiF(i)/mksi(i)
c ... viscosidade efetiva
          viscEf1= viscosity + eddyVisc(idCell)  
          viscEf2= k(1,i)    + eddyVisc(i) 
c ... media harmonica  visosidade molevular
          kf    = alpha/viscEf1 + (1.0d0-alpha)/viscEf2 
          kf    = 1.0d0/kf
c .....................................................................
          wf(1) = (1.0d0-alpha)*w(1,idcell)    + alpha*w(1,i)
          wf(2) = (1.0d0-alpha)*w(2,idcell)    + alpha*w(2,i)
          gf(1) = (1.0d0-alpha)*grad1(1,idcell) + alpha*grad1(1,i)
          gf(2) = (1.0d0-alpha)*grad1(2,idcell) + alpha*grad1(2,i)
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
          rof   = k(2,idCell)
          cv    = rof*wfn*meta(i)
c ... fluxo convectivo de ordem superior
c
c ... central-differencing scheme
          if(icod .eq. 9) then
            uc0 = u(i)      
     .          + grad1(1,i)*(xm(1,i) - xc(1,i))
     .          + grad1(2,i)*(xm(2,i) - xc(2,i))  
            uc1 = u(idCell) 
     .          + grad1(1,idCell)*d(1,i) + grad1(2,idCell)*d(1,i)   
            du  = 0.5d0*(uc0 + uc1)              
            if(cv.lt.0.0d0) then
              cvc = du - u(idCell)
            else
              cvc = du - u(i)
            endif
c .....................................................................
c
c ... interpolcao TVD
          else
            if(cv.lt.0.0d0) then
              du    = u(idcell) - u(i) + eps 
              gpk   = grad1(1,i)*(xc(1,idcell)-xc(1,i))
     .              + grad1(2,i)*(xc(2,idcell)-xc(2,i))
              alpha = mKsiF(i)/mksi(i)
              r     = 2.0d0*gpk/du - 1.0d0
              cvc   = (1.0d0-alpha)*limit(r,icod)*du
            else  
              du    = u(i) - u(idcell) + eps 
              gpk   = (grad1(1,idCell)*ksi(1,i)
     .              + grad1(2,idCell)*ksi(2,i))*mksi(i)
              alpha = mKsiF(i)/mksi(i)   
              r     = 2.0d0*gpk/du - 1.0d0
              cvc   = alpha*limit(r,icod)*du
            endif
c ... interpolacao unidirecional  
            cvc   = cvc +  gf(1)*km(1,i) + gf(2)*km(2,i)
c .....................................................................
          endif
c .....................................................................
c
c ...
          a(i) = a(i) + max(-cv,0.0d0)
          sp   = sp   + cv
          p    =  p   - cv*cvc
c .....................................................................
c
c ... gradiente da pressao com resconstrucao de segunda ordem
          alpha = mKsiF(i)/mksi(i)
          if(sPressure) then
            p1   = u0(idCell) + gradP(1,idCell)*df(1,i) 
     .                        + gradP(2,idCell)*df(2,i)
            p2   = u0(i)      + gradP(1,i)*(xm(1,i)-xc(1,i))
     .                        + gradP(2,i)*(xm(2,i)-xc(2,i))
            pFace = 0.5d0*(p1 + p2)
          else
            p1    = u0(idCell) 
            p2    = u0(i)  
            pFace = (1-alpha)*p1 + alpha*p2
          endif
          pf   = pf + pFAce*meta(i)*n(iws2,i)
c .....................................................................
c
c ... termos nao lineares da viscosidade
          alpha = mKsiF(i)/mksi(i)
          gf(1) = (1.0d0-alpha)*grad1(iws2,idcell) + alpha*grad1(iws2,i)
          gf(2) = (1.0d0-alpha)*grad2(iws2,idcell) 
     .          +        alpha *grad2(iws2,i)
          if( iws2 .eq. 1) then
            p     = p + kf*meta(i)*(gf(1)*n(1,i)+ gf(2)*n(2,i))
          else if(iws2 .eq. 2 ) then
            p     = p + kf*meta(i)*(gf(1)*n(2,i) + gf(2)*n(1,i))
          endif
        endif
      enddo loopCell
c .....................................................................
c
c ... 
      if( pedge(idCell) .eq. -2) sp = sp + const
c .....................................................................
c
c ...
      a(idcell) = sp
      do i = 1, nshared
         a(idCell)= a(idCell) + a(i)
      enddo
c ......................................................................
c
c ... derivada temporal
      m    = k(2,idCell)*cc*area(idCell)/dt
c ... Backward de segunda ordem 
      if(bs) m = 1.5d0*m
      a(idCell) = a(idCell) + m
c .....................................................................
c
c ... under-relaxation(simple)
      a(idCell) = a(idCell)/underU 
      p         = p + (1-underU)*a(idCell)*u(idCell)
c .....................................................................
c
c ...
      rm  = 0.0d0
      aNb = 0.0d0
      do i = 1, nshared
        rm  = rm  + a(i)*u(i)
        aNb = aNb + a(i)
      enddo
c .....................................................................
c
c ... interpolacao
      iM(iws2  ,idCell)  = rm + p 
      iM(iws2+2,idCell)  = a(idCell)
c ... campo D para velocidade de correcao
      if(simpleC) then               
        d(iws2,idCell)     = area(idCell) / (a(idCell) -aNb) 
      else
        d(iws2,idCell)     = area(idCell) /  a(idCell) 
      endif
c .....................................................................
c
c ... gradP
      p   = p - pf 
c      p   = p - gradP(iws2,idCell)*area(idCell)
c .....................................................................
c
c ... residou da celula
      rm  = p + rm - a(idCell)*u(idCell)
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
      grad1(1,1) = 0.0d0
      grad1(2,1) = 0.0d0
      do i = 1, nshared
        uf = um(i)
        grad1(1,1) = grad1(1,1) + uf*n(1,i)
        grad1(2,1) = grad1(2,1) + uf*n(2,i)
      enddo
      grad1(1,1) = grad1(1,1) / area(1)
      grad1(2,1) = grad1(2,1) / area(1)
c ... funcao de limatacao do fluxo para o termo advectivo 
      fluxl(1) = 1.0d0
      return
c .....................................................................
c
c .. F = dt*u0*V
  300 continue
      m    = k(2,idCell)*cc*area(idCell)/dt
      p    = m*u(idCell)
c ... Backward de segunda ordem
      m    = k(2,idCell)*cc*area(idcell)/dt 
      if(bs) p = 2*p - 0.5d0*m*u0(idcell)
      return  
c .....................................................................
c
c .....................................................................
c
c ...
  400 continue
      return
  500 continue
      return
  600 continue
      return
c .....................................................................
c
c .. CFL e Reynalds
  700 continue
c ... cfl number
      cfl = (dt*vm)/dsqrt(area(idcell))
      Mp(1) = cfl 
c .....................................................................
c
c ... Reynolds number
      re = specificMassA*vm*dsqrt(area(idcell))
      re = re/(viscosity+eddyVisc(idCell))
      Mp(2) = re               
c .....................................................................
c
c ...
      Mp(3) = area(idCell)
c .....................................................................
c
c ...
      Mp(4) = vm          
c .....................................................................
c
c ... 
      Mp(5) = 0.0d0              
c .....................................................................
c
c ... 
      Mp(6) = 0.0d0
c ......................................................................
c
c ... massa 
      Mp(7) = k(2,idCell)*area(idCell)
c ......................................................................
c 
c ...
      p   = 0.0d0
      dfd = 0.0d0
      do i = 1, nshared
        viznel = viz(i)
c ... contorno      
        if( viznel .lt. 0) then
c ... velocidade prescrita
          if(pedge(i) .eq. 2 ) then
c ... fluxo advectivo de primeira ordem
            wfn = sedge(1,i)*n(1,i) + sedge(2,i)*n(2,i)
            cv  = specificMassA*wfn*meta(i)
            dfd = dfd + cv
c ... pressao prescrita
          else if(pedge(i) .eq. 3) then
            wfn = w(1,idCell)*n(1,i) + w(2,idCell)*n(2,i)
            if( wfn .gt. 0.d0 ) then 
              cv  = specificMassA*wfn*meta(i)
              p   = p + cv
            else if(wfn .lt. 0.d0) then
c               wfn = -dsqrt(2.d0*(sedge(3,i)-u0(idCell))/specificMassA)
              cv  =  specificMassA*wfn*meta(i)
              dfd =  dfd + cv
            endif
c ...................................................................... 
          endif 
        endif 
      enddo
c ......................................................................
c
c ...
      Mp(8) = dfd     ! entrada de massa 
      Mp(9) = p       ! saida de massa
c ......................................................................
c
c ... Energia cinetica 
      Mp(10) = 0.5d0*specificMassA*vm*vm*area(idCell)
c ......................................................................
c
c ... dissipacao de energia da escala resolvivel
c ... S:S - contracao do tensor de tensoes
      s11     =   grad1(1,idCell)                 ! du1/dx1                 
      s22     =   grad2(2,idCell)                 ! du2/dx2
      s12     =   grad1(2,idCell)+grad2(1,idCell) ! du1/dx2 + du2/dx1
      modS   =   2.0d0*(s11*s11+s22*s22)+s12*s12         ! S:S
      Mp(11) =   (viscosity/specificMassA)*modS
c ......................................................................
c
c ... transferencia de energa da escala resolvivel para escala nao
c     resolvivel
c ... S:S - contracao do tensor de tensoes
c     s11     =   grad1(1,idCell)                 ! du1/dx1                 
c     s22     =   grad2(2,idCell)                 ! du2/dx2
c     s12     =   grad1(2,idCell)+grad2(1,idCell) ! du1/dx2 + du2/dx1
c     modS   =   s11*s11+s22*s22+0.5d0*s12*s12         ! S:S
      Mp(12) =  (eddyVisc(idCell)/specificMassA)*modS
c ......................................................................
      return
c ......................................................................
c
c ...              
  800 continue
c ... vetor que une os centroides dos vizinho 
c     do i = 1, ndm
c       do l = 1, nshared
c ... centroide l  
c         eta(i,l) = xcg(i,sn(2,l)) - xcg(i,sn(1,l))    
c       enddo
c     enddo
c ... vetor normal com o modulo igual a arestas
c     do l = 1, nshared
c       n(1,l) = eta(2,l)
c       n(2,l) =-eta(1,l)
c     enddo
c .....................................................................
c
c ... area                 
c     area(1) = areacell(eta,acod)
c .....................................................................
c
c ... smagorinsky (grad1 -> grad2; grad2 -> gradU2)
      if(lesModel .eq. 1 .or. lesModel .eq. 0) then 
c ... |S|
        s11     =   grad1(1,idCell)                 ! du1/dx1                 
        s22     =   grad2(2,idCell)                 ! du2/dx2
        s12     =   grad1(2,idCell)+grad2(1,idCell) ! du2/dx1 + du1/dx2
        modS   = dsqrt(2.0d0*(s11*s11+s22*s22)+s12*s12) 
c ... Cs*filtro
        filtro = dsqrt(area(idCell))
        lMin   = Cs*filtro
c ... near wall Van Driest
        if(wall) then
c        lMin= min(lMin,lMin*(1.0d0-dexp(-yPlusMax/vanDriest)))
          lMin = min(vonKarman*dMin,lMin)
          lMin = lMin*(1.0d0-dexp(-yPlusMax/vanDriest))
        endif
c ... viscosidade turbulenta
        eddyVisc(idCell) = specificMassA*(lMin)*(lMin)*modS
c .....................................................................
c
c ... Wall-Adpating Local Eddy-Viscosity (WALE) model
      elseif(lesModel .eq. 2) then
c ... S:S - contracao do tensor de tensoes
        s11     =   grad1(1,idCell)                 ! du1/dx1                 
        s22     =   grad2(2,idCell)                 ! du2/dx2
        s12     =   grad1(2,idCell)+grad2(1,idCell) ! du1/dx2 + du2/dx1
        modS   =   s11*s11+s22*s22+0.5d0*s12*s12         ! S:S
c .....................................................................
c
c ... Sd:Sd
        g11    = grad1(1,idCell)*grad1(1,idCell) ! du1/dx1*du1/dx1
        g22    = grad2(2,idCell)*grad2(2,idCell) ! du2/dx2*du2/dx2
        g12    = grad1(2,idCell)*grad1(2,idCell) ! du1/dx2*du1/dx2
        g21    = grad2(1,idCell)*grad2(1,idCell) ! du2/dx1*du2/dx1
        s12     = r13*(g11+g22)
        s11     = g11 - s12
        s22     = g22 - s12
        s12     = g12 + g21
        modSd  = s11*s11+s22*s22+0.5d0*s12*s12           !Sd:Sd
c .....................................................................
c
c ...
        modS   = (modSd**(1.5d0))/(modS**(2.5d0) + modSd**(1.25)+eps)
c .....................................................................
c
c ... Cs*filtro
        filtro = dsqrt(area(idCell))
        lMin   = Cw*filtro
c ... near wall Van Driest
c       yPlusMax = 0.0d0
        if(wall) then
c        lMin= min(lMin,lMin*(1.0d0-dexp(-yPlusMax/vanDriest)))
          lMin = min(vonKarman*dMin,lMin)
c         lMin = lMin*(1.0d0-dexp(-yPlusMax/vanDriest))
        endif
c ... viscosidade turbulenta
        eddyVisc(idCell) = specificMassA*(lMin)*(lMin)*modS
c        print*,nel,modS,modSd,g11,g22,g12,g21
c .....................................................................
c 
c ... Vreman's Model                                   
      elseif(lesModel .eq. 3) then
c ... filtro*filtro
        filtro = area(idCell)
c ... alpha
        g11    = grad1(1,idCell)                  ! du1/dx1
        g22    = grad2(2,idCell)                  ! du2/dx2
        g21    = grad1(2,idCell)                  ! du1/dx2
        g12    = grad2(1,idCell)                  ! du2/dx1
c ...
        s11     = filtro*(g11*g11+g21*g21)
        s22     = filtro*(g12*g12+g22*g22)
        s12     = filtro*(g11*g12+g21*g22)
        modS   = s11*s22 - s12*s12                             !B           
        modSd  = g11*g11 + g22*g22 + g12*g12 + g21*g21     !alpha:alpha 
c .....................................................................
c
c ...
        modS   = dsqrt((modS+eps)/(modSd+eps))
c .....................................................................
c
c ... viscosidade turbulenta
        eddyVisc(idCell) = specificMassA*2.5d0*Cs*Cs*modS
c       print*,nel,modS,modSd,g11,g22,g12,g21 
c       print*,nel, eddyVisc(idCell)
c ......................................................................
c
c ... germando - lilly                                
      elseif(lesModel .eq. 4) then
        volW   = area(idCell)
        volWt  = volW           
c ... termo |S|Sij para o filtro de teste
        s11    =   grad1(1,idCell)                    ! du1/dx1                 
        s22    =   grad2(2,idCell)                    ! du2/dx2
        s12    = 0.5d0*(grad1(2,idCell)+grad2(1,idCell))! du1/dx2 + du2/dx1
        modS   = dsqrt(s11*s11+s22*s22+2.0d0*s12*s12)
        modSd  = modS
c ... vol*modS*S11
        tFilterModSS(1) = modS*s11*volW 
c ... vol*modS*S22       
        tFilterModSS(2) = modS*s22*volW
c ... vol*modS*S12 = vol*modS*S21
        tFilterModSS(3) = modS*s12*volW  
c ... termo |S| para o flitro de teste 
        tFilterModS     = modS*volW             ! vol*modS 
c ... termo Sij para o flitro de teste 
c ... vol*S11 
        tFilterS(1)   = s11*volW               
c ... vol*S22 
        tFilterS(2)   = s22*volW               
c ... vol*S12 = vol*S21 
        tFilterS(3)   = s12*volW               
c ... termo vv para para o filtro de teste
        vel(1)        = w(1,idCell)
        vel(2)        = w(2,idCell)
c ... vol*v1*v1
        tFilterVv(1)  = vel(1)*vel(1)*volW      
c ... vol*v2*v2
        tFilterVv(2)  = vel(2)*vel(2)*volW      
c ... vol*v1*v2 = vol*v2*v1
        tFilterVv(3)  = vel(1)*vel(2)*volW      
c ... termo v para para o filtro de teste
c ... vol*v1
        tFilterV(1)   = vel(1)*volW                   
c ... vol*v2
        tFilterV(2)   = vel(2)*volW      
        
        do i = 1, nshared
          viznel = viz(i)
          if( viznel .gt. 0) then
            volW   = area(i)
            volWt  = volWt + volW           
c ... termo |S|Sij para o filtro de teste
            s11    =   grad1(1,i)                   ! du1/dx1                 
            s22    =   grad2(2,i)                   ! du2/dx2
            s12    = 0.5d0*(grad1(2,i)+grad2(1,i))  ! du1/dx2 + du2/dx1
            modS   = dsqrt(s11*s11+s22*s22+2.0d0*s12*s12)
c ... vol*modS*S11
            tFilterModSS(1) = tFilterModSS(1) + modS*s11*volW 
c ... vol*modS*S22       
            tFilterModSS(2) = tFilterModSS(2) + modS*s22*volW
c ... vol*modS*S12 = vol*modS*S21
            tFilterModSS(3) = tFilterModSS(3) + modS*s12*volW  
c ... termo |S| para o flitro de teste 
            tFilterModS     = tFilterModS + modS*volW    ! vol*modS 
c ... termo Sij para o flitro de teste 
c ... vol*S11 
            tFilterS(1)   = tFilterS(1) + s11*volW               
c ... vol*S22 
            tFilterS(2)   = tFilterS(2) + s22*volW               
c ... vol*S12 = vol*S21 
            tFilterS(3)   = tFilterS(3) + s12*volW               
c ... termo vv para para o filtro de teste
            vel(1)        = w(1,i)
            vel(2)        = w(2,i)
c ... vol*v1*v1
            tFilterVv(1)  = tFilterVv(1) + vel(1)*vel(1)*volW      
c ... vol*v2*v2
            tFilterVv(2)  = tFilterVv(2) + vel(2)*vel(2)*volW      
c ... vol*v1*v2 = vol*v2*v1
            tFilterVv(3)  = tFilterVv(3) + vel(2)*vel(1)*volW      
c ... termo v para para o filtro de teste
c ... vol*v1
            tFilterV(1)   = tFilterV(1) + vel(1)*volW                  
c ... vol*v2
            tFilterV(2)   = tFilterV(2) + vel(2)*volW
          endif 
        enddo
c ......................................................................
c
c ...
        tFilterModSS(1:3) =  tFilterModSS(1:3)/volWt 
        tFilterModS       =  tFilterModS/volWt 
        tFilterS(1:3)     =  tFilterS(1:3)/volWt 
        tFilterVv(1:3)    =  tFilterVv(1:3)/volWt 
        tFilterV(1:2)     =  tFilterV(1:2)/volWt
c ......................................................................
c
c ... 
c       tfiltro = dsqrt(volWt) 
c       filtro  = dsqrt(aread(idCell)) 
        filtro  = area(idCell)             ! flitro *flitro
        tfiltro = volWt                    ! tflitro*tflitro
c ......................................................................
c
c ...  
c .. L11
        leornadS(1)       = tFilterVv(1) - tFilterV(1)*tFilterV(1)  
c .. L22
        leornadS(2)       = tFilterVv(2) - tFilterV(2)*tFilterV(2) 
c .. L12 = L21
        leornadS(3)       = tFilterVv(3) - tFilterV(1)*tFilterV(2)
c ...
        g11   = tfiltro*tFilterModS
c .. M11
        mLilly(1) = filtro*tFilterModSS(1) - g11*tFilterS(1) + 1.0d-60
c .. M22
        mLilly(2) = filtro*tFilterModSS(2) - g11*tFilterS(2) + 1.0d-60
c .. M12 = M21
        mLilly(3) = filtro*tFilterModSS(3) - g11*tFilterS(3) + 1.0d-60
c ......................................................................
c
c ... LijMij 
        g11 = leornadS(1)*mLilly(1) + leornadS(2)*mLilly(2)
     .      + 2.0d0*leornadS(3)*mLilly(3)
c ... MijMij 
        g22 = mLilly(1)*mLilly(1) + mLilly(2)*mLilly(2)
     .      + 2.0d0*mLilly(3)*mLilly(3)
c ......................................................................
c 
c ... cs² = LijMij/MijMij
        g11 = 0.5d0*(g11/(g22))
c ... limita o calor de Cs entre 0.0 - Cs
        g11 = min(g11,0.04d0)
        g11 = max(g11,0.0d0)
c ... media dos Cs²
        u1(idCell) = g11
        g11        = g11*area(idCell)
        do i = 1, nshared
          viznel = viz(i)
          if( viznel .gt. 0) then
            g11        = g11 + u1(i)*area(i)
          endif
        enddo
        g11 = g11/volWt
c ...
        eddyVisc(idCell) = specificMassA*g11*filtro*modSd
c ......................................................................
      endif
c ......................................................................
c
c ...                                      
      u0(1) = yPlusMax 
c .....................................................................
      return  
      end
c **********************************************************************
c
c **********************************************************************
c * CELLVEL :                                                          *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * a      - nao definido                                              *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine cellVel(wfn,w,grad,df,d,iM,ksi,n,p,mksi,meta,area,alpha
     .                  ,ndm,nshared,iws,viz)
      implicit none
      include 'simple.fi'
      real*8 vmf(2),w(ndm,*),df(*),ksi(3,*),n(3,*),mksi(*),iM(4,*)
      real*8 grad(ndm,*),p(*),alpha,gf(2),gfp(2),area(*)
      real*8 nk,ks(3),wfn,d(2,nshared+1),meta(*)
      integer iws,ndm,nShared,idCell,viz
      idCell = nShared + 1 
      goto (100,200,300,400,500,600) iws 
c ... interpola linear
  100 continue
        vmf(1) = (1-alpha)*w(1,idCell) + alpha*w(1,viz)
        vmf(2) = (1-alpha)*w(2,idCell) + alpha*w(2,viz)
         wfn   = vmf(1)*n(1,viz) + vmf(2)*n(2,viz)
         wfn   = wfn*meta(viz)
      return
c .....................................................................
c
c ... interpolacao Rhie-Chow
  200 continue
        vmf(1) = (1-alpha)*w(1,idCell) + alpha*w(1,viz)
        vmf(2) = (1-alpha)*w(2,idCell) + alpha*w(2,viz)
        wfn    = vmf(1)*n(1,viz) + vmf(2)*n(2,viz)
c ... interpolacao linear dos gradientes das pressoes
        gf(1)  = (1-alpha)*grad(1,idcell) + alpha*grad(1,viz)
        gf(2)  = (1-alpha)*grad(2,idcell) + alpha*grad(2,viz)
        gfp(1) = ((p(viz) - p(idcell))/mksi(viz))
        gfp(2) = gf(1)*ksi(1,viz) + gf(2)*ksi(2,viz)
c ... 
        wfn    = wfn +0.5d0*(df(1)+df(2))*(gfp(2)-gfp(1))
        wfn    = wfn*meta(viz)
      return
c .....................................................................
c
c ... intepolacao pela equaco de momentum    
  300  continue
c          gfp(1) = p(viz) - p(idCell)
c ... vetor de area         
c          vA(1) = n(1,viz)*meta(viz)
c          vA(2) = n(2,viz)*meta(viz)
c ... distancia entre os dois vetores
c          ks(1) = ksi(1,viz)*mksi(viz)
c          ks(2) = ksi(2,viz)*mksi(viz)
c ...
c          ndd   = gfp(1)*(df(1)*vA(1)*vA(1)+df(2)*vA(2)*vA(2))
c          gfk   = ndd/(ks(1)*vA(1)+ks(2)*vA(2))
c ... interpolacao linear dos gradientes das pressoes
c          gf(1) = (1-alpha)*grad(1,idCell) + alpha*grad(1,viz)
c          gf(2) = (1-alpha)*grad(2,idCell) + alpha*grad(2,viz)
c          nk    = 1/(ksi(1,viz)*n(1,viz) + ksi(2,viz)*n(2,viz))
c ... k = S - E = Sn - (S/ne)e
c          k(1)  = (n(1,viz) - nk*ksi(1,viz))*meta(viz)
c          k(2)  = (n(2,viz) - nk*ksi(2,viz))*meta(viz)
c ...         
c          gfn   = df(1)*gf(1)*k(1) + df(2)*gf(2)*k(2)
c ... momemtum interpolacao
c          aU1   = (1-alpha)*(iM(1,idCell)/iM(3,idCell)) 
c     .         +    alpha *(iM(1,viz   )/iM(3,viz   ))
c          aU2   = (1-alpha)*(iM(2,idCell)/iM(4,idCell)) 
c     .         +    alpha *(iM(2,viz   )/iM(4,viz   ))
c ...
c        aU1   = aU1*vA(1)
c        aU2   = aU2*vA(2)
c ...
c        wfn   = aU1+aU2 - (gfk+gfn)
c ... interpolacao das pseudo-velocidades
c      else if(iws.eq.4) then 
c        aU1    = (1-alpha)*(iM(1,idCell)/iM(3,idCell)) 
c     .         +    alpha *(iM(1,viz   )/iM(3,viz   ))
c        aU2    = (1-alpha)*(iM(2,idCell)/iM(4,idCell)) 
c     .         +    alpha *(iM(2,viz   )/iM(4,viz   ))
c        vmf(1) = aU1 
c        vmf(2) = aU2
c         wfn   = vmf(1)*n(1,viz) + vmf(2)*n(2,viz)
c         wfn   = wfn*meta(viz)
      return
c .....................................................................
c
c ...
  400 continue
      return
c .....................................................................
c
c ... Ferziger-Precic
  500 continue
        vmf(1) = (1-alpha)*w(1,idCell) + alpha*w(1,viz)
        vmf(2) = (1-alpha)*w(2,idCell) + alpha*w(2,viz)
         wfn   = vmf(1)*n(1,viz) + vmf(2)*n(2,viz)        
c ... interpolacao linear dos gradientes das pressoes
        gf(1)  = (1-alpha)*grad(1,idcell) + alpha*grad(1,viz)
        gf(2)  = (1-alpha)*grad(2,idcell) + alpha*grad(2,viz)
        gfp(1) = p(viz) - p(idCell)
        ks(1)  = mksi(viz)*ksi(1,viz) 
        ks(2)  = mksi(viz)*ksi(2,viz)
        nk     = ks(1)*n(1,viz) + ks(2)*n(2,viz)
        gfp(2) = gf(1)*ks(1)  + gf(2)*ks(2)
        gfp(2) = (gfp(2)-gfp(1))/nk
c ... 
        wfn = wfn + 0.5d0*(df(1)+df(2))*gfp(2)
        wfn = wfn*meta(viz)
      return
c .....................................................................
c
c ... Verteeg - Malalasekera 
  600 continue 
c        vmf(1) = (1-alpha)*w(1,idCell) + alpha*w(1,viz)
c        vmf(2) = (1-alpha)*w(2,idCell) + alpha*w(2,viz)
c         wfn   = vmf(1)*n(1,viz) + vmf(2)*n(2,viz)        
c ... 
c        gf(1)  = (0.5d0*(df(1)+df(2)))*((p(idCell)-p(viz))/mksi(viz))
c ...        
c        aU1    = 0.5d0*(d(1,idCell)+d(2,idCell))
c        aU2    = 0.5d0*(   d(1,viz)+   d(2,viz))
c        gfp(1) = grad(1,idCell)*ksi(1,idCell) 
c     .         + grad(2,idCell)*ksi(2,idCell)
c        gfp(2) = grad(1,viz)*ksi(1,viz) 
c     .         + grad(2,viz)*ksi(2,viz)
c        gf(2)  = (1-alpha)*aU1*gfp(1) + alpha*aU2*gfp(2)
c ...        
c        wfn = wfn + gf(1) - gf(2)
      return
c .....................................................................
c
c ...
      return
      end
c **********************************************************************                           
