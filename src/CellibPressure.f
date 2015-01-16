c **********************************************************************
c * CELL_P  : Celula 2D de para equacao de correcao de pressao do      *
c * metodo simple para gas ideal incompressivel                        *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * a      - nao definido                                              *
c * x      - coordenadas dos vertices da celula central e do seus      *
c * vizinhos                                                           *
c * u      - temperatura por celula                                    *
c * u0     - campo de pressao estimado                                 *
c * u1     - temperatura da celula do passo anterior                   *
c * w      - campo de velociade conhecido                              *
c * k      - propriedades da celula e do seus vizinhos                 *
c * grad   - gradiente reconstruido na celula                          *
c * gradP  - gradiente de pressao reconstruido na celula               *
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
c * iws1   - 1 - pressao de correcao (reconstrucao de gradiente)       *
c *        - 2 - pressao             (reconstrucao de gradiente)       *
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
      subroutine cell_si_p_idg_ggl(a,x,u,u0,u1,ro,w,d,grad,gradP,iM
     .                            ,fluxl,k,rm,p,sedge,dt,pedge,viz
     .                            ,nshared,ndm,iws,iws1
     .                            ,nel,sn,acod,bs)
      implicit none
      include 'simple.fi'
c ... variavel externas
      real*8 a(*),rm,p,sedge(3,*),ap,grad(ndm,*),fluxl(nshared+1)
      real*8 dt,iM(4,nshared+1),ro(3,*),w(ndm,*),gradP(ndm,*) 
      integer pedge(*),viz(*),iws,iws1
      integer nel,acod,nshared,ndm,sn(2,*)
c ... variavel internas
      real*8 xc(3,5),x(ndm,nshared,nshared+1),d(2,nshared+1)
      real*8 gfk,umin,umax,r,specificMass,m,cc,du,dfd,gpn,kf,nk,ndd
      real*8 um(4),kis(3,4),ksi(3,4),eta(3,4),n(3,4),meta(4),mksi(4)
      real*8 u(5),u0(5),u1(5),ca(4),k(10,*),xm(3,4),gf(2),gfp(2),area(5)
      real*8 areacell,ug(4),xcg(3,4),alpha,dd(2),vm(2,4),vmf(2),uf
      real*8 limitv,eps,df(3,4),par(2),aux,pl,ddf(2),pre(nshared+1)
      real*8 wfn,sp,cvc,cv,cp,cf,specificMassA,mKsiF(4)
      real*8 const,km(3,4),ZERO,dVirtualP(2),dVirtualViz(2),gama
      integer i,j,l,idcell,viznel,icod,tRo1,tRo2,cor
      parameter (ZERO  = 1.0d-14)
      logical bs
      parameter (const = 1.0d60)
      external limitv,areacell
c ...      
      cor = intVel
c ...
      icod   = 1
      idCell = nshared + 1
      tRo1   = 2
      tRo2   = 3
c ......................................................................
c
c ...
       call cellGeom2D(x,mksi,meta,ksi,eta,xc,area,n,xm,df,ca
     .                ,viz,nshared,ndm,sn,acod)
c .....................................................................
c
c ... propriedeade da celula
      specificMass = rO(tRo2,idCell)
      specificMassA=  k(2,idCell)
      cc           =  k(6,idCell)
      icod         =  k(7,idCell)
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
c ... condicao de controno para as celulas fantasma
          if(pedge(i) .eq. 3) then
            if(iws1 .eq. 2) then
              ug(i) = 0.0d0
            else    
              ug(i) = sedge(3,i)
            endif
c ...                              
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
c ....................................................................
  100 continue
      p  = 0.d0
      sp = 0.d0
      do i = 1, nshared
        viznel = viz(i)
c ... condicao de contorno      
        if( viznel .lt. 0) then
          a(i) = 0.d0
c ... velocidade prescrita
          if(pedge(i) .eq. 2 ) then
c ... fluxo advectivo de primeira ordem
            wfn = sedge(1,i)*n(1,i) + sedge(2,i)*n(2,i)
            cv  = specificMassA*wfn*meta(i)
            p   =  p  - cv
c ... pressao prescrita
          else if(pedge(i) .eq. 3) then
            ddf(1) = d(1,idcell)  
            ddf(2) = d(2,idcell) 
c ... produtos interno
            if( cor .eq. 2 .or. cor .eq. 3 ) then
              dd(1) = ddf(1)*n(1,i)
              dd(2) = ddf(2)*n(2,i)
              ndd   = dd(1)*n(1,i) + dd(2)*n(2,i)         
c ...                
              dfd  = (specificMass*ndd*meta(i))/ca(i)     
               sp  = sp + dfd
            else if( cor .eq. 5) then  
              ndd   = 0.5d0*(ddf(1) + ddf(2))         
c ...              
              dfd   = (specificMass*meta(i)*ndd)/ca(i)
              sp  = sp + dfd
            endif   
c ... velocidade
            wfn = w(1,idcell)*n(1,i) + w(2,idcell)*n(2,i)
            if(dabs(wfn) .gt. ZERO) then
              if( wfn .gt. 0.d0  ) then 
                cv  = specificMass*wfn*meta(i)
              elseif( wfn .lt. 0.d0) then
c                wfn =-dsqrt(2.d0*(sedge(3,i)-u1(idCell))/specificMassA)
                cv  = specificMassA*wfn*meta(i)
              endif
            else
              cv = 0.0d0
            endif
            p     = p - cv    
c ... localmente parabolica
          else if(pedge(i) .eq. 4) then
c .....................................................................
          endif
        else
          if(cor .eq. 1 .or. cor .eq. 2 .or. cor .eq. 3) then      
c ... interpolacao da propriedades
            alpha = mKsiF(i)/mksi(i)
            kf = (1.0d0-alpha)*ro(tRo2,idCell) + alpha*ro(tRo2,i)
c ... aritimetica
            ddf(1) = (1.0d0-alpha)*d(1,idcell)  + alpha*d(1,i) 
            ddf(2) = (1.0d0-alpha)*d(2,idcell)  + alpha*d(2,i)
c ... produtos interno
            dd(1) = ddf(1)*n(1,i)
            dd(2) = ddf(2)*n(2,i)
            ndd   = dd(1)*n(1,i)    + dd(2)*n(2,i)         
            nk    = n(1,i)*ksi(1,i) + n(2,i)*ksi(2,i)
c ... difusao direta
            dfd   = (kf*meta(i)*ndd)/(nk*mksi(i))
c ...
            a(i)  = dfd
c ... interpolacao linear das velociadades
c ... Ferziger Peric            
          else if(cor .eq. 5) then
c ... interpolacao da propriedades
            alpha = mKsiF(i)/mksi(i)
            kf = (1.0d0-alpha)*ro(tRo2,idCell) + alpha*ro(tRo2,i)
c ... aritimetica
            ddf(1) = (1.0d0-alpha)*d(1,idcell)  + alpha*d(1,i) 
            ddf(2) = (1.0d0-alpha)*d(2,idcell)  + alpha*d(2,i)
c ... produtos interno
            ndd   = 0.5d0*(ddf(1) + ddf(2))         
            nk    = n(1,i)*ksi(1,i) + n(2,i)*ksi(2,i)            
c ... difusao direta
            dfd   = (kf*meta(i)*ndd)/(nk*mksi(i))
            a(i)  = dfd
          endif
c ... interpolacao linear das velociadades
          call cellVel(wfn,w,gradP,ddf,d,iM,ksi,n,u1,mksi,meta,area
     .                ,alpha,ndm,nshared,intVel,i)
c ...
          p   = p - kf*wfn
c ...           
c          cf = (1-alpha)*cp(u0(idCell),1) + alpha*cp(u0(i),1)
c          a(i) = a(i) -     alpha*cf*wfn*meta(i)
c          sp   = sp   + (1-alpha)*cf*wfn*meta(i)
        endif
      enddo
c .....................................................................
c
c ... fixa o valor da pressao de correcao igual a zero
      if( pedge(idCell) .eq. -1) sp = sp + const
c .....................................................................
c
c ...
c      p         = p*dt
c      sp        = sp*dt
      a(idcell) = sp
      do i = 1, nshared
c        a(i)      = a(i)*dt
        a(idcell) = a(idcell) + a(i)
      enddo
c .....................................................................
c
c ...
c      m = cp(u0(idCell),1)*area(idCell)
c      a(idcell) = a(idcell) + m
      if(bs) then
        cvc = (1.5d0*rO(3,idCell)-2.0d0*rO(2,idCell) 
     .      +  0.5d0*rO(1,idCell))
      else
        cvc = (rO(3,idCell)-rO(2,idCell))
      endif
      p = p - area(idCell)*cvc/dt
c .....................................................................
c
c ... residou da celula
c      rm  = 0.0d0
c       do i = 1, nshared
c         rm = rm + a(i)*u(i)
c       enddo
c ... 
c       rm  = p + rm - a(idCell)*u(idCell)
       rm = dabs(p)
c ..................................................................... 
      return
c ......................................................................
c
c ... reconstrucao da gradiente de pressao de correcao/pressao
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
c ... Green-Gauss
        grad(1,1) = 0.0d0
        grad(2,1) = 0.0d0
        do i = 1, nshared
          uf = um(i)
          grad(1,1) = grad(1,1) + uf*n(1,i)
          grad(2,1) = grad(2,1) + uf*n(2,i)
        enddo
        grad(1,1) = grad(1,1) / area(1)
        grad(2,1) = grad(2,1) / area(1)
c ......................................................................
c
c ... funcao de limatacao do fluxo para o termo advectivo 
      fluxl(1) = 1.0d0
      return
c .....................................................................
c      
c ...  
  300 continue
      return
c ......................................................................
c
c ... skewnessCorrection    
  400 continue
      p  = 0.d0
      do i = 1, nshared
        viznel = viz(i)
c ... condicao de contorno      
        if( viznel .gt. 0) then
          if(cor .eq. 1 .or. cor .eq. 2 .or. cor .eq. 3) then  
          elseif(cor .eq. 5) then  
            alpha = mKsiF(i)/mksi(i)
               kf = (1.0d0-alpha)*k(2,idCell) + alpha*k(2,i)
c ... aritimetica
            ddf(1) = (1.0d0-alpha)*d(1,idcell)  + alpha*d(1,i)
            ddf(2) = (1.0d0-alpha)*d(2,idcell)  + alpha*d(2,i)
            ndd    = 0.5d0*(ddf(1) + ddf(2))
c ... calculo dos pontos virtuais
            gama   = (ksi(1,i)*n(1,i) + ksi(2,i)*n(2,i))*mksi(i)
c ... celula central        
            aux          = df(1,i)*n(1,i)+df(2,i)*n(2,i)
            dVirtualP(1) = xm(1,i) - aux*n(1,i) - xc(1,idCell)
            dVirtualP(2) = xm(2,i) - aux*n(2,i) - xc(2,idCell)
            gf(1)        = grad(1,idCell)*dVirtualP(1)
     .                   + grad(2,idCell)*dVirtualP(2)
c ... celula vizinho        
            aux            =(xm(1,i)-xc(1,i))*n(1,i)
     .                     +(xm(2,i)-xc(2,i))*n(2,i)
            dVirtualViz(1) = xm(1,i) - aux*n(1,i) - xc(1,i)
            dVirtualViz(2) = xm(2,i) - aux*n(2,i) - xc(2,i)
            gf(2)          = grad(1,i)*dVirtualViz(1)
     .                     + grad(2,i)*dVirtualViz(2)
c ... 
            p   = p + ((kf*meta(i)*ndd)*(gf(2) -gf(1)))/gama
          endif  
        endif
      enddo
c .....................................................................      
      return
      end
c **********************************************************************
c
c **********************************************************************
c * CELL_si_p_i_ggl : Celula 2D de para equacao de correcao de pressao *
c * do metodo simple incompressivel                                    *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * a      - nao definido                                              *
c * x      - coordenadas dos vertices da celula central e do seus      *
c * vizinhos                                                           *
c * u      - temperatura por celula                                    *
c * u0     - campo de pressao estimado                                 *
c * u1     - temperatura da celula do passo anterior                   *
c * w      - campo de velociade conhecido                              *
c * k      - propriedades da celula e do seus vizinhos                 *
c * grad   - gradiente reconstruido na celula                          *
c * gradP  - gradiente de pressao reconstruido na celula               *
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
      subroutine cell_si_p_i_ggl(a,x,un,u,u0,u1,w,d,grad,gradP,iM,fluxl
     .                          ,k,rm
     .                          ,p,sedge,dt,pedge,viz,nshared,ndm,iws
     .                          ,nel,sn,acod,bs)
      implicit none
      include 'simple.fi'
c ... variavel externas
      real*8 a(*),rm,p,sedge(3,*),ap,fluxl(nshared+1)
      real*8 dt,iM(4,nshared+1),w(ndm,*),grad(ndm,*),gradP(ndm,*) 
      integer pedge(*),viz(*),iws
      integer nel,acod,nshared,ndm,sn(2,*)
c ... variavel internas
      real*8 xc(3,5),x(ndm,nshared,nshared+1),d(2,nshared+1),un(*)
      real*8 gfk,umin,umax,r,specificMass,dfd,kf,nk,ndd,cd
      real*8 um(4),kis(3,4),ksi(3,4),eta(3,4),n(3,4),meta(4),mksi(4)
      real*8 u(5),u0(5),u1(5),ca(4),k(10,*),xm(3,4),gf(2),gfp(2),area(5)
      real*8 areacell,ug(4),xcg(3,4),alpha,dd(2),vm(2,4),vmf(2),uf
      real*8 limitv,eps,df(3,4),par(2),aux,pl,ddf(2),pre(nshared+1)
      real*8 wfn,sp,cvc,cv,specificMassA,mKsiF(4),gama
      real*8 const,km(3,4),ZERO
      real*8  dVirtualP(2),dVirtualViz(2)
      integer i,j,l,idcell,viznel,icod,cor
      logical bs
      parameter (ZERO  = 1.0d-60)
      parameter (const = 1.0d60)
      external limitv,areacell
c ...
      icod   = 1
      idCell = nshared + 1
      cor   =  intVel
c ......................................................................
c
c ...
       call cellGeom2D(x,mksi,meta,ksi,eta,xc,area,n,xm,df,ca
     .                ,viz,nshared,ndm,sn,acod)
c .....................................................................
c
c ... propriedeade da celula
      specificMassA=  k(2,idCell)
      icod         =  k(7,idCell)
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
c ... condicao de controno para as celulas fantasma
          if(pedge(i) .eq. 3) then
            if(iws .eq. 2) then
              ug(i) = 0.0d0
            else    
              ug(i) = sedge(3,i)
            endif
c ... variavel prescrita prescrita 
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
c ....................................................................
c
c ... correcao de pressao  
  100 continue
      p  = 0.d0
      sp = 0.d0
      do i = 1, nshared
        viznel = viz(i)
c ... condicao de contorno      
        if( viznel .lt. 0) then
          a(i) = 0.d0
c ... velocidade prescrita
          if(pedge(i) .eq. 2 ) then
c ... fluxo advectivo de primeira ordem
            wfn = sedge(1,i)*n(1,i) + sedge(2,i)*n(2,i)
            cv  = specificMassA*wfn*meta(i)
            p   =  p  - cv
c ... pressao prescrita
          else if(pedge(i) .eq. 3) then
            ddf(1) = d(1,idcell)  
            ddf(2) = d(2,idcell) 
c ... produtos interno
            dd(1) = ddf(1)*n(1,i)
            dd(2) = ddf(2)*n(2,i)
            ndd   = dd(1)*n(1,i) + dd(2)*n(2,i)         
c ... difusao direta
            dfd  = (specificMassA*ndd*meta(i))/ca(i)     
             sp  = sp + dfd
c ... velocidade 
            wfn = w(1,idcell)*n(1,i) + w(2,idcell)*n(2,i)
            if(dabs(wfn) .gt. ZERO) then
              if(wfn .gt. 0.0) then
                cv  = specificMassA*wfn*meta(i)
              elseif( wfn .lt. 0.0) then
                cv  = specificMassA*wfn*meta(i)
              endif
            else
                cv = 0.0d0
            endif
            p     = p - cv    
c ... localmente parabolica
          else if(pedge(i) .eq. 4) then
            wfn = sedge(1,i)*n(1,i) + sedge(2,i)*n(2,i)
            cv  = specificMassA*wfn*meta(i)
            p   =  p  - cv
c .....................................................................
          endif
        else
c ... M. Darwish, I. Sraj, F. Maukalled (2009)       
          if(cor .eq. 1 .or. cor .eq. 2 .or. cor .eq. 3) then
c ... interpolacao da propriedades
            alpha = mKsiF(i)/mksi(i)
            kf = k(2,idCell)
c ... aritimetica
            ddf(1) = (1.0d0-alpha)*d(1,idcell)  + alpha*d(1,i) 
            ddf(2) = (1.0d0-alpha)*d(2,idcell)  + alpha*d(2,i)
c ... produtos interno
            dd(1) = ddf(1)*n(1,i)
            dd(2) = ddf(2)*n(2,i)
            ndd   = dd(1)*n(1,i)    + dd(2)*n(2,i)         
            nk    = n(1,i)*ksi(1,i) + n(2,i)*ksi(2,i)
c ... difusao direta
            dfd   = (kf*meta(i)*ndd)/(nk*mksi(i))
            a(i)  = dfd
c ... Ferziger Peric            
          else if(cor .eq. 5) then
c ... interpolacao da propriedades
            alpha = mKsiF(i)/mksi(i)
            kf = (1.0d0-alpha)*k(2,idCell) + alpha*k(2,i)
c ... aritimetica
            ddf(1) = (1.0d0-alpha)*d(1,idcell)  + alpha*d(1,i) 
            ddf(2) = (1.0d0-alpha)*d(2,idcell)  + alpha*d(2,i)
c ... produtos interno
            ndd   = 0.5d0*(ddf(1) + ddf(2))         
            nk    = n(1,i)*ksi(1,i) + n(2,i)*ksi(2,i)            
c ... difusao direta
            dfd   = (kf*meta(i)*ndd)/(nk*mksi(i))
            a(i)  = dfd
          endif
c ... interpolacao linear das velociadades
          call cellVel(wfn,w,gradP,ddf,d,iM,ksi,n,u1,mksi,meta,area
     .                ,alpha,ndm,nshared,intVel,i)
c ...
          p   = p - kf*wfn
        endif
      enddo
c .....................................................................
c
c ... fixa o valor da pressao de correcao igual a zero
      if( pedge(idCell) .eq. -1) sp = sp + const
c .....................................................................
c
c ...
      a(idcell) = sp
      do i = 1, nshared
        a(idcell) = a(idcell) + a(i)
      enddo
c .....................................................................
c
c .....................................................................
c
c ... residou de massa por celula
      rm = p
c ..................................................................... 
      return
c ......................................................................
c
c ... reconstrucao da gradiente de pressao 
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
c .....................................................................
c
c ... funcao de limatacao do fluxo para o termo advectivo 
      fluxl(1) = 1.0d0
      return
c .....................................................................
c
c ... equacao da pressao                   
  300 continue
c     p  = 0.d0
c     sp = 0.d0
c     do i = 1, nshared
c       viznel = viz(i)
c ... condicao de contorno      
c       if( viznel .lt. 0) then
c         a(i) = 0.d0
c ... velocidade prescrita
c         if(pedge(i) .eq. 2 ) then
c ... fluxo advectivo de primeira ordem
c           wfn = sedge(1,i)*n(1,i) + sedge(2,i)*n(2,i)
c           cv  = specificMassA*wfn*meta(i)
c           p   =  p  - cv
c ... pressao prescrita
c         else if(pedge(i) .eq. 3) then
c           ddf(1) = d(1,idcell)  
c           ddf(2) = d(2,idcell) 
c ... produtos interno
c           dd(1) = ddf(1)*n(1,i)
c           dd(2) = ddf(2)*n(2,i)
c           ndd   = dd(1)*n(1,i) + dd(2)*n(2,i)         
c ... difusao direta
c           dfd  = (specificMassA*ndd*meta(i))/ca(i)     
c            sp  = sp + dfd
c             p  = p  + dfd*sedge(3,i)
c ... velocidade 
c           wfn = w(1,idcell)*n(1,i) + w(2,idcell)*n(2,i)
c           if(dabs(wfn) .gt. ZERO) then
c             if(wfn .gt. 0.0) then
c               cv  = specificMassA*wfn*meta(i)
c             elseif( wfn .lt. 0.0) then
c               cv  = specificMassA*wfn*meta(i)
c             endif
c           else
c               cv = 0.0d0
c           endif
c           p     = p - cv    
c ... localmente parabolica
c         else if(pedge(i) .eq. 4) then
c           wfn = sedge(1,i)*n(1,i) + sedge(2,i)*n(2,i)
c           cv  = specificMassA*wfn*meta(i)
c           p   =  p  - cv
c .....................................................................
c         endif
c       else
c ... interpolacao da propriedades
c         alpha = mKsiF(i)/mksi(i)
c         kf    = (1.0d0-alpha)*k(2,idCell) + alpha*k(2,i)
c ... aritimetica
c         ddf(1) = (1.0d0-alpha)*d(1,idCell)  + alpha*d(1,i) 
c         ddf(2) = (1.0d0-alpha)*d(2,idCell)  + alpha*d(2,i)
c ... produtos interno
c         dd(1) = ddf(1)*n(1,i)
c         dd(2) = ddf(2)*n(2,i)
c         ndd   = dd(1)*n(1,i)    + dd(2)*n(2,i)         
c         nk    = n(1,i)*ksi(1,i) + n(2,i)*ksi(2,i)
c ... difusao direta
c         dfd   = (kf*meta(i)*ndd)/(nk*mksi(i))
c         a(i)  = dfd
c ... interpolacao linear das velociadades
c         call cellVel(wfn,w,gradP,ddf,d,iM,ksi,n,u1,mksi,area,alpha,ndm
c    .                ,nshared,4,i)
c ...
c         p   = p - kf*wfn*meta(i)
c       endif
c     enddo
c .....................................................................
c
c ... fixa o valor da pressao de correcao igual a zero
c     if( pedge(idCell) .eq. -1) sp = sp + const
c .....................................................................
c
c ...
c     a(idcell) = sp
c     do i = 1, nshared
c       a(idcell) = a(idcell) + a(i)
c     enddo
c .....................................................................
c
c ... under-relaxation(simple)
c     a(idCell) = a(idCell)/underP 
c     p         = p + (1-underP)*a(idCell)*u(idCell)
c .....................................................................
c
c ... residou da celula
c     rm  = p - a(idCell)*u(idCell)
c     do i = 1, nshared
c       rm  = rm  + a(i)*u(i)
c     enddo
c .....................................................................
c
c ...      
c     return
c .....................................................................
c  
c ... skewnessCorrection    
  400 continue
      p  = 0.d0
      do i = 1, nshared
        viznel = viz(i)
c ... condicao de contorno      
        if( viznel .gt. 0) then
          if(cor .eq. 1 .or. cor .eq. 2 .or. cor .eq. 3) then
          elseif(cor .eq. 5) then  
            alpha = mKsiF(i)/mksi(i)
               kf = (1.0d0-alpha)*k(2,idCell) + alpha*k(2,i)
c ... aritimetica
            ddf(1) = (1.0d0-alpha)*d(1,idcell)  + alpha*d(1,i)
            ddf(2) = (1.0d0-alpha)*d(2,idcell)  + alpha*d(2,i)
            ndd    = 0.5d0*(ddf(1) + ddf(2))
c ... calculo dos pontos virtuais
            gama   = (ksi(1,i)*n(1,i) + ksi(2,i)*n(2,i))*mksi(i)
c ... celula central        
            aux          = df(1,i)*n(1,i)+df(2,i)*n(2,i)
            dVirtualP(1) = xm(1,i) - aux*n(1,i) - xc(1,idCell)
            dVirtualP(2) = xm(2,i) - aux*n(2,i) - xc(2,idCell)
            gf(1)        = grad(1,idCell)*dVirtualP(1)
     .                   + grad(2,idCell)*dVirtualP(2)
c ... celula vizinho        
            aux            =(xm(1,i)-xc(1,i))*n(1,i)
     .                     +(xm(2,i)-xc(2,i))*n(2,i)
            dVirtualViz(1) = xm(1,i) - aux*n(1,i) - xc(1,i)
            dVirtualViz(2) = xm(2,i) - aux*n(2,i) - xc(2,i)
            gf(2)          = grad(1,i)*dVirtualViz(1)
     .                     + grad(2,i)*dVirtualViz(2)
c ... 
            p   = p + ((kf*meta(i)*ndd)*(gf(2) -gf(1)))/gama
          endif  
        endif
      enddo
c .....................................................................
c
c ...
      return  
      end
c **********************************************************************