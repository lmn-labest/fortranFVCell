c *********************************************************************
c * CELLLIB: biblioteca de celulas                                    *
c * ------------------------------------------------------------------*
c *  a      - nao definido                                            *
c *  x      - coordenadas nodais                                      *
c *  un     - volores nodais                                          *
c *  u      - valor por celula                                        *
c *  u0     - valor por celula                                        *
c *  u1     - valor por celula                                        *
c *  ro     - massa especifica da celula e sua vizinhas em t e t+1    *
c *  grad   - gradiente nas celulas                                   *
c *  grad1  - gradiente nas celulas                                   *
c *  grad2  - gradiente nas celulas                                   *
c *  div    - divergente da velocidade                                *
c *  iM     - interpola para velocidades nas faces                    *
c *  fluxl  - limitador de fluxo na celula                            *
c *  w      - velocidade nas celulas                                  *
c *  d      - campo d do simple                                       *
c *  k      - propriededas da celula e sua vizinhas                   *
c *  sp     - nao definido                                            *
c *   p     - nao definido                                            *
c * sedge   - condicoes de contorno nas arestas                       *
c * pedge   - tipo de condicao de contorno                            *
c * viz     - vizinhos                                                *
c * ndm     - numero de dimensoes                                     *
c * nen     - numero de nos por celula                                *
c * nshared - numero de celula compartilhados pela celula central     *
c * type    - tipo de elemento                                        *
c * iws     - 1 - motagem do sistema (KU=F)                           *
c *           2 - calculo do gradiente                                *
c * iws1    - 1 - direcao 1                                           *
c *           2 - direcao 2                                           *
c * lib     -                                                         *
c * nel     - numero da celula                                        *
c * dt      - passo de tempo                                          *
c * alpha   - paramentro alpha                                        *
c * mP      - parametros da equacao de momentos(cfl,Reynalds)         *
c * eddyVisc- viscosidade tubulenca                                   *
c * bs      - Euller Backward de segunda ordem                        *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c *  a     - valores fora da diagonal de coeficientes                 *
c *  p     - residuo                                                  *
c * sp     - residou por elemento ( bi-(Ai)*u(viz),Aiu(i))            *
c *  iM    - interpola para velocidades nas faces                     *
c * ------------------------------------------------------------------*
c *********************************************************************
      subroutine celllib(a    ,x       ,un     ,u      ,u0   
     .                  ,u1   ,ro      ,grad   ,grad1  ,grad2
     .                  ,div  ,iM      ,fluxl  ,w      ,d    
     .                  ,k    ,sp      ,p      ,sedge  ,pedge
     .                  ,viz  ,lls     ,nen    ,nshared,ndm  
     .                  ,type ,iws    ,iws1    ,lib    ,nel  
     .                  ,dt   ,alpha  ,mP      ,eddyVisc,bs)
      implicit none
      include 'error.fi'
      real*8 a(*),x(*)
      real*8 u(*),u0(*),u1(*)
      real*8 grad(*),grad1(*),grad2(*),k(*)
      real*8 lls(*),w(*),dt,alpha,fluxl(*),d,sp(*),p(*),sedge(*),div(*)
      real*8 iM(*),Mp(*),ro(*),un(*),eddyVisc(*)
      integer pedge(*),viz(*),nen,nshared,ndm,type,iws,iws1,nel
      integer acod,sn(2,4),lib
      logical bs
c ... variaveis da celula que dependem do numero de faces     
      call setcell(type,sn,acod)
c .....................................................................
c
c ... cell 2D - difusao (GGl)
      if(type .eq. 3 .or. type .eq. 5) then
c        call cell_d_ggl(a,x,u,u0,grad,fluxl,k,sp,p,sedge,dt,pedge,viz
c     .                 ,nshared,ndm,iws,nel,alpha,sn,acod,bs) 
c ......................................................................
c
c ... cell 2D - difusao (WLS)
      elseif(type .eq. 7) then 
c        call cell_d_wls(a,x,u,u0,grad,fluxl,k,sp,p,sedge,dt,pedge,viz
c     .                 ,lls,nshared,ndm,iws,nel,alpha,sn,acod,bs)
c ......................................................................      
c
c ... cell 2D advecao-difusao (Edge-based-GGl-TVD)
      elseif(type .eq. 13 .or. type .eq. 15) then 
c        call cell_ad_eb_ggl(a,x,u,u0,w,grad,fluxl,k,sp,p,sedge
c     .                     ,dt,pedge,viz,nshared,ndm,iws,nel,sn,alpha
c     .                     ,acod,bs)
c ......................................................................
c
c ... cell 2D advecao-difusao (advecao- Volume-based-GGl-SOU - difusao)
      elseif(type .eq. 23 .or. type .eq. 25) then 
c        call cell_ad_vb_ggl(a,x,u,u0,w,grad,fluxl,k,sp,p,sedge
c     .                     ,dt,pedge,viz,nshared,ndm,iws,nel,alpha
c     .                     ,sn,acod,bs)
c ......................................................................
c
c ... cell 2D adveccao pura(Edge-based-GGl-TVD)
      elseif(type .eq.  33 .or. type .eq. 35) then 
        call cell_a_eb_ggl(a,x,u,u0,w,grad,fluxl,k,sp,p
     .       ,sedge,dt,pedge,viz,nshared,ndm,iws,nel,alpha,sn,acod,bs)
c ......................................................................
c
c ... cell 2D adveccao pura(Volume-based-GGl-SOU)
      elseif(type .eq.  43 .or. type .eq. 45) then 
        call cell_a_vb_ggl(a,x,u,u0,w,grad,fluxl,k,sp,p,sedge,dt,pedge
     .                     ,viz,nshared,ndm,iws,nel,alpha,sn,acod,bs)
c .....................................................................
c
c ... cell 2D simple - ideal gas (Edge-based-GGl-TVD)
      elseif(type .eq.  53 .or. type .eq. 55) then
        if( lib .eq. 1) then
          call cell_si_eb_idg_ggl(a,x,u,u0,u1,rO,w,d,grad,grad1,grad2
     .                           ,div,iM,fluxl,k,sp,p,sedge,dt
     .                       ,pedge,viz,nshared,ndm,iws,iws1,nel,sn,acod
     .                       ,Mp,bs)
        elseif(lib .eq.2) then
          call cell_si_p_idg_ggl(a,x,u,u0,u1,rO,w,d,grad,grad1,iM,fluxl
     .                          ,k,sp,p,sedge,dt,pedge,viz,nshared,ndm
     .                          ,iws,iws1,nel,sn,acod,bs)
        elseif(lib .eq.3) then
          call cell_e_eb_idg_ggl(a,x,u,u0,rO,w,grad,fluxl,k,sp,p,sedge
     .                      ,dt,pedge,viz,nshared,ndm,iws,nel,sn,acod
     .                      ,bs)
        endif 
c .....................................................................
c
c ... cell 2D simple - incompressivel (Edge-based-GGl-TVD)
      elseif(type .eq.  63 .or. type .eq. 65) then
        if( lib .eq. 1) then
          call cell_si_eb_i_ggl(a,x,u,u0,u1,w,d,grad,grad1,grad2
     .                         ,iM,fluxl,k,sp,p,sedge,dt
     .                         ,pedge,viz,nshared,ndm,iws,iws1,nel,sn
     .                         ,acod,Mp,bs)
        elseif(lib .eq.2) then
          call cell_si_p_i_ggl(a,x,un,u,u0,u1,w,d,grad,grad1,iM,fluxl,k
     .                        ,sp,p,sedge,dt,pedge,viz,nshared,ndm,iws
     .                        ,nel,sn,acod,bs)
        elseif(lib .eq. 4) then
          call cell_t_eb_ggl(a,x,u,u0,rO,w,grad,fluxl,k,sp,p,sedge
     .                      ,dt,pedge,viz,nshared,ndm,iws,nel,sn,acod
     .                      ,bs)
        endif 
c .....................................................................
c
c ... cell 2D simple - incompressivel-LES (Wall-Model-Edge-based-GGl-TVD)
      elseif(type .eq.  73 .or. type .eq. 75) then
        if( lib .eq. 1) then
          call cell_si_wm_les_eb_i_ggl(a,x,u,u0,u1,w,d,grad,grad1,grad2
     .                                ,iM,fluxl,k,sp,p,sedge,dt
     .                                ,pedge,viz,nshared,ndm,iws,iws1
     .                                ,nel,sn,acod,Mp,eddyVisc,bs)
        elseif(lib .eq.2) then
          call cell_si_p_i_ggl(a,x,un,u,u0,u1,w,d,grad,grad1,iM,fluxl,k
     .                        ,sp,p,sedge,dt,pedge,viz,nshared,ndm,iws
     .                        ,nel,sn,acod,bs)
        elseif(lib .eq. 4) then
        endif 
c .....................................................................
      else
        write(*,2000)'celllib','Celllib.f',
     .               'Tipo de celula nao existente'
        stop
      endif 
      return
      end
c **********************************************************************
c
c *********************************************************************
c * BOUNDQUAD: transforma condicoes de contorno nodais para arestas    *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * lpn    - restricoes nodais (1-temperatura;3-fluxo;4-fluxo e temp   *
c * lpe    - nao definido                                              *
c * ltn    - temperatura por no                                        *
c * lfn    - fluxo por no                                              *
c * lse    - nao definido                                              *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * lse    - valor das condicoes de contorno                           *
c * lpe    - arestas com restricoes de contorno (1-temp;0-flux)        *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine boundquad(lpn,lpe,ltn,lfn,lse) 
      implicit none
      integer lpn(*),lpe(*)
      real*8  ltn(*),lfn(*),lse(*)
      integer isnod(2,4),i,no1,no2,aux
      data isnod/1,2,2,3,3,4,4,1/
c      print*,'nnode'
c      print*,lpn(1),lpn(2),lpn(3),lpn(4)
c      write(*,'(4f6.2)')ltn(1),ltn(2),ltn(3),ltn(4)
c      write(*,'(4f6.2)')lfn(1),lfn(2),lfn(3),lfn(4)
c ... passando as restricoes nodais para as arestas     
      do i = 1, 4
        no1 = isnod(1,i)
        no2 = isnod(2,i)
        aux = lpn(no1) + lpn(no2)
        if( (aux .eq. 2) .or. (aux .eq. 5) ) then
c          print*,ltn(no1),ltn(no2)
          lpe(i) = 1
          lse(i) = 0.5d0*(ltn(no1) + ltn(no2)) 
        endif
        if( (aux .eq. 7) .or. (aux .eq. 6) ) then
          lse(i) = 0.5d0*(lfn(no1) + lfn(no2))
        endif
      enddo
c .....................................................................
c
c ...
c      print*,'el'
c      print*,lpe(1),lpe(2),lpe(3),lpe(4)
c      write(*,'(4f6.2)')lse(1),lse(2),lse(3),lse(4)
c      print*,'-------------------------------------'
      return
      end
c *********************************************************************
c
c *********************************************************************
c * BOUNDTRIA: transforma condicoes de contorno nodais para arestas    *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * lpn    - restricoes nodais (1-temperatura;3-fluxo;4-fluxo e temp   *
c * lpe    - nao definido                                              *
c * ltn    - temperatura por no                                        *
c * lfn    - fluxo por no                                              *
c * lse    - nao definido                                              *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * lse    - valor das condicoes de contorno                           *
c * lpe    - arestas com restricoes de contorno (1-temp;0-flux)        *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine boundtria(lpn,lpe,ltn,lfn,lse) 
      implicit none
      integer lpn(*),lpe(*)
      real*8  ltn(*),lfn(*),lse(*)
      integer isnod(2,3),i,no1,no2,aux
      data isnod/1,2,2,3,3,1/
c      print*,'nnode'
c      print*,lpn(1),lpn(2),lpn(3),lpn(4)
c      write(*,'(4f6.2)')ltn(1),ltn(2),ltn(3),ltn(4)
c      write(*,'(4f6.2)')lfn(1),lfn(2),lfn(3),lfn(4)
c ... passando as restricoes nodais para as arestas     
      do i = 1, 3
        no1 = isnod(1,i)
        no2 = isnod(2,i)
        aux = lpn(no1) + lpn(no2)
        if( (aux .eq. 2) .or. (aux .eq. 5) ) then
c          print*,ltn(no1),ltn(no2)
          lpe(i) = 1
          lse(i) = 0.5d0*(ltn(no1) + ltn(no2)) 
        endif
        if( (aux .eq. 7) .or. (aux .eq. 6) ) then
          lse(i) = 0.5d0*(lfn(no1) + lfn(no2))
        endif
      enddo
c .....................................................................
c
c ...
c      print*,'el'
c      print*,lpe(1),lpe(2),lpe(3),lpe(4)
c      write(*,'(4f6.2)')lse(1),lse(2),lse(3),lse(4)
c      print*,'-------------------------------------'
      return
      end
c *********************************************************************
c
c **********************************************************************
c * AREATRIACELL: Calculo da area do triangulo                         *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * eta    - vetor nas aresta                                          *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * -------------------------------------------------------------------*
c **********************************************************************
      real*8 function areatriacell(eta)
      implicit none
      real*8 eta(3,3),c(3),v1(3),v2(3),a,dot_local
c ...
      v1(1) = eta(1,1)
      v1(2) = eta(2,1)
      v1(3) = eta(3,1)
      v2(1) = eta(1,2)
      v2(2) = eta(2,2)
      v2(3) = eta(3,2)
      call vet(v1,v2,c)
      a = 0.5d0*dsqrt(dot_local(c,c,3))
c ......................................................................
      areatriacell = a
      return
      end
c **********************************************************************
c * AREAQUADCELL: Calculo da area do quadrilatero                      *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * eta    - vetor nas aresta                                          *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * -------------------------------------------------------------------*
c **********************************************************************
      real*8 function areaquadcell(eta)
      implicit none
      real*8 eta(3,4),c(3),v1(3),v2(3),a,dot_local
c ...      
      v1(1) = eta(1,1)
      v1(2) = eta(2,1)
      v1(3) = eta(3,1)
      v2(1) = eta(1,2)
      v2(2) = eta(2,2)
      v2(3) = eta(3,2)
      call vet(v1,v2,c)
      a = 0.5d0*dsqrt(dot_local(c,c,3))
c .....................................................................      
c
c ...      
      v1(1) = eta(1,3)
      v1(2) = eta(2,3)
      v1(3) = eta(3,3)
      v2(1) = eta(1,4)
      v2(2) = eta(2,4)
      v2(3) = eta(3,4)
      call vet(v1,v2,c)
      a = a + 0.5d0*dsqrt(dot_local(c,c,3))
c .....................................................................      
c
c ...      
      v1(1) = eta(1,1)
      v1(2) = eta(2,1)
      v1(3) = eta(3,1)
      v2(1) = eta(1,4)
      v2(2) = eta(2,4)
      v2(3) = eta(3,4)
      call vet(v1,v2,c)
      a = a + 0.5d0*dsqrt(dot_local(c,c,3))
c .....................................................................      
c
c ...      
      v1(1) = eta(1,2)
      v1(2) = eta(2,2)
      v1(3) = eta(3,2)
      v2(1) = eta(1,3)
      v2(2) = eta(2,3)
      v2(3) = eta(3,3)
      call vet(v1,v2,c)
      a = a + 0.5d0*dsqrt(dot_local(c,c,3))
c .....................................................................      
c
c ...
      areaquadcell = 0.5d0*a
      return
      end
c *********************************************************************
c
c *********************************************************************
c     subroutine errosol(u,error,x,ix,numel,nen,ndf,ndm,cod)
c     implicit none
c     logical cod
c     real*8 u(ndf,*),error(*),du(3),f(3),dot,x(ndm,*),xc(3)
c     integer ix(nen+1,*),numel,nen,ndf,ndm,nel,j,no,k
c     do nel = 1, numel
c        do k = 1, ndm
c          xc(k) = 0.d0
c        enddo  
c ... calculo do centroide            
c         do j = 1, nen
c           no = ix(j,nel)
c           do k = 1, ndm 
c             xc(k) = xc(k) + x(k,no)
c           enddo
c         enddo
c ... Triangulo          
c         if( nen .eq. 3) then
c           do j = 1, ndm 
c             xc(j) = xc(j)/3.0d0
c           enddo
c         endif
c ... Quadrilatero          
c         if( nen .eq. 4) then
c           do j = 1, ndm 
c             xc(j) = xc(j)*0.25d0
c           enddo
c         endif  
c .....................................................................    
c       call exectsol(f,xc,cod)
c       do j = 1, ndf
c         du(j) = u(j,nel) - f(j)
c       enddo
c       error(nel) = dsqrt(dot(du,du,ndf))/dsqrt(dot(f,f,ndf))
c     enddo
c     return
c     end
c *********************************************************************
c 
c *********************************************************************
      subroutine exectsol(f,x,cod)
      implicit none
      logical cod
      real*8 x(3),f(3),a,r
      r    = dsqrt(x(1)**2+x(2)**2)
c ... solucao      
      if(cod) then
        f(1) = 500.d0/dlog(0.1d0)*dlog(r)
c ... gradiente        
      else
        a    = 500.d0/dlog(0.1d0)*(1.0d0/r)
        f(1) = a*x(1)/r
        f(2) = a*x(2)/r
      endif
      return
      end
c *********************************************************************
c
c **********************************************************************
c * LIMIT: funcao limitadora para o termo convectivo                   *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * par  - para metros utilizados                                      *
c * icod - tipo da funcao limitadora                                   *
c **********************************************************************
      real*8 function limit(par,icod)
      implicit none
      real*8 par(*),r,a,b,c,eps
      integer icod
c ... UD - upwind de primeira ordem
      if(icod.eq.1) then
        limit = 0.0d0
c .....................................................................
c
c ... LUD - upwind de segunda ordem
      elseif(icod.eq.2) then
        limit = 1.0d0
c .....................................................................
c
c ... Van Leer - TVD
      elseif(icod.eq.3) then
c ... razao entre o fluxo upwind e downwind
        r = par(1)
        limit = (r + dabs(r)) / (1.0d0 + dabs(r))
c .....................................................................
c
c ... Van Albada - TVD
      elseif(icod.eq.4) then
c ... razao entre o fluxo upwind e downwind
        r = par(1)
        limit = 0.0d0
        if( r .gt. 0.d0 ) then
          limit = (r + r*r) / (1.0d0 + r*r)
        endif
c .....................................................................
c
c ... Mid - Mod - TVD
      elseif(icod.eq.5) then
c ... razao entre o fluxo upwind e downwind
        r = par(1)
        limit = max(0.0d0,min(r,1.0d0))
c .....................................................................
c
c ... OSHER - TVD
      elseif(icod.eq.6) then
c ... razao entre o fluxo upwind e downwind
        r = par(1)
        limit = max(0.0d0,min(2.0d0,r))
c .....................................................................
c
c ... SUPERBEE - TVD
      elseif(icod.eq.7) then
c ... razao entre o fluxo upwind e downwind
        r     = par(1)
        a     = min(2.0d0*r,1.0d0)
        b     = min(r,2.0d0)
        c     = max(a,b) 
        limit = max(0.0d0,c) 
c .....................................................................
c
c ... Albada - Van Leer
      elseif(icod.eq.8) then
c ... gradiente upwind 
        a = par(1)
c ... gradiente central
        b = par(2)
c ... parametro de controle 
        eps  = par(3)
c ...
        limit = 0.0d0
        if(a*b .gt. 0.0d0) then 
          limit = ((a*a+eps)*b+(b*b+eps)*a)/(a*a+b*b+2*eps)
        endif
c .....................................................................
c
c ...
      else
        print*,'Funcao limitadora nao existente!!'
        stop
      endif
c ......................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
c * LIMITV: funcao limitadora para o termo convectivo volume           *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * r    - x                                                           *
c * par  - parametros                                                  *
c * icod - tipo da funcao limitadora                                   *
c **********************************************************************
      real*8 function limitv(r,par,icod)
      implicit none
      real*8 r,par(*),e
      integer icod
c ... UP                        
      if(icod.eq.1) then
        limitv = 0.0d0
c ... SOU
      elseif(icod.eq.2) then
        limitv= 1.d0
c ...  Barth                                
      elseif(icod.eq.3) then
        limitv= min(1.d0,r)
c ...  modified Batth                                
      elseif(icod.eq.4) then
        e     = par(1)*par(2)
        e     = e*e*e
        limitv= (r*r+2*r+e)/(r*r+r+2+e)
c ...                   
      else
        print*,'Funcao limitadora nao existente!!'
        stop
      endif
      return
      end
c **********************************************************************
c
c **********************************************************************
c * AREACELL: Calulo da area de um poligono                            *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * eta    - vetor nas aresta                                          *
c * icod   - 3 triagulos                                               *
c *          4 quadrilateros                                           *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * retorna o valor da area do poligono                                *
c * -------------------------------------------------------------------*
c **********************************************************************
      real*8 function areacell(eta,icod)
      implicit none
      integer icod
      real*8 eta(*),area,areatriacell,areaquadcell
      external areatriacell,areaquadcell
      area = 0.0d0
      if(icod .eq. 3) then
        area = areatriacell(eta)
      elseif(icod .eq. 4) then
        area = areaquadcell(eta)
      endif
      areacell = area
      return
      end
c **********************************************************************
c
c **********************************************************************
c * SETCELL: varais que depende do numero de faces                     *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * type   - tipo de celula                                            *
c * sn(2,i)- nos da aresta i                                           *
c * acod   - 3 triagulos                                               *
c *          4 quadrilateros                                           *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * acod   -                                                           *
c * sn     -                                                           *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine setcell(type,sn,acod)
      implicit none
      integer acod,type,sn(2,*)
c ... celulas triagulares
      if(type .eq. 3  .or. type .eq. 13 .or. type .eq. 23 .or.
     .   type .eq. 33 .or. type .eq. 43 .or. type .eq. 53 .or. 
     .   type .eq. 63 .or. type .eq. 73) then
        sn(1,1) = 1
        sn(2,1) = 2
        sn(1,2) = 2
        sn(2,2) = 3
        sn(1,3) = 3
        sn(2,3) = 1
           acod = 3
c ... celulas quadrilateras
      elseif(type .eq. 5  .or. type .eq. 7  .or. type .eq. 15 .or. 
     .       type .eq. 25 .or. type .eq. 35 .or. type .eq. 45 .or. 
     .       type .eq. 55 .or. type .eq. 65 .or. type .eq. 75) then
        sn(1,1) = 1
        sn(2,1) = 2
        sn(1,2) = 2
        sn(2,2) = 3
        sn(1,3) = 3
        sn(2,3) = 4
        sn(1,4) = 4
        sn(2,4) = 1
        acod    = 4
      endif
      return
      end
c **********************************************************************
c
c **********************************************************************
c * VECTORKM2D: vetor que une a intersecão da reta entre os dois       *
c * centroides ao ponto médio da aresta                                *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * xcg    - coordenadas dos centroides das celulas vizinha/ fantasma  *
c * xc     - coordenadas dos centroides das celulas vizinhas e central *
c * x      - coordenadas dos vertices da celula central e do seus      *
c * xm     - coordenadas do ponto medio das arestas da celula central  *
c * km     -                                                           *
c * mkm    -                                                           *
c * sn(2,i)- nos da aresta i                                           *
c * nshared- numeros de faces por celulas                              *
c * ndm    - numeros de dimesoes                                       *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * km     - vetor que une a intersecão da reta entre os dois          *
c * centroides ao ponto médio da aresta                                *
c * mkm    - distancia do ksi ate a face                               *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine vectorKm2d(xcg,xc,x,xm,km,mkm,sn,nshared,ndm)
      implicit none
      real*8 a1,a2,b1,b2,c1,c2,det
      real*8 xcg(3,*),xc(3,*),x(ndm,nshared,nshared+1),km(3,*),xm(3,*)
      real*8 a(2,2),xi(2),f(2),mkm(*),ksif(2,4)
      integer i,idcell,ndm,nshared,no1,no2,sn(2,*)
      idcell = nshared+1
c ... reta 1 a1x + b1y + c1 = 0
c     reta 2 a2x + b2y + c2 = 0
c
c     | a1  b1 ||x|=|-c1| 
c     | a2  b2 ||y| |-c2|
c .....................................................................
c     a1 = i e b1 = idcell
c     a2 e b2 pontos da aresta da idcell
      do i = 1, nshared
c ... ya - yb
        a1 = xcg(2,i) - xc(2,idcell)
c ... xb - xa
        b1 = xc(1,idcell) - xcg(1,i)
c ... xayb - xbya
        c1 = xcg(1,i)*xc(2,idcell) - xc(1,idcell)*xcg(2,i)
c .....................................................................
        no1= sn(1,i)
        no2= sn(2,i)
c ... ya - yb
        a2 = x(2,no1,idcell) - x(2,no2,idcell)
c ... xb - xa
        b2 = x(1,no2,idcell) - x(1,no1,idcell)
c ... xayb - xbya
        c2 = x(1,no1,idcell)*x(2,no2,idcell) 
     .     - x(1,no2,idcell)*x(2,no1,idcell)
c ... verifica se ha interseco entre as dua retas
        det = a1*b2 - a2*b1
        if(det .eq. 0.0d0) then
          print*, 'As retas nao sao concorrentes'
          stop
        endif
c ... matrix inversa A
        a(1,1) =  (1.0d0/det)*b2
        a(1,2) = -(1.0d0/det)*b1
        a(2,1) = -(1.0d0/det)*a2
        a(2,2) =  (1.0d0/det)*a1
c ... f
        f(1) = -c1
        f(2) = -c2
c ... x=A-1f
        xi(1) = a(1,1)*f(1) + a(1,2)*f(2)
        xi(2) = a(2,1)*f(1) + a(2,2)*f(2)
c ... vetor km
        km(1,i) = xm(1,i) - xi(1)
        km(2,i) = xm(2,i) - xi(2) 
        km(3,i) = 0.0d0
c ... distancia da aresta reta que une os centroides ate a aresta        
        ksif(1,i) = xi(1) - xc(1,idcell)
        ksif(2,i) = xi(2) - xc(2,idcell)
        mkm(i)    = dsqrt(ksif(1,i)*ksif(1,i) + ksif(2,i)*ksif(2,i)) 
      enddo
      return
      end
c **********************************************************************
c
c **********************************************************************
c * CELLGEOM2D: calcula informarcoes geometricas basicas da molecula   *
c * computacional 2D                                                   *
c * -------------------------------------------------------------------*
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * x      - coordenadas dos vertices da celula central e do seus viz  *
c * mksi   - indefinido                                                *
c * meta   - indefinido                                                *
c * ksi    - indefinifo                                                *
c * eta    - indefinido                                                *
c * xc     - indefinado                                                *
c * area   - indefinido                                                *
c * n      - indefindo                                                 *
c * df     - indefinido                                                *
c * ca     - indefinido                                                *
c * viz    - vizinhos da celula                                        *
c * nshared- numeros de faces por celulas                              *
c * ndm    - numeros de dimesoes                                       *
c * sn(2,i)- nos da aresta i                                           *
c * acod   - codigo para o calculo da area                             *
c *         ( 3 - triangulo; 4 - quadrilatero)                         *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * mksi   - distancia entre os centroide                              *
c * meta   - comprimento das arestas da celula central                 *
c * ksi    - vetor que une os centroisdes(Central -> Vizinho)          *
c * eta    - vetors paralelo as arestas da celula central              *
c * xc     - coordenadas dos centroides das celulas vizinhas e central *
c * area   - area das celulas                                          *
c * n      - vetor normal a aresta da celula central                   *
c * xm     - coordenadas do ponto medio das arestas da celula central  *
c * df     - vetor que une o centroide ao centro da aresta             *
c * ca     - menor distancia entre o centroide e a aresta do lado i    *
c * -------------------------------------------------------------------*
c **********************************************************************
      subroutine cellGeom2D(x,mksi,meta,ksi,eta,xc,area,n,xm,df,ca
     .                     ,viz,nshared,ndm,sn,acod)
      implicit none
      real*8 areacell
      real*8 x(ndm,nshared,nshared+1),xc(3,*)
      real*8 mksi(*),meta(*),ksi(3,*),eta(3,*)
      real*8 n(3,*),xm(3,*),df(3,*),ca(*),area(*)
      integer i,j,l,idcell,acod,sn(2,*),ndm,nshared,viz(*)
c ...
      idcell = nshared + 1
c ...      
      do i = 1, nshared
        mksi(i) = 0.d0
        meta(i) = 0.d0
         ksi(1,i) = 0.d0
         ksi(2,i) = 0.d0
         ksi(3,i) = 0.d0
         eta(1,i) = 0.d0
         eta(2,i) = 0.d0
         eta(3,i) = 0.d0
      enddo
c      do i = 1, nshared + 1
c        area(i)  = 0.d0
c        xc(1,i)  = 0.d0
c        xc(2,i)  = 0.d0
c        xc(3,i)  = 0.d0 
c      enddo
c .....................................................................
c
c ... calculo do centro geometrico das celulas
      do i = 1, nshared + 1
        do j = 1, ndm
          xc(j,i) = 0.0d0
          do l = 1, nshared
            xc(j,i) = xc(j,i) + x(j,l,i)
          enddo
          xc(j,i) = (1.0d0/nshared) * xc(j,i)
        enddo
      enddo
c .....................................................................
c
c ... areas das celulas
       do j = 1, nshared + 1
         do i = 1, ndm
           do l = 1, nshared
c ... aresta l     
             eta(i,l) = x(i,sn(2,l),j) - x(i,sn(1,l),j)    
           enddo
         enddo     
         area(j) = areacell(eta,acod)
      enddo
c .....................................................................
c
c ... vetor que une os centroides (ksi)
c
c ...
      do i = 1, nshared
        if(viz(i) .gt. 0) then
          do j = 1, ndm
            ksi(j,i) = xc(j,i) - xc(j,idcell)   
            mksi(i)  = ksi(j,i)*ksi(j,i) + mksi(i)
          enddo
          mksi(i)    = dsqrt(mksi(i))
        endif  
      enddo
c ... distancia do vetor que une os vertices da celula central (eta) 
      do i = 1, nshared
        do j = 1, ndm
          meta(i) = meta(i) + eta(j,i)*eta(j,i)
        enddo
        meta(i)    = dsqrt(meta(i))
      enddo
c .....................................................................
c
c ... vetores unitarios
      do j = 1, ndm
        do i = 1, nshared
          eta(j,i) = eta(j,i)/meta(i)
          if(viz(i) .gt. 0) ksi(j,i) = ksi(j,i)/mksi(i)
        enddo  
      enddo  
c ... vetor normal a arestas
      do l = 1, nshared
        n(1,l) = eta(2,l)
        n(2,l) =-eta(1,l)
      enddo
c .....................................................................      
c
c ... pontos medios das aresta
      do i = 1, nshared
        xm(1,i) = 0.5d0*(x(1,sn(2,i),idcell) + x(1,sn(1,i),idcell))
        xm(2,i) = 0.5d0*(x(2,sn(2,i),idcell) + x(2,sn(1,i),idcell))
        df(1,i) = xm(1,i) - xc(1,idcell)
        df(2,i) = xm(2,i) - xc(2,idcell)
c ... distancia normal a aresta i      
        ca(i)   = df(1,i)*n(1,i) + df(2,i)*n(2,i)
      enddo
c ...
      return
      end
c **********************************************************************
c
c **********************************************************************
c * VORT2D : calculo da vorticidade de um campo vetorial               *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * gradU1       - gradiente da variavel u1                            *
c * gradU2       - gradiente da variavel u2                            *
c * rot          - indefinido                                          *
c * n            - dimensoes do vetores                                *
c * ndm          - numero dimensoes espaciais                          *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * rot          - rotacional                                          *
c **********************************************************************
      subroutine vort2D(gradU1,gradU2,rot,n,ndm)
      implicit none
      real*8 gradU1(ndm,*),gradU2(ndm,*),rot(*)
      integer ndm,n,i
      do i = 1, n
        rot(i) = gradU2(1,i) - gradU1(2,i)
      enddo
      return
      end
c ***********************************************************************   
c
c **********************************************************************
c * DEVITORICSTRESS : calculo do tensor desviador                      *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * gradU1       - gradiente da variavel u1                            *
c * gradU2       - gradiente da variavel u2                            *
c * dStress      - indefinido                                          *
c * n            - dimensoes do vetores                                *
c * ndm          - numero dimensoes espaciais                          *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * dStress      - tensor desviador das tensoes                        *
c **********************************************************************
      subroutine devitoricStress(gradU1,gradU2,dStress,n,ndm)
      implicit none
      real*8 gradU1(ndm,*),gradU2(ndm,*),dStress(ndm*ndm,*)
      integer ndm,n,i
      do i = 1, n
        dStress(1,i) = gradU1(1,i)                    ! du1/dx1
        dStress(2,i) = 0.5d0*(gradU1(2,i)+gradU2(1,i))!du1/dx2+du2/dx1
        dStress(3,i) = gradU2(2,i)                    ! du2/dx2 
        dStress(4,i) = dStress(2,i)                   !du2/dx1+du1/dx2
      enddo
      return
      end
c ***********************************************************************   
c
c **********************************************************************
c * STRESS : calculo do tensor                                         *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * p            - campo de pressao                                    *
c * gradU1       - gradiente da variavel u1                            *
c * gradU2       - gradiente da variavel u2                            *
c * s            - indefinido                                          *
c * n            - dimensoes do vetores                                *
c * ndm          - numero dimensoes espaciais                          *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * s            - tensor das tensoes                                  *
c **********************************************************************
      subroutine stress(p,gradU1,gradU2,s,n,ndm)
      implicit none
      real*8 gradU1(ndm,*),gradU2(ndm,*),s(ndm*ndm,*),p(*)
      integer ndm,n,i
      do i = 1, n
        s(1,i) = p(i) + gradU1(1,i)             !p + du1/dx1
        s(2,i) = 0.5d0*(gradU1(2,i)+gradU2(1,i))!du1/dx2+du2/dx1
        s(3,i) = p(i) + gradU2(2,i)             !p + du2/dx2 
        s(4,i) = s(2,i)                         !du2/dx1+du1/dx2
      enddo
      return
      end
c ***********************************************************************   
c **********************************************************************
c * WALLMODEL : Modelo de parede                                       *
c * ------------------------------------------------------------------ *
c * Parametros de entrada:                                             *
c * -------------------------------------------------------------------*
c * yPLus        - indefinido                                          *
c * stressW      - indefinido                                          *
c * viscosity    - viscosidade molecular                               *
c * specificMass - massa especifica                                    *
c * vt           - modulo da velocidade tagencial                      *
c * dy           - distancia da pareda ao centro da celula             *
c * nel          - numero da celula                                    *
c * -------------------------------------------------------------------*
c * Parmetros de saida:                                                *
c * -------------------------------------------------------------------*
c * yPLus        - distancia nao dimensional da parede                 *
c * stressW      - tensao na parede                                    *
c **********************************************************************
      subroutine wallModel(yPlus,stressW,viscosity,specificMass,vt,dy
     .                    ,nel,codWall)
      implicit none
      real*8 stressW,stressW0
      real*8 temp1,temp2,temp3,onePlusB
      real*8 yPlus3,yPlus4
      real*8 uPlusL,uPlusT
      real*8 fu0,fu,uPlus,yPlus,gama,gamai
      real*8 viscosity,vt,dy,specificMass
      real*8 vonKarmani,Ei
      real*8 , parameter :: vonKarman = 0.4187d0
      real*8 , parameter :: r16       = 0.16666666666666d0
      real*8 , parameter :: E         = 9.793d0
c ... Kader - 1981 
      real*8 , parameter :: aKader    = 0.01d0 
      real*8 , parameter :: bKader    = 5.0d0   
      real*8 , parameter :: aWerner   = 8.3d0 
      real*8 , parameter :: bWerner   = 0.142857143d0   
      integer, parameter :: maxIt     =  5000
      integer i,nel,codWall
c ... universal near wall model
      if( codWall .eq. 1 ) then
c ... wall shear stress (viscosidade)
        stressW  =  viscosity*vt/dy
c ... friction velocity
        fu     = dsqrt(stressW/specificMass)
        do i = 1, maxIt
          fu0      = fu 
c ...
          yPlus    = specificMass*fu*dy/viscosity
          uPlus    = vt/fu
c .....................................................................
c 
c ... 
          if( yPlus .le. 11.81d0) then
            fu    = 2.0d0*fu*(uPlus/(uPlus+uPlus))
          else 
            fu    = fu*((1.0d0 
     .            + 2.0d0*vonKarman*uPlus
     .            -dlog(E*yPlus))/ (1.0d0 + vonKarman*uPlus))
          endif
          if(dabs(fu-fu0) .lt. 1.0d-7) goto 100 
        enddo
        print*,'Funcao de parede nao convergiu !!!'
        print*,'Tipo : ',codWall                    
        print*,i,dabs(fu-fu0),nel
        stop
  100   continue
        yPlus = specificMass*fu*dy/viscosity
        stressW  = specificMass*(vt/uPLus)*(vt/uPLus)          
c .....................................................................
c
c ... universal near wall model kader 1981
      else if( codWall .eq. 2 ) then
        vonKarmani = 1.0d0/vonKarman
c ... wall shear stress (viscosidade)
        stressW  =  viscosity*vt/dy
c ... friction velocity(chute inicial)
        fu       = dsqrt(stressW/specificMass)
c .....................................................................
        do i = 1, maxIt
          fu0    = fu 
c ...
          yPlus    = specificMass*fu*dy/viscosity
          uPlus    = vt/fu
c ... 
          onePlusB = 1.0d0 + bKader*yPlus
          yPlus3   = yPlus*yPlus*yPlus
          yPlus4   = yPlus3*yPlus
c ... u+ laminar (linear)
          uPlusL   = yPlus
c ... u+ trubulento (log)
          uPlusT   = vonKarmani*dlog(E*yPlus)
c ... parametro T = -(a(y+)**4)/(1+by+)
          gama  = - (aKader*yPlus4) / onePlusB
          gamai = 1.0d0/gama
c ... derivada dT/dy+
          temp1 = aKader*bKader*yPlus4-4.0d0*aKader*yPlus3*onePlusB
          temp1 = temp1/(onePlusB*onePlusB)
c ... f = e(T)*(y+) + e(1/T)*(1/k)ln(Ey+) - u+
          temp2 = dexp(gama)*uPlusL + dexp(gamai)*uPlusT - uPlus
c ... derivada de f
          temp3  = (temp1*yPlus+1.0d0)*dexp(gama)*yPlus   
     .           + vonKarmani*dexp(gamai)*( 1.0d0 
     .           - (temp1/(gama*gama))*yPlus*dlog(E*yPLus))
     .           + uPlus
c .....................................................................
c
c ...
c         print*,yPlus
          fu     = fu*(1.0d0 - temp2/temp3)
c .....................................................................
          if(dabs(fu-fu0) .lt. 1.0d-7) goto 200 
        enddo
        print*,'Funcao de parede nao convergiu !!!'
        print*,'Tipo : ',codWall                    
        print*,i,dabs(temp2/temp3),nel
        stop
  200   continue
        yPlus = specificMass*fu*dy/viscosity
        stressW  =  specificMass*fu*fu
c .....................................................................
c
c ... The Werner-Wengle - 1993
      else if( codWall .eq. 3) then
        gama = 2.0d0/(1.0d0-bWerner) 
        temp1 = 0.5d0*(viscosity/(specificMass*dy))
        temp1 = temp1*(aWerner**gama)
        if( vt .le. temp1) then
          stressW = (2.0d0*viscosity*vt)/dy
        else
          temp3    = viscosity/(specificMass*dy)
          onePlusB = 1.0d0 + bwerner
c
          gama     = onePlusB /(1.0d0 - aWerner)
          temp1    = 0.5d0*(1.0d0-bWerner)*(aWerner**gama)
          temp1    = temp1*temp3**(onePlusB) 
c
          temp2    = (onePlusB/aWerner)*(temp3**bwerner)*vt
c 
          temp2    = (temp1+temp2)**(2.0d0/onePlusB)
c
          stressW  =  specificMass*temp2       
        endif 
c ... friction velocity
        fu     = dsqrt(stressW/specificMass)
c ...
        yPlus  = specificMass*fu*dy/viscosity 
c .....................................................................
c
c ... Spaldings law - Villiers - 2006
      else if( codWall .eq. 4) then
        Ei = 1.0d0/E
        stressW  =  viscosity*vt/dy
        fu       = dsqrt(stressW/specificMass)
        do i = 1, maxIt
          fu0    = fu 
c                              
          yPlus = specificMass*fu*dy/viscosity
          uPlus = vt/fu
c ... friction velocity
          gama   = vonKarman*uPlus
c
          temp1  = uPlus - yPlus 
     .           + (Ei)*(dexp(gama) - 1.0d0  
     .           - 0.5d0*gama*gama - r16*gama*gama*gama)
c
          temp2  = -uPlus - yPlus 
     .           + (Ei)*(-gama*dexp(gama) + gama  
     .           + gama*gama + 0.5d0*gama*gama*gama)
          fu     = fu*(1.0d0 - temp1/temp2)
          if(dabs(fu-fu0) .lt. 1.0d-7) goto 300 
        enddo
        print*,'Funcao de parede nao convergiu !!!'
        print*,'Tipo : ',codWall                    
        print*,i,dabs(fu-fu0),nel
        stop
  300   continue
        yPlus = specificMass*fu*dy/viscosity
        stressW  =  specificMass*fu*fu
      endif
c ......................................................................
c
c ...
      return
      end
c ***********************************************************************   
