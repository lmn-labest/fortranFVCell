c *********************************************************************
c * PCG: Solucao de sistemas de equacoes pelo metodo dos gradientes   *
c * conjugados com precondicionador diagonal para matrizes simetricas.*  
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * neq    - numero de equacoes                                       *
c * nad    - numero de elementos nao nulos fora da diagonal principal *
c * ia     - ponteiro do arranjo do csr da matrix A                   *
c * ja     - arranjo csr da matriz A                                  *
c * ad     - diagonal principal de A                                  *
c * au(nad)- coeficiantes fora da diagonal principal de A             *
c * al(nad)- coeficiantes fora da diagonal principal de A             *
c * m(neq) - precondicianado diagonal                                 *
c * b(neq) - termo indepedente                                        *
c * x(neq) - nao definido                                             *
c * z(neq) - auxiliar                                                 *
c * r(neq) - auxiliar                                                 *
c * tol    - tolerancia do solver                                     *
c * maxit  - numero maximo de iteracao do solver                      *
c * nlit     - numero da iteraçao nao lienar                          *
c * init     - inicializacao do vetor do resultados com nulo          *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * x    - solucao do sistema                                         *
c * ----------------------------------------------------------------- *
c *********************************************************************      
      subroutine pcg(neq,nad,ia,ja,ad,au,al,m,b,x,z,r,tol,maxit,matvec
     .              ,dot,nlit,init)
      implicit none
      include 'time.fi'
      real*8 au(*),al(*),ad(*),x(*),b(*)
      real*8 r(*),z(*),m(*)
      real*8 d,conv,tol,alpha,beta,energy
      real*8 dot
      integer neq,ia(*),ja(*)
      integer i,j,maxit,nad,nlit
      logical init
      external matvec,dot
c ...  
      time0 = get_time()
      if(init) then      
        do i = 1, neq
          x(i) = 0.d0
        enddo
      endif
c ..................................................................... 
c
c ...
      call matvec(ad,au,al,x,z,ia,ja,neq)
      do i = 1, neq
        r(i) = b(i) - z(i)
        z(i) = r(i) / m(i)
        b(i) = z(i)
      enddo
      d    = dot(r,z,neq)
      conv = tol*dsqrt(dabs(d))
c ..................................................................... 
      do j = 1, maxit
         call matvec(ad,au,al,b,z,ia,ja,neq)
         alpha = d/dot(b,z,neq)
         do i = 1, neq
            x(i) = x(i) + alpha * b(i)
            r(i) = r(i) - alpha * z(i)
            z(i) = r(i) / m(i)
         enddo
         beta = dot(r,z,neq) / d
         do i = 1, neq
            b(i) = z(i) + beta * b(i)
         enddo   
         d = beta*d
         if (dsqrt(dabs(d)) .lt. conv) goto 300
      enddo
c ----------------------------------------------------------------------
      write(*,1200) maxit,dsqrt(dabs(d)),conv
      stop
  300 continue
c
c ... Energy norm:
c
      call matvec(ad,au,al,x,z,ia,ja,neq)
      energy   = dot(x,z,neq)
c ......................................................................
      time = get_time()
      time = time-time0
c ----------------------------------------------------------------------
      if(nlit .eq. 1) write(*,1100)tol,neq,nad,j,energy,time
c ......................................................................
c     Controle de flops
      write(13,'(a,a,i9,a,d20.10,a,d20.10)')
     .      "PCG:       ","it",j, " energy norm ",energy,
     .      " tol ",tol
c ......................................................................
      return
c ======================================================================
 1100 format(' (PCG) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING:  PCG no convergence reached after ',i9,
     .        ' iterations !',d20.6,' Res ',d20.6,' Conv'/)
      end
c *********************************************************************
c    
c *********************************************************************
c * PBICGSTAB: Solucao de sistemas de equacoes pelo metodo            *
c * dos gradientes biconjugados com precondicionador diagonal para    *
c * matrizes nao-simetricas.                                          * 
c * ------------------------------------------------------------------*
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * neq    - numero de equacoes                                       *
c * nad    - numero de elementos nao nulos fora da diagonal principal *
c * ia     - ponteiro do arranjo do csr da matrix A                   *
c * ja     - arranjo csr da matriz A                                  *
c * ad(neq)- diagonal principal de A                                  *
c * au(nad)- coeficiantes fora da diagonal principal de A             *
c * al(nad)- coeficiantes fora da diagonal principal de A             *
c * m(neq) - precondicianado diagonal                                 *
c * b(neq) - termo indepedente                                        *
c * x(neq) - nao definido                                             *
c * t(neq) - nao definido                                             *
c * r(neq) - nao definido                                             *
c * p(neq) - nao definido                                             *
c * z(neq) - nao definido                                             *
c * tol    - tolerancia do solver                                     *
c * maxit  - numero maximo de iteracao do solver                      *
c * nlit     - numero da iteraçao nao lienar                          *
c * init     - inicializacao do vetor do resultados com nulo          *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * x    - solucao do sistema                                         *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine pbicgstab(neq,nad,ia,ja,ad,au,al,m,b,x,t,v,r,p,z,tol
     .                    ,maxit,matvec,dot,nlit,init) 
      implicit none
      include 'time.fi'
c .....................................................................   
      integer neq,maxit,nad,i,j,nlit
      integer ja(*),ia(*)
      real*8  ad(*),au(*),al(*),m(*),x(*),r(*),p(*),b(*),t(*),v(*),z(*)
      real*8  dot,tol,conv,energy,d,alpha,beta,rr0,w
      logical init
      external dot,matvec
c ======================================================================
c ......................................................................
      time0 = get_time()
c ......................................................................
c
c ... Chute inicial:
c
      if(init) then
        do i = 1, neq
           x(i) = 0.d0
        enddo
      endif
c ----------------------------------------------------------------------
      call matvec(ad,au,al,x,z,ia,ja,neq)
      do i = 1, neq
         r(i) = b(i) - z(i)
         p(i) = r(i)
         b(i) = p(i)
         z(i) = p(i)/m(i) 
      enddo   
      d    = dot(r,z,neq)
      conv = tol*dsqrt(dabs(d))
c ----------------------------------------------------------------------
      do j = 1, maxit
         call matvec(ad,au,al,z,v,ia,ja,neq)
         rr0 = dot(b,r,neq)
         alpha = rr0/dot(v,r,neq)
         do i = 1, neq
            x(i) = x(i) + alpha * z(i)
            b(i) = b(i) - alpha * v(i)
            z(i) = b(i) / m(i)
         enddo    
         call matvec(ad,au,al,z,t,ia,ja,neq)
         w = dot(t,b,neq) / dot(t,t,neq)
         do i = 1, neq
            x(i) = x(i) + w*z(i)
            b(i) = b(i) - w*t(i)
         enddo   
         beta = dot(b,z,neq)/d
         if (dsqrt(dabs(beta)) .lt. conv) goto 300
         beta = (dot(r,b,neq) / rr0)*(alpha/w)
         do i = 1, neq
             p(i) = b(i) + beta*(p(i)-w*v(i))
             z(i) = p(i)/m(i)
         enddo
      enddo   
c ----------------------------------------------------------------------
      write(*,1200) maxit,d,conv
      stop
  300 continue
c
c ... Energy norm:
c
      call matvec(ad,au,al,x,z,ia,ja,neq)
      energy   = dot(x,z,neq)
c ......................................................................
c     time = MPI_Wtime()
      time = get_time()-time0
c ----------------------------------------------------------------------
      if(nlit.eq.1) then
        write(*,1100)tol,neq,nad,j,energy,time
      endif
c ......................................................................
c
c     Controle de flops
      write(13,'(a,a,i9,a,d20.10,a,d20.10)')
     .         "PBICGSTAB: ","it",j, " energy norm ",energy,
     .         " tol ",tol
c ......................................................................
      return
c ======================================================================
 1100 format(' (PBICGSTAB) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1150 format(5x,'Number of iterations = ',i20/
     .       5x,'Energy norm          = ',d20.10/
     .       5x,'CPU time (s)         = ',f20.2/)
 1200 format(' *** WARNING: BiPCGSTAB no convergence reached after ',i9,
     .       ' iterations !',d20.6,' Res ',d20.6,' Conv'/)
      end
c **********************************************************************
c *   GMRES: Solucao iterativa de sistemas simetricos e nao-simetricos *
c *          pelo metodo GMRES com precondicionador diagonal.          *
c * ------------------------------------------------------------------ *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c *   neq    - numero de equacoes                                      *
c *   ia(*)  - ponteiro do formato CSR                                 *
c *   ja(*)  - ponteiro das colunas no formato CSR                     *
c *   ad(neq)- diagonal da matriz A                                    *
c *   au(*)  - parte triangular superior de A                          *
c *   al(*)  - parte triangular inferior de A                          *
c *   m(*)   - precondicionador diagonal                               *
c *   b(neq) - vetor de forcas                                         *
c *   x(neq) - chute inicial                                           *
c *   k      - numero de bases                                         *
c *   z(neq) - arranjo local de trabalho                               *
c *   r(neq) - arranjo local de trabalho                               *
c *   tol    - tolerancia de convergencia                              *
c *   maxit  - numero maximo de iteracoes                              *
c *   matvec - nome da funcao externa para o produto matrix-vetor      *
c *   dot    - nome da funcao externa para o produto escalar           *
c *                                                                    *
c *   Arranjos locais de trabalho:                                     *
c *                                                                    *
c *      g(neq+1,k+1)                                                  *
c *      h(k+1,k)                                                      *
c *      y(k)                                                          *
c *      c(k)                                                          *
c *      s(k)                                                          *
c *      e(k+1)                                                        *
c * ------------------------------------------------------------------ *
c * Parametros de saida :                                              *
c * ------------------------------------------------------------------ *
c * x    - solucao do sistema                                          *
c **********************************************************************                 
      subroutine gmres(neq,nad,ia,ja,ad,au,al,m,b,x,k,g,h,y,c,s,e,
     .                 tol,maxit,matvec,nlit)
      implicit none
      include 'time.fi'
c .....................................................................      
      integer neq,k,maxit,ia(*),ja(*),nit,i,j,l,ni,ic,nad
      real*8  ad(*),au(*),al(*),m(*),b(*),x(*)
      real*8  g(neq,k+1),h(k+1,k),y(k),c(k),s(k),e(k+1),tol
      real*8  energy,econv,norm,dot,ddot,r,aux1,aux2,beta
      integer nii(maxit),nlit
      external matvec,dot
c ......................................................................
      time0 = get_time()
c ......................................................................
c
c.... Chute inicial:
c
      do 10 i = 1, neq
         x(i) = 0.d0
   10 continue
c ----------------------------------------------------------------------
c
c ... Limite de convergencia:
c
      norm  = dsqrt(dot(b,b,neq))
      econv = tol*norm
c ----------------------------------------------------------------------
c
c ... Ciclos GMRES:
c
      nit = 0
      do 1000 l = 1, maxit
c
c ...... Residuo g(1) = b - A x:
c 
         call matvec(ad,au,al,x,g,ia,ja,neq)
c
c ...... Residuo com precondicionador diagonal:
c
         do 200 i = 1, neq
            g(i,1) = (b(i) - g(i,1))/m(i)
  200    continue
c
c ...... Norma do residuo:
c
         e(1) = dsqrt(dot(g,g,neq))
c
c ...... Normalizacao de g1:
c
         do 210 i = 1, neq
            g(i,1) = g(i,1)/e(1)
  210    continue
c
c ...... Iteracoes GMRES:
c
         ni = 0
         do 400 i = 1, k
            nit = nit + 1
            ni  = ni  + 1
c
c ......... Produto g(i+1) = A.g(i):
c
            call matvec(ad,au,al,g(1,i),g(1,i+1),ia,ja,neq)
c
c ......... Precondicionador diagonal:
c
            do 300 j = 1, neq
               g(j,i+1) = g(j,i+1)/m(j)
  300       continue
c
c ......... Ortogonalizacao (Gram-Schmidt modificado):
c
            do 320 j = 1, i
               beta = dot(g(1,i+1),g(1,j),neq)
               do 310 ic = 1, neq
                  g(ic,i+1) = g(ic,i+1) - beta * g(ic,j)
  310          continue
               h(j,i) = beta
  320       continue
c
c ......... Norma de g(i+1):
c
            norm = dsqrt(dot(g(1,i+1),g(1,i+1),neq))
c
            h(i+1,i) = norm
c
c ......... Normalizacao de g(i+1):
c
            do 330 ic = 1, neq
               g(ic,i+1) = g(ic,i+1)/norm
  330       continue
c
            do 340 j = 1, i-1
               aux1 =  c(j) * h(j,i) + s(j) * h(j+1,i)
               aux2 = -s(j) * h(j,i) + c(j) * h(j+1,i)
               h(j,i)   = aux1
               h(j+1,i) = aux2
  340       continue
            r = dsqrt(h(i,i)*h(i,i) + h(i+1,i)*h(i+1,i))
            c(i) = h(i,i)/r
            s(i) = h(i+1,i)/r
            h(i,i)   = r
            h(i+1,i) = 0.d0
            e(i+1) = -s(i) * e(i)
            e(i)   =  c(i) * e(i)
            if (dabs(e(i+1)) .le. econv) goto 500
  400    continue
  500    continue
c
c ...... Resolve o sistema h y = e :
c
         y(ni) = e(ni) / h(ni,ni)
         do 520 i = ni-1, 1, -1
            y(i) = 0.d0
            do 510 j = i+1, ni
               y(i) = y(i) - h(i,j)*y(j)
  510       continue
            y(i) = (y(i) + e(i)) / h(i,i)
  520    continue
c
c ...... Atualizacao de x:
c
         do 610 i = 1, neq
            do 600 j = 1, ni
               x(i) = x(i) + y(j) * g(i,j)
  600       continue
  610    continue
c
c ...... Verifica a convergencia:
c
         nii(l)=ni
         if (dabs(e(ni+1)) .le. econv) goto 1100
c         print*,dsqrt(dabs(e(ni+1)))
 1000 continue
      write(*,2100) maxit
      stop 
c ----------------------------------------------------------------------
 1100 continue
c
c ... Norma de energia da solucao
c
      energy = dot(x,b,neq)
c ......................................................................
      time = get_time()-time0
c ----------------------------------------------------------------------
c      write(*,2000) tol,neq,l,nit,dabs(e(ni+1)),energy,time
c     if (dabs(e(ni+1)) .gt. econv) then
c        write(*,2100) maxit
c        stop
c     endif
      if(nlit.eq.1) then
        write(*,2000)tol,neq,nad
        write(*,2050)l,nit,abs(e(ni+1)),energy,time
      else
        write(*,2050)l,nit,abs(e(ni+1)),energy,time
      endif
c ......................................................................
c     Controle de flops
      write(13,'(a,a,i9,a,d20.10,a,d20.10)')
     .         "GMRES: ","it",l, " energy norm ",energy,
     .         " tol ",tol
c ......................................................................
      return
c ----------------------------------------------------------------------
 2000 format(' (GMRES) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Nad                  = ',i20)
 2050 format(5x,'Number of cycles     = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Norm                 = ',d20.6/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 2100 format(' *** WARNING: no convergence reached for ',i9,' cycles !',
     . /)
      end
