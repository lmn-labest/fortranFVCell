c *********************************************************************
c * PCG_OMP:Solucao de sistemas de equacoes pelo metodo dos gradientes*
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
c * thread_y - buffer                                                 *
c * nlit     - numero da iteraçao nao lienar                          *
c * init     - inicializacao do vetor do resultados com nulo          *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * x    - solucao do sistema                                         *
c * ----------------------------------------------------------------- *
c *********************************************************************      
      subroutine pcg_omp(neq,nad,ia,ja,ad,au,al,m,b,x,z,r,tol,maxit
     .                 ,matvec,dot,thread_y,nlit,init)
      implicit none
      include 'openmp.fi'
      include 'time.fi'
      real*8 au(*),al(*),ad(*),x(*),b(*)
      real*8 r(*),z(*),m(*),thread_y(*)
      real*8 d,conv,tol,alpha,beta,energy
      real*8 dot
      integer neq,ia(*),ja(*)
      integer i,j,maxit,nad,nlit
      logical init
      external matvec,dot
c ...  
      time0 = get_time()
c$omp parallel private(i,j,d,conv,beta,alpha)  
c$omp.shared(neq,nad,ia,ja,al,ad,au,b,x,m,z,r,tol,maxit,thread_y,nlit)
c$omp.num_threads(nThreadsSolver)
c
c ...
      if(init) then
c$omp do      
        do i = 1, neq
          x(i) = 0.d0
        enddo
c$omp end do
      endif
c ..................................................................... 
c
c ...
      call matvec(ad,au,al,x,z,ia,ja,neq,thread_y)
c$omp do
      do i = 1, neq
        r(i) = b(i) - z(i)
        z(i) = r(i) / m(i)
        b(i) = z(i)
      enddo
c$omp end do
      d    = dot(r,z,neq)
      conv = tol*dsqrt(dabs(d))
c ..................................................................... 
      do j = 1, maxit
         call matvec(ad,au,al,b,z,ia,ja,neq,thread_y)
         alpha = d/dot(b,z,neq)
c$omp do
         do i = 1, neq
            x(i) = x(i) + alpha * b(i)
            r(i) = r(i) - alpha * z(i)
            z(i) = r(i) / m(i)
         enddo
c$omp end do
         beta = dot(r,z,neq) / d
c$omp do
         do i = 1, neq
            b(i) = z(i) + beta * b(i)
         enddo
c$omp end do            
         d = beta*d
         if (dsqrt(dabs(d)) .lt. conv) goto 300
       enddo
c ----------------------------------------------------------------------
c$omp single      
      write(*,1200) maxit,dsqrt(dabs(d)),conv
      stop
c$omp end single  
  300 continue
c
c ... Energy norm:
c
      call matvec(ad,au,al,x,z,ia,ja,neq,thread_y)
      energy   = dot(x,z,neq)
c ......................................................................
c$omp single
      time = get_time()
      time = time-time0
c ----------------------------------------------------------------------
      if(nlit .eq. 1) write(*,1100)tol,neq,nad,j,energy,time
c ......................................................................
c     Controle de flops
      write(13,'(a,a,i9,a,d20.10,a,d20.10)')
     .      "PCG OMP      : ","it",j, " energy norm ",energy,
     .      " tol ",tol
c ......................................................................
c$omp end single
c$omp end parallel 
c ......................................................................
      return
c ======================================================================
 1100 format(' (PCG_OMP) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING:  PCG OMP no convergence reached after ',i9,
     .        ' iterations !',d20.6,' Res ',d20.6,' Conv'/)
      end
c *********************************************************************
c
c *********************************************************************
c * PBICGSTAB_OMP: Solucao de sistemas de equacoes pelo metodo        *
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
      subroutine pbicgstab_omp(neq,nad,ia,ja,ad,au,al,m,b,x,t,v,r,p,z
     .                        ,tol,maxit,matvec,dot,nlit,init) 
      implicit none
      include 'openmp.fi'
      include 'time.fi'
c .....................................................................   
      integer neq,maxit,nad,i,j,nlit
      integer ja(*),ia(*)
      real*8  ad(*),au(*),al(*),m(*),x(*),r(*),p(*),b(*),t(*),v(*),z(*)
      real*8  dot,tol,conv,energy,d,alpha,beta,rr0,w,dum
      logical init
      external dot,matvec
c ======================================================================
c ......................................................................
      time0 = get_time()
c ......................................................................
c$omp parallel private(i,j,d,conv,beta,alpha,rr0,w,energy)
c$omp.shared(neq,nad,ia,ja,al,ad,au,b,x,m,t,v,r,p,z,tol,maxit,nlit)
c$omp.num_threads(nThreadsSolver)
c
c ... Chute inicial:
c
      if(init) then
c$omp do      
        do i = 1, neq
          x(i) = 0.d0
        enddo
c$omp end do
      endif
c ----------------------------------------------------------------------
      call matvec(ad,au,al,x,z,ia,ja,neq,dum)
c$omp do
      do i = 1, neq
         r(i) = b(i) - z(i)
         p(i) = r(i)
         b(i) = p(i)
         z(i) = p(i)/m(i) 
      enddo
c$omp enddo         
      d    = dot(r,z,neq)
      conv = tol*dsqrt(dabs(d))
c ----------------------------------------------------------------------
      do j = 1, maxit
         call matvec(ad,au,al,z,v,ia,ja,neq,dum)
         rr0 = dot(b,r,neq)
         alpha = rr0/dot(v,r,neq)
c$omp do
         do i = 1, neq
            x(i) = x(i) + alpha * z(i)
            b(i) = b(i) - alpha * v(i)
            z(i) = b(i) / m(i)
         enddo
c$omp enddo
         call matvec(ad,au,al,z,t,ia,ja,neq,dum)
         w = dot(t,b,neq) / dot(t,t,neq)
c$omp do
         do i = 1, neq
            x(i) = x(i) + w*z(i)
            b(i) = b(i) - w*t(i)
         enddo
c$omp enddo            
         beta = dot(b,z,neq)/d
         if (dsqrt(dabs(beta)) .lt. conv) goto 300
         beta = (dot(r,b,neq) / rr0)*(alpha/w)
c$omp do
         do i = 1, neq
             p(i) = b(i) + beta*(p(i)-w*v(i))
             z(i) = p(i)/m(i)
         enddo
c$omp enddo
      enddo   
c ----------------------------------------------------------------------
c$omp single
      write(*,1200) maxit,beta,conv
      stop
c$omp end single
  300 continue
c
c ... Energy norm:
c
      call matvec(ad,au,al,x,z,ia,ja,neq,dum)
      energy   = dot(x,z,neq)
c ......................................................................
c$omp single
      time = get_time()-time0
c ----------------------------------------------------------------------
      if(nlit.eq.1) then
        write(*,1100)tol,neq,nad,j,energy,time
      endif
c ......................................................................
c
c     Controle de flops
      write(13,'(a,a,i9,a,d20.10,a,d20.10)')
     .         "PBICGSTAB OMP: ","it",j, " energy norm ",energy,
     .         " tol ",tol
c ......................................................................
c$omp end single
c$omp end parallel
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
c1150 format(5x,'Number of iterations = ',i20/
c    .       5x,'Energy norm          = ',d20.10/
c    .       5x,'CPU time (s)         = ',f20.2/)
 1200 format(' *** WARNING: BiPCGSTAB OMP no convergence reached after '
     .      ,i9,' iterations !',d20.6,' Res ',d20.6,' Conv'/)
      end    
