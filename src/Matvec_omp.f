c *********************************************************************
c * MATVEC_CSRD_OMP: produto matriz vetor no formato csr ( grafo de A *
c * nao simetrico; coeficiente de A nao simetrico; ad - diagonal,     *
c * a - superior e inferior)                                          *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ad     - diagonal princiapal de A                                 *
c * a      - vaolores nao nulos fora da diagonal principal            *
c * x      - valores do vetor a ser multiplicado por A                *
c * y      - nao definido                                             *
c * ia     - ponteiro do arranjo do csr da matrix A                   *
c * ja     - arranjo csr da matriz A                                  *
c * neq    - numero de equacoes                                       *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * y    - vetor com a operacao Ax                                    *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine matvec_csrd_omp(ad,a,dum,x,y,ia,ja,neq)
      implicit none
      include 'time.fi'
      real*8 ad(*),a(*),x(*),y(*),dum,t
      integer ja(*),neq,i,j,ia(*)
c ......................................................................
c
c ...
c$omp single      
      matvectime = get_time() - matvectime
c$omp end single
c ......................................................................
c
c ...
c$omp do private(i,j,t)      
      do i = 1, neq
c ... diagonal principal
        t = ad(i)*x(i)
c ... produto da linha
        do j = ia(i) , ia(i+1) - 1
          t  = t + a(j)*x(ja(j))
        enddo
        y(i) = t
      enddo
c$omp end do
c .....................................................................
c
c ...
c$omp single      
      matvectime = get_time() - matvectime
c$omp end single
      return
      end
c *********************************************************************
c
c *********************************************************************
c * MATVEC_CSRD_ILU2_OMP: produto matriz vetor no formato csr         * 
c *  ( grafo de A nao simetrico; coeficiente de A nao simetrico       *
c *   ad - diagonal,  a - superior e inferior)                        *
c *                   (loop interno desenrolado)                      *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ad     - diagonal princiapal de A                                 *
c * a      - vaolores nao nulos fora da diagonal principal            *
c * x      - valores do vetor a ser multiplicado por A                *
c * y      - nao definido                                             *
c * ia     - ponteiro do arranjo do csr da matrix A                   *
c * ja     - arranjo csr da matriz A                                  *
c * neq    - numero de equacoes                                       *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * y    - vetor com a operacao Ax                                    *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine matvec_csrd_ilu2_omp(ad,a,dum,x,y,ia,ja,neq)
      implicit none
      include 'time.fi'
      real*8 ad(*),a(*),x(*),y(*),dum,t
      integer ja(*),neq,i,j,ia(*),k1,k2,n
c ......................................................................
c
c ...
c$omp single      
      matvectime = get_time() - matvectime
c$omp end single
c ......................................................................
c
c ...
c$omp do private(i,j,t,k1,k2)      
      do i = 1, neq
c ... diagonal principal
        t = ad(i)*x(i)
        k1 = ia(i)
        k2 = ia(i+1)
        if( k2 .eq. k1) goto 120
        n = mod(k2-k1,2)
        if(n .eq. 0) goto 100
        t  = t + a(k1)*x(ja(k1))
        k1 = k1 + 1
c ... produto da linha
  100   do j = k1 , k2 - 1, 2
          t  = t + a(j)*x(ja(j)) + a(j+1)*x(ja(j+1))
        enddo
  120   y(i) = t
      enddo
c$omp end do
c .....................................................................
c
c ...
c$omp single      
      matvectime = get_time() - matvectime
c$omp end single
      return
      end
c **********************************************************************
c
c *********************************************************************
c * MATVEC_CSRD_ILU4_OMP: produto matriz vetor no formato csr         * 
c *  ( grafo de A nao simetrico; coeficiente de A nao simetrico       *
c *   ad - diagonal,  a - superior e inferior)                        *
c *                   (loop interno desenrolado)                      *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ad     - diagonal princiapal de A                                 *
c * a      - vaolores nao nulos fora da diagonal principal            *
c * x      - valores do vetor a ser multiplicado por A                *
c * y      - nao definido                                             *
c * ia     - ponteiro do arranjo do csr da matrix A                   *
c * ja     - arranjo csr da matriz A                                  *
c * neq    - numero de equacoes                                       *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * y    - vetor com a operacao Ax                                    *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine matvec_csrd_ilu4_omp(ad,a,dum,x,y,ia,ja,neq)
      implicit none
      include 'time.fi'
      real*8 ad(*),a(*),x(*),y(*),dum,t
      integer ja(*),neq,i,j,ia(*),k1,k2,n
c ......................................................................
c
c ...
c$omp single      
      matvectime = get_time() - matvectime
c$omp end single
c ......................................................................
c
c ...
c$omp do private(i,j,t,k1,k2)      
      do i = 1, neq
c ... diagonal principal
        t = ad(i)*x(i)
        k1 = ia(i)
        k2 = ia(i+1)
        if( k2 .eq. k1) goto 120
        n = mod(k2-k1,4)
        if(n .eq. 0) goto 100
        t  = t + a(k1)*x(ja(k1))
        k1 = k1 + 1
        if(k2 .eq. k1) goto 120
        n = n - 1
        if(n .eq. 0) goto 100
        if(n .eq. 1) then
          t  = t + a(k1)*x(ja(k1))
          k1 = k1 + 1
          if(k2 .eq. k1) goto 120
        else
          t = t + a(k1)*x(ja(k1)) + a(k1+1)*x(ja(k1+1))
          k1 = k1 + 2
          if( k2 .eq. k1) goto 120
        endif
c ... produto da linha
  100   do j = k1 , k2 - 1, 4
          t  = t + a(j)*x(ja(j))     + a(j+1)*x(ja(j+1)) 
     .           + a(j+2)*x(ja(j+2)) + a(j+3)*x(ja(j+3))
        enddo
  120   y(i) = t
      enddo
c$omp end do
c .....................................................................
c
c ...
c$omp single      
      matvectime = get_time() - matvectime
c$omp end single
      return
      end
c **********************************************************************
c
c *********************************************************************
c * MATVEC_CSRD_iLO2_ILU2_OMP: produto matriz vetor no formato csr    * 
c *  ( grafo de A nao simetrico; coeficiente de A nao simetrico       *
c *   ad - diagonal,  a - superior e inferior)                        *
c *                   (loop externo/interno desenrolado)              *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ad     - diagonal princiapal de A                                 *
c * a      - vaolores nao nulos fora da diagonal principal            *
c * x      - valores do vetor a ser multiplicado por A                *
c * y      - nao definido                                             *
c * ia     - ponteiro do arranjo do csr da matrix A                   *
c * ja     - arranjo csr da matriz A                                  *
c * neq    - numero de equacoes                                       *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * y    - vetor com a operacao Ax                                    *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine matvec_csrd_ilo2_ilu2_omp(ad,a,dum,x,y,ia,ja,neq)
      implicit none
      include 'time.fi'
      real*8 ad(*),a(*),x(*),y(*),dum,t,t1
      integer ja(*),neq,i,j,ia(*),k1,k2,k3,n,nn,i0
c ......................................................................
c
c ...
c$omp single      
      matvectime = get_time() - matvectime
c$omp end single
c ......................................................................
c
c ...
      i0 = 1
      n  = mod(neq,2)
      if(n.eq.0) goto 10
      y(1) = ad(1)*x(1)
      i0 = 2
   10 continue
c$omp do private(i,j,t,t1,k1,k2,k3,nn)  
      do i = i0, neq,2
c ... diagonal principal
        t  =   ad(i)*x(i)
        t1 = ad(i+1)*x(i+1)
        k1 = ia(i)
        k2 = ia(i+1)
        k3 = ia(i+2)
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
        if( k2 .eq. k1) goto 120
        nn = mod(k2-k1,2)
        if(nn .eq. 0) goto 100
        t  = t + a(k1)*x(ja(k1))
        k1 = k1 + 1
c ... produto da linha
  100   do j = k1 , k2 - 1, 2
          t  = t + a(j)*x(ja(j)) + a(j+1)*x(ja(j+1))
        enddo
c
c ...    Armazena o resultado em y(i):
c 
  120   y(i) = t
c -------------------------------------------------------------
c
c ...    Loop nos coeficientes nao nulos da linha i+1:
c      
        if( k3 .eq. k2) goto 220
        nn = mod(k3-k2,2)
        if(nn .eq. 0) goto 200
        t1  = t1 + a(k2)*x(ja(k2))
        k2  = k2 + 1
c ... produto da linha
  200   do j = k2 , k3 - 1, 2
          t1 = t1 + a(j)*x(ja(j)) + a(j+1)*x(ja(j+1))
        enddo     
c
c ...    Armazena o resultado em y(i):
c
  220   y(i+1) = t1      
      enddo
c$omp end do
c .....................................................................
c
c ...
c$omp single      
      matvectime = get_time() - matvectime
c$omp end single
      return
      end
c *********************************************************************
c
c *********************************************************************
c * MATVEC_CSRD_SYM_OMP: produto matriz-vetor y = Ax  (A simetrica),  *
c *                   coef. de A no formato CSR e armazenamento       *
c *                   da parte superior, com loops aninhados.         *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ad       - diagonal princiapal de A                               *
c * al(nad)  - parte triangular inferior de A no formato CSR          *
c * x        - valores do vetor a ser multiplicado por A              *
c * y        - nao definido                                           *
c * ia       - ponteiro do arranjo do csr da matrix A                 *
c * ja       - arranjo csr da matriz A(inferior)                      *
c * neq      - numero de equacoes                                     *
c * thread_y - buffer                                                 *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * y    - vetor com a operacao Ax                                    *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine matvec_csrd_sym_omp(ad,dum,al,x,y,ia,ja,neq,thread_y)
      implicit none
      include 'openmp.fi'
      include 'time.fi'
      real*8 ad(*),al(*),x(*),y(*),dum,t,thread_y(*),xi,s
      integer ja(*),neq,i,j,k,ia(*),jak,inc
c ......................................................................
c$omp single      
      matvectime = get_time() - matvectime
c$omp end single
c
c ... inciando buffer      
c$    thread_id = omp_get_thread_num() + 1
      do i = 1, nThreadsSolver
         inc = (i-1)*neq
c$omp do
         do j = thread_height(i)+inc, thread_begin(i)+inc-1
            thread_y(j) = 0.d0
         enddo
c$omp end do
      enddo
c .....................................................................
c
c ... matVec
      inc = (thread_id - 1)*neq
c$omp barrier
      do i = thread_begin(thread_id), thread_end(thread_id)
         y(i) = 0.d0
         xi   = x(i)
c
c ...    Produto da diagonal de A por x:
c
         t  = ad(i)*xi
c
         do k = ia(i), ia(i+1)-1
            jak = ja(k)
            s   = al(k)
c
c ...    Parte superior de A, linha i:
c
            t = t + s*x(jak)
c
c ...    Parte inferior de A, coluna i:
c
            jak = jak + inc   
            thread_y(jak) = thread_y(jak) + s*xi
         enddo   
c
c ...    Armazena o resultado em y(i):
c
         thread_y(i+inc) =  t
      enddo
c$omp barrier
c
c ... Resgatando e acumulando os valores thread_y -> y
      do i = 1,nThreadsSolver
         inc = (i-1)*neq
c$omp do
         do j = thread_height(i), thread_end(i)
            y(j) = y(j) + thread_y(j+inc)
         end do
c$omp end do
      end do
c .....................................................................
c
c ... 
c$omp single          
      matvectime = get_time() - matvectime
c$omp end single
      return
      end
c **********************************************************************
c
c *********************************************************************
c * MATVEC_CSRD_SYM_ILU2_OMP: produto matriz-vetor                    *
c *                   y = Ax (A simetrica)                            *
c *                   coef. de A no formato CSR e armazenamento       *
c *                   da parte inferio, com loops aninhados.          *
c *                   (loop interno desenrolado)                      *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ad       - diagonal princiapal de A                               *
c * al(nad)  - parte triangular inferior de A no formato CSR          *
c * x        - valores do vetor a ser multiplicado por A              *
c * y        - nao definido                                           *
c * ia       - ponteiro do arranjo do csr da matrix A                 *
c * ja       - arranjo csr da matriz A(inferior)                      *
c * neq      - numero de equacoes                                     *
c * thread_y - buffer                                                 *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * y    - vetor com a operacao Ax                                    *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine matvec_csrd_sym_ilu2_omp(ad,dum,al,x,y,ia,ja,neq
     .                                   ,thread_y)
      implicit none
      include 'openmp.fi'
      include 'time.fi'
      real*8 ad(*),al(*),x(*),y(*),dum,t,thread_y(*),xi,s,s1
      integer ja(*),neq,i,j,k,ia(*),jak,inc,k1,k2,n,jak1
c ......................................................................
c$omp single      
      matvectime = get_time() - matvectime
c$omp end single
c
c ... inciando buffer      
c$    thread_id = omp_get_thread_num() + 1
      do i = 1, nThreadsSolver
         inc = (i-1)*neq
c$omp do
         do j = thread_height(i)+inc, thread_begin(i)+inc-1
            thread_y(j) = 0.d0
         enddo
c$omp end do
      enddo
c .....................................................................
c
c ... matVec
      inc = (thread_id - 1)*neq
c$omp barrier
      do i = thread_begin(thread_id), thread_end(thread_id)
         y(i) = 0.d0
         xi   = x(i)
c
c ...    Produto da diagonal de A por x:
c
         t  = ad(i)*xi
c
         k1 = ia(i)
         k2 = ia(i+1)
         if(k2 .eq. k1) goto 102
         n = mod(k2-k1,2)
         if(n .eq. 0) goto 101
         jak1 = ja(k1)
         t    =    t + al(k1)*x(jak1)
         jak1 = jak1 + inc
         thread_y(jak1) = thread_y(jak1) + al(k1)*xi
         k1 = k1 + 1
c
  101    do k = k1, k2-1,2
            jak  = ja(k)
            jak1 = ja(k+1)
            s    = al(k)
            s1   = al(k+1)
c
c ...    Parte superior de A, linha i:
c
            t = t + s*x(jak) + s1*x(jak1)
c
c ...    Parte inferior de A, coluna i:
c
            jak  = jak  + inc
            jak1 = jak1 + inc
            thread_y(jak)  = thread_y(jak ) +  s*xi     
            thread_y(jak1) = thread_y(jak1) + s1*xi
         enddo   
c
c ...    Armazena o resultado em y(i):
c
  102    thread_y(i+inc) =  t
      enddo
c$omp barrier
c
c ... Resgatando e acumulando os valores thread_y -> y
      do i = 1,nThreadsSolver
         inc = (i-1)*neq
c$omp do
         do j = thread_height(i), thread_end(i)
            y(j) = y(j) + thread_y(j+inc)
         end do
c$omp end do
      end do
c .....................................................................
c
c ... 
c$omp single          
      matvectime = get_time() - matvectime
c$omp end single
      return
      end
c *********************************************************************
c
c *********************************************************************
c * MATVEC_CSRD_SYM_ILU4_OMP: produto matriz-vetor                    *
c *                   y = Ax (A simetrica)                            *
c *                   coef. de A no formato CSR e armazenamento       *
c *                   da parte inferio, com loops aninhados.          *
c *                   (loop interno desenrolado)                      *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ad       - diagonal princiapal de A                               *
c * al(nad)  - parte triangular inferior de A no formato CSR          *
c * x        - valores do vetor a ser multiplicado por A              *
c * y        - nao definido                                           *
c * ia       - ponteiro do arranjo do csr da matrix A                 *
c * ja       - arranjo csr da matriz A(inferior)                      *
c * neq      - numero de equacoes                                     *
c * thread_y - buffer                                                 *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * y    - vetor com a operacao Ax                                    *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine matvec_csrd_sym_ilu4_omp(ad,dum,al,x,y,ia,ja,neq
     .                                   ,thread_y)
      implicit none
      include 'openmp.fi'
      include 'time.fi'
      real*8 ad(*),al(*),x(*),y(*),dum,t,thread_y(*),xi,s,s1,s2,s3
      integer ja(*),neq,i,j,k,ia(*),jak,inc,k1,k2,n,jak1,jak2,jak3
c ......................................................................
c$omp single      
      matvectime = get_time() - matvectime
c$omp end single
c
c ... inciando buffer      
c$    thread_id = omp_get_thread_num() + 1
      do i = 1, nThreadsSolver
         inc = (i-1)*neq
c$omp do
         do j = thread_height(i)+inc, thread_begin(i)+inc-1
            thread_y(j) = 0.d0
         enddo
c$omp end do
      enddo
c .....................................................................
c
c ... matVec
      inc = (thread_id - 1)*neq
c$omp barrier
      do i = thread_begin(thread_id), thread_end(thread_id)
         y(i) = 0.d0
         xi   = x(i)
c
c ...    Produto da diagonal de A por x:
c
         t  = ad(i)*xi
c
         k1 = ia(i)
         k2 = ia(i+1)
         if(k2 .eq. k1) goto 102
         n = mod(k2-k1,4)
         if(n .eq. 0) goto 101
         jak1 = ja(k1)
         t    =    t + al(k1)*x(jak1)
         jak1 = jak1 + inc
         thread_y(jak1) = thread_y(jak1) + al(k1)*xi
         k1  = k1 + 1
         if(k2 .eq. k1) goto 102
         n = n - 1
         if(n .eq. 0) goto 101
         if(n .eq. 1) then
           jak1 = ja(k1)
           t    = t + al(k1)*x(jak1)
           jak1 = jak1 + inc
           thread_y(jak1) = thread_y(jak1) + al(k1)*xi
           k1 = k1 + 1
           if (k2 .eq. k1) goto 102
         else
           jak2 = ja(k1+1)
           t = t + al(k1)*x(jak1) + al(k1+1)*x(jak2)
           jak1 = jak1 + inc
           jak2 = jak2 + inc
           thread_y(jak1) = thread_y(jak1) + al(k1)*xi
           thread_y(jak2) = thread_y(jak2) + al(k1+1)*xi
           k1 = k1 + 2
           if (k2 .eq. k1) goto 102 
         endif
c
  101    do k = k1, k2-1,4
            jak  = ja(k)
            jak1 = ja(k+1)
            jak2 = ja(k+2)
            jak3 = ja(k+3) 
            s    = al(k)
            s1   = al(k+1)
            s2   = al(k+2)
            s3   = al(k+3)
c
c ...    Parte superior de A, linha i:
c
            t = t + s*x(jak) + s1*x(jak1) + s2*x(jak2)  + s3*x(jak3)
c
c ...    Parte inferior de A, coluna i:
c
            jak  = jak  + inc
            jak1 = jak1 + inc
            jak2 = jak2 + inc
            jak3 = jak3 + inc
            thread_y(jak)  = thread_y(jak ) +  s*xi     
            thread_y(jak1) = thread_y(jak1) + s1*xi
            thread_y(jak2) = thread_y(jak2) + s2*xi
            thread_y(jak3) = thread_y(jak3) + s3*xi
         enddo   
c
c ...    Armazena o resultado em y(i):
c
  102    thread_y(i+inc) =  t
      enddo
c$omp barrier
c
c ... Resgatando e acumulando os valores thread_y -> y
      do i = 1,nThreadsSolver
         inc = (i-1)*neq
c$omp do
         do j = thread_height(i), thread_end(i)
            y(j) = y(j) + thread_y(j+inc)
         end do
c$omp end do
      end do
c .....................................................................
c
c ... 
c$omp single          
      matvectime = get_time() - matvectime
c$omp end single
      return
      end
c **********************************************************************
c
c *********************************************************************
c * MATVEC_CSRD_SYM_OLI2_ILU2_OMP: produto matriz-vetor               *
c *                   y = Ax (A simetrica)                            *
c *                   coef. de A no formato CSR e armazenamento       *
c *                   da parte inferio, com loops aninhados.          *
c *                   (loop externo/interno desenrolado)              *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ad       - diagonal princiapal de A                               *
c * al(nad)  - parte triangular inferior de A no formato CSR          *
c * x        - valores do vetor a ser multiplicado por A              *
c * y        - nao definido                                           *
c * ia       - ponteiro do arranjo do csr da matrix A                 *
c * ja       - arranjo csr da matriz A(inferior)                      *
c * neq      - numero de equacoes                                     *
c * thread_y - buffer                                                 *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * y    - vetor com a operacao Ax                                    *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine matvec_csrd_sym_ilo2_ilu2_omp(ad,dum,al,x,y,ia,ja,neq
     .                                   ,thread_y)
      implicit none
      include 'openmp.fi'
      include 'time.fi'
      real*8 ad(*),al(*),x(*),y(*),dum,t,t1,thread_y(*),xi,xi1,s,s1
      integer ja(*),neq,i,i0,j,k,ia(*),jak,inc,k1,k2,k3,nrow,n,jak1
c ......................................................................
c$omp single      
      matvectime = get_time() - matvectime
c$omp end single
c
c ... inciando buffer      
c$    thread_id = omp_get_thread_num() + 1
      do i = 1, nThreadsSolver
         inc = (i-1)*neq
c$omp do
         do j = thread_height(i)+inc, thread_begin(i)+inc-1
            thread_y(j) = 0.d0
         enddo
c$omp end do
      enddo
c .....................................................................
c
c ... matVec
      inc = (thread_id - 1)*neq
c$omp barrier
      i0   = thread_begin(thread_id)
      nrow = thread_end(thread_id) - thread_begin(thread_id) + 1
      n    = mod(nrow,2)
      if( n .eq. 0 ) goto 10
c ... produto para caso de numero de equacoes impar
      xi = x(i0)
      t  =ad(i0)*xi
      do k = ia(i0), ia(i0+1) - 1
        jak = ja(k)
        s   = al(k)
        t   = t + s*x(jak)
        jak  = jak  + inc
        thread_y(jak) = thread_y(jak) + s*xi
      enddo 
      thread_y(i0+inc) = t
      i0=i0+1
  10  do i = i0, thread_end(thread_id),2
         y(i)   = 0.d0
         y(i+1) = 0.d0
         xi    = x(i)
         xi1   = x(i+1)
c
c ...    Produto da diagonal de A por x:
c
         t  =   ad(i)*xi
         t1 = ad(i+1)*xi1
         k1  = ia(i)
         k2  = ia(i+1)
         k3  = ia(i+2)
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         if(k2 .eq. k1) goto 120
         n = mod(k2-k1,2)
         if(n .eq. 0) goto 100
         jak1 = ja(k1)
         t    =    t + al(k1)*x(jak1)
         jak1 = jak1 + inc
         thread_y(jak1) = thread_y(jak1) + al(k1)*xi
         k1 = k1 + 1
c
  100    do k = k1, k2-1,2
            jak  = ja(k)
            jak1 = ja(k+1)
            s    = al(k)
            s1   = al(k+1)
c
c ...    Parte superior de A, linha i:
c
            t = t + s*x(jak) + s1*x(jak1)
c
c ...    Parte inferior de A, coluna i:
c
            jak  = jak  + inc
            jak1 = jak1 + inc
            thread_y(jak)  = thread_y(jak ) +  s*xi     
            thread_y(jak1) = thread_y(jak1) + s1*xi
         enddo   
c
c ...    Armazena o resultado em y(i):
c
  120    thread_y(i+inc) =  t
c .....................................................................
c
c ...    Loop nos coeficientes nao nulos da linha i+1:
c
         if (k3 .eq. k2) goto 220
         n = mod(k3-k2,2)
         if (n .eq. 0)   goto 200
         jak = ja(k2) 
         t1  = t1 + al(k2)*x(jak)
         jak = jak + inc
         thread_y(jak) = thread_y(jak) + al(k2)*xi1
         k2 = k2 + 1
  200    do 210 k = k2, k3-1, 2
            jak  = ja(k)
            jak1 = ja(k+1)
            s    = al(k)
            s1   = al(k+1)
            t1   = t1 + s*x(jak) + s1*x(jak1)
            jak  = jak  + inc
            jak1 = jak1 + inc
            thread_y(jak)  = thread_y(jak)  +  s*xi1
            thread_y(jak1) = thread_y(jak1) + s1*xi1
  210    continue
  220    thread_y(i+inc+1) = t1
      enddo
c$omp barrier
c
c ... Resgatando e acumulando os valores thread_y -> y
      do i = 1,nThreadsSolver
         inc = (i-1)*neq
c$omp do
         do j = thread_height(i), thread_end(i)
            y(j) = y(j) + thread_y(j+inc)
         end do
c$omp end do
      end do
c .....................................................................
c
c ... 
c$omp single          
      matvectime = get_time() - matvectime
c$omp end single
      return
      end
c **********************************************************************
c *********************************************************************
c * MATVEC_CSRD_SYM_OLI2_ILU4_OMP: produto matriz-vetor               *
c *                   y = Ax (A simetrica)                            *
c *                   coef. de A no formato CSR e armazenamento       *
c *                   da parte inferio, com loops aninhados.          *
c *                   (loop externo/interno desenrolado)              *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ad       - diagonal princiapal de A                               *
c * al(nad)  - parte triangular inferior de A no formato CSR          *
c * x        - valores do vetor a ser multiplicado por A              *
c * y        - nao definido                                           *
c * ia       - ponteiro do arranjo do csr da matrix A                 *
c * ja       - arranjo csr da matriz A(inferior)                      *
c * neq      - numero de equacoes                                     *
c * thread_y - buffer                                                 *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * y    - vetor com a operacao Ax                                    *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine matvec_csrd_sym_ilo2_ilu4_omp(ad,dum,al,x,y,ia,ja,neq
     .                                   ,thread_y)
      implicit none
      include 'openmp.fi'
      include 'time.fi'
      real*8 ad(*),al(*),x(*),y(*),dum,t,t1,thread_y(*)
      real*8 xi,xi1,s,s1,s2,s3
      integer jak,jak1,jak2,jak3
      integer ja(*),neq,i,i0,j,k,ia(*),inc,k1,k2,k3,nrow,n
c ......................................................................
c$omp single      
      matvectime = get_time() - matvectime
c$omp end single
c
c ... inciando buffer      
c$    thread_id = omp_get_thread_num() + 1
      do i = 1, nThreadsSolver
         inc = (i-1)*neq
c$omp do
         do j = thread_height(i)+inc, thread_begin(i)+inc-1
            thread_y(j) = 0.d0
         enddo
c$omp end do
      enddo
c .....................................................................
c
c ... matVec
      inc = (thread_id - 1)*neq
c$omp barrier
      i0   = thread_begin(thread_id)
      nrow = thread_end(thread_id) - thread_begin(thread_id) + 1
      n    = mod(nrow,2)
      if( n .eq. 0 ) goto 10
c ... produto para caso de numero de equacoes impar
      xi = x(i0)
      t  =ad(i0)*xi
      do k = ia(i0), ia(i0+1) - 1
        jak = ja(k)
        s   = al(k)
        t   = t + s*x(jak)
        jak  = jak  + inc
        thread_y(jak) = thread_y(jak) + s*xi
      enddo 
      thread_y(i0+inc) = t
      i0=i0+1
  10  do i = i0, thread_end(thread_id),2
         y(i)   = 0.d0
         y(i+1) = 0.d0
         xi    = x(i)
         xi1   = x(i+1)
c
c ...    Produto da diagonal de A por x:
c
         t  =   ad(i)*xi
         t1 = ad(i+1)*xi1
         k1  = ia(i)
         k2  = ia(i+1)
         k3  = ia(i+2)
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         if(k2 .eq. k1) goto 120
         n = mod(k2-k1,4)
         if(n .eq. 0) goto 100
         jak1 = ja(k1)
         t    =    t + al(k1)*x(jak1)
         jak1 = jak1 + inc
         thread_y(jak1) = thread_y(jak1) + al(k1)*xi
         k1 = k1 + 1
         if (k2 .eq. k1) goto 120
         n  = n  - 1
         if( n.eq. 0) goto 100
         if( n .eq. 1) then
           jak1 = ja(k1)
           t = t + al(k1)*x(jak1)
           jak1 = jak1 + inc 
           thread_y(jak1) = thread_y(jak1) + al(k1)*xi
           k1 = k1 + 1
           if (k2 .eq. k1) goto 120 
         else
           jak1 = ja(k1)
           jak2 = ja(k1+1)
           t = t + al(k1)*x(jak1) + al(k1+1)*x(jak2)
           jak1 = jak1 + inc
           jak2 = jak2 + inc
           thread_y(jak1) = thread_y(jak1) + al(k1)*xi
           thread_y(jak2) = thread_y(jak2) + al(k1+1)*xi
           k1 = k1 + 2
           if (k2 .eq. k1) goto 120
         endif
c
  100    do k = k1, k2-1,4
            jak  = ja(k)
            jak1 = ja(k+1)
            jak2 = ja(k+2)
            jak3 = ja(k+3)
            s    = al(k)
            s1   = al(k+1)
            s2   = al(k+2)
            s3   = al(k+3)
c
c ...    Parte superior de A, linha i:
c
            t    = t + s*x(jak) + s1*x(jak1) + s2*x(jak2) + s3*x(jak3)
c
c ...    Parte inferior de A, coluna i:
c
            jak  = jak  + inc
            jak1 = jak1 + inc
            jak2 = jak2 + inc
            jak3 = jak3 + inc
            thread_y(jak)  = thread_y(jak)  +  s*xi
            thread_y(jak1) = thread_y(jak1) + s1*xi
            thread_y(jak2) = thread_y(jak2) + s2*xi
            thread_y(jak3) = thread_y(jak3) + s3*xi
         enddo   
c
c ...    Armazena o resultado em y(i):
c
  120    thread_y(i+inc) =  t
c .....................................................................
c
c ...    Loop nos coeficientes nao nulos da linha i+1:
c
         if (k3 .eq. k2) goto 220
           n = mod(k3-k2,4)
           if (n .eq. 0) goto 200
           jak1 = ja(k2)
           t1   = t1 + al(k2)*x(jak1)
           jak1 = jak1 + inc
           thread_y(jak1) = thread_y(jak1) + al(k2)*xi1
           k2 = k2 + 1
           if (k3 .eq. k2) goto 220
           n  = n  - 1
           if( n.eq. 0) goto 200
           if( n .eq. 1) then
             jak1 = ja(k2)
             t1 = t1 + al(k2)*x(jak1)
             jak1 = jak1 + inc 
             thread_y(jak1) = thread_y(jak1) + al(k2)*xi1
             k2 = k2 + 1
             if (k3 .eq. k2) goto 220 
           else
             jak1 = ja(k2)
             jak2 = ja(k2+1)
             t1   = t1 + al(k2)*x(jak1) + al(k2+1)*x(jak2)
             jak1 = jak1 + inc
             jak2 = jak2 + inc
             thread_y(jak1) = thread_y(jak1) + al(k2)*xi1
             thread_y(jak2) = thread_y(jak2) + al(k2+1)*xi1
             k2 = k2 + 2
             if (k3 .eq. k2) goto 220
           endif
  200      do k = k2, k3-1, 4
            jak  = ja(k)
            jak1 = ja(k+1)
            s    = al(k)
            s1   = al(k+1)
            t1   = t1 + s*x(jak) + s1*x(jak1)
            jak  = jak  + inc
            jak1 = jak1 + inc
            thread_y(jak)  = thread_y(jak)  +  s*xi1
            thread_y(jak1) = thread_y(jak1) + s1*xi1
         enddo   
  220    thread_y(i+inc+1) = t1
      enddo
c$omp barrier
c
c ... Resgatando e acumulando os valores thread_y -> y
      do i = 1,nThreadsSolver
         inc = (i-1)*neq
c$omp do
         do j = thread_height(i), thread_end(i)
            y(j) = y(j) + thread_y(j+inc)
         end do
c$omp end do
      end do
c .....................................................................
c
c ... 
c$omp single          
      matvectime = get_time() - matvectime
c$omp end single
      return
      end
c **********************************************************************
c
c **********************************************************************
      real*8 function dot_omp(a,b,n)
c **********************************************************************
c *                                                                    *
c *   DOT_OMP: Produto escalar a.b                                     *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'time.fi'
      include 'openmp.fi'
      integer n,i
      real*8  a(*),b(*)
c ......................................................................
c$omp single      
      dottime = get_time() - dottime
      omp_dot    = 0.d0
c$omp end single
c ......................................................................
c
c$omp do reduction(+:omp_dot)     
      do i = 1, n
         omp_dot = omp_dot + a(i)*b(i)
      enddo
c$omp enddo
c$omp single 
      dottime = get_time() - dottime
c$omp end single
      dot_omp = omp_dot
c$omp barrier
c ......................................................................    
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine dot_omp2(a,b,dot,n)
c **********************************************************************
c *                                                                    *
c *   DOT_OMP: Produto escalar a.b                                     *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'time.fi'
      integer n,i
      real*8  a(*),b(*),dot
c ......................................................................
c$omp single      
c      dottime = get_time() - dottime
      dot     = 0.d0
c$omp end single
c ......................................................................
c
c$omp do reduction(+:dot)     
      do i = 1, n
         dot = dot + a(i)*b(i)
      enddo
c$omp enddo
c$omp single 
c      dottime = get_time() - dottime
c$omp end single
c$omp barrier
c ......................................................................    
      return
      end
c **********************************************************************