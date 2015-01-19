c *********************************************************************
c * MATVEC_CSRD: produto matriz vetor no formato csr ( grafo de A     *
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
      subroutine matvec_csrd(ad,a,dum,x,y,ia,ja,neq)
      implicit none
      include 'time.fi'
      real*8 ad(*),a(*),x(*),y(*),dum,t
      integer ja(*),neq,i,j,ia(*)
c ......................................................................
      matvectime = get_time() - matvectime
      do i = 1, neq
c ... diagonal principal
        t = ad(i)*x(i)
c ... produto da linha
        do j = ia(i) , ia(i+1) - 1
          t  = t + a(j)*x(ja(j))
        enddo
        y(i) = t
      enddo
      matvectime = get_time() - matvectime
      return
      end
c **********************************************************************
c
c *********************************************************************
c * MATVEC_CSR: produto matriz vetor no formato csr ( grafo de A      *
c * nao simetrico; coeficiente de A nao simetrico; a - inferior,      *
c * superior e diagonal)                                              *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * a      - matriz no formato CSR                                    *
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
      subroutine matvec_csr(a,dum1,dum2,x,y,ia,ja,neq)
      implicit none
      include 'time.fi'
      real*8 a(*),x(*),y(*),dum1,dum2,t
      integer ja(*),neq,i,j,ia(*)
c ......................................................................
       matvectime = get_time() - matvectime
       do i = 1, neq
         t = 0.0d0
c ... produto da linha
         do j = ia(i) , ia(i+1) - 1
           t  = t + a(j)*x(ja(j))
         enddo
         y(i) = t
       enddo
      matvectime = get_time() - matvectime
      return
      end
c **********************************************************************
c
c *********************************************************************
c * MATVEC_CSRC: produto matriz vetor no formato csrc( grafo de A     *
c * simetrico; coeficiente de A nao simetrico; a - inferior,          *
c * superior e diagonal)                                              *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ad     - diagonal da matriz A                                     *
c * x      - valores do vetor a ser multiplicado por A                *
c * y      - nao definido                                             *  
c * al(nad)- parte triangular inferior de A no formato CSR            *
c * au(*)  - parte triangular superior de A no formato CSC            *                                             *
c * ia     - ponteiro do arranjo do csr da matrix A                   *
c * ja     - arranjo csr da matriz A                                  *
c * la     - valores locais da matriz A da celula P                   *
c * neq    - numero de equacoes                                       *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * y    - vetor com a operacao Ax                                    *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine matvec_csrc(ad,au,al,x,y,ia,ja,neq)
      implicit none
      include 'time.fi'
      real*8 ad(*),x(*),y(*),al(*),au(*),t,xi
      integer ja(*),neq,i,k,ia(*),jak
c ......................................................................
       matvectime = get_time() - matvectime
       do i = 1, neq
c ...    Produto da diagonal de A por x:
c
         xi = x(i)
         t  = ad(i)*xi
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         do k = ia(i), ia(i+1)-1
            jak = ja(k)
c
c ...       Produto da linha i pelo vetor x (parte triangular inferior):
c
            t   = t + al(k)*x(jak)
c
c ...       Produto dos coef. da parte triangular superior por x(i):
c
            y(jak) = y(jak) + au(k)*xi
         enddo
c
c ...    Armazena o resultado em y(i):
c
         y(i) = t
      enddo         
      matvectime = get_time() - matvectime
      return
      end
c **********************************************************************
c
c *********************************************************************
c * MATVEC_CSRC_SYM: produto matriz vetor no formato csrc( grafo de A *
c * simetrico; coeficiente de A simetrico; ad-diagonal;al-inferior)   *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ad     - diagonal da matriz A                                     *
c * x      - valores do vetor a ser multiplicado por A                *
c * y      - nao definido                                             *  
c * al(nad)- parte triangular inferior de A no formato CSR            *
c * au(*)  - nao utilizado                                            *  
c * ia     - ponteiro do arranjo do csr da matrix A                   *
c * ja     - arranjo csr da matriz A                                  *
c * la     - valores locais da matriz A da celula P                   *
c * neq    - numero de equacoes                                       *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * y    - vetor com a operacao Ax                                    *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine matvec_csrc_sym(ad,dum,al,x,y,ia,ja,neq)
      implicit none
      include 'time.fi'
      real*8 ad(*),x(*),y(*),al(*),dum,t,xi,s
      integer ja(*),neq,i,k,ia(*),jak
c ......................................................................
       matvectime = get_time() - matvectime
       y(1) = ad(1)*x(1)
       do i = 2, neq
c ...    Produto da diagonal de A por x:
c
         xi = x(i)
         t  = ad(i)*xi
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         do k = ia(i), ia(i+1)-1
            jak = ja(k)
            s   = al(k)
c
c ...       Produto da linha i pelo vetor x (parte triangular inferior):
c
            t   = t + s*x(jak)
c
c ...       Produto dos coef. da parte triangular superior por x(i):
c
            y(jak) = y(jak) + s*xi
         enddo
c
c ...    Armazena o resultado em y(i):
c
         y(i) = t
      enddo         
      matvectime = get_time() - matvectime
      return
      end
c **********************************************************************
c
c *********************************************************************
c  * MATVEC_CSRD_SYM: produto matriz-vetor y = Ax  (A simetrica),     *
c *                   coef. de A no formato CSR e armazenamento       *
c *                   da parte superior, com loops aninhados.         *                                  
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ad     - diagonal da matriz A                                     *
c * x      - valores do vetor a ser multiplicado por A                *
c * y      - nao definido                                             *  
c * al(nad)- parte triangular inferior de A no formato CSR            *
c * dum    - nao utilizado                                            * 
c * ia     - ponteiro do arranjo do csr da matrix A                   *
c * ja     - arranjo csr da matriz A                                  *
c * la     - valores locais da matriz A da celula P                   *
c * neq    - numero de equacoes                                       *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * y    - vetor com a operacao Ax                                    *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine matvec_csrd_sym(ad,dum,al,x,y,ia,ja,neq)
      implicit none
      include 'time.fi'
      real*8 ad(*),x(*),y(*),al(*),dum,t,xi,s
      integer ja(*),neq,i,k,ia(*),jak
c ......................................................................
      matvectime = get_time() - matvectime
      y(1) = ad(1)*x(1)
      do i = 2, neq
         xi = x(i)
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
            y(jak) = y(jak) + s*xi
         enddo   
c
c ...    Armazena o resultado em y(i):
c
         y(i) = t
      enddo    
      matvectime = get_time() - matvectime
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine vet(a,b,c)
c **********************************************************************
c *                                                                    *
c *   VET: Produto vetorial axb                                        *
c *                                                                    *
c **********************************************************************
      implicit none
      real*8  c(3),a(3),b(3)
c ......................................................................
      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
c      print*,c(1),c(2),c(3)
c ......................................................................    
      return
      end
c **********************************************************************
c
c *********************************************************************
      subroutine vsmul(a,x,n,c)
c **********************************************************************
c *                                                                    *
c *   VSUM:                                                            *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    a(n) - vetor de dimensao n                                      *
c *    x    - escalar                                                  *
c *    n    - dimensao                                                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   c = x*a                                                          *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,i
      real*8 a(*),x,c(*)
c ......................................................................
      do 100 i = 1, n
         c(i) = a(i)*x
  100 continue
      return
      end  
c *********************************************************************
c
c *********************************************************************
      subroutine asumb(a,b,n,c)
c **********************************************************************
c *                                                                    *
c *   VSUM:                                                            *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    a(n) - vetor de dimensao n                                      *
c *    x    - escalar                                                  *
c *    n    - dimensao                                                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   c = x*a                                                          *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,i
      real*8 a(*),b(*),c(*)
c ......................................................................
      do 100 i = 1, n
         c(i) = a(i)+b(i)
  100 continue
      return
      end  
c *********************************************************************
      real*8 function dot(a,b,n)
c **********************************************************************
c *                                                                    *
c *   DOT: Produto escalar a.b                                         *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'time.fi'
      integer n,i
      real*8  a(*),b(*)
c ......................................................................
      dottime = get_time() - dottime
      dot = 0.d0
      do i = 1, n
         dot = dot + a(i)*b(i)
      enddo
      dottime = get_time() - dottime
c ......................................................................    
      return
      end
c **********************************************************************
      real*8 function dot_local(a,b,n)
c **********************************************************************
c *                                                                    *
c *   DOT_LOCAL: Produto escalar a.b para vetores pequenos             *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'time.fi'
      integer n,i
      real*8  a(*),b(*)
c ......................................................................
      dot_local = 0.d0
      do i = 1, n
         dot_local = dot_local + a(i)*b(i)
      enddo
c ......................................................................    
      return
      end
c **********************************************************************
c
      subroutine aequalb(a,b,n)
c **********************************************************************
c *                                                                    *
c *   AEQUALB:                                                         *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    a(n) - vetor de dimensao n                                      *
c *    b(n) - vetor de dimensao n                                      *
c *    n    - dimensao                                                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   a = b                                                            *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,i
      real*8 a(*),b(*)
c ......................................................................
      do 100 i = 1, n
         a(i) = b(i)
  100 continue
      return
      end
c **********************************************************************
c
c **********************************************************************
      subroutine aequalbVetor(a,b,n,ndf)
c **********************************************************************
c *                                                                    *
c *   AEQUALB:                                                         *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    a(n) - vetor de dimensao n                                      *
c *    b(n) - vetor de dimensao n                                      *
c *    n    - dimensao                                                 *
c *    ndf  - graus de liberdade                                       *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   a = b                                                            *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,ndf,i,j
      real*8 a(ndf,*),b(ndf,*)
c ......................................................................
      do i = 1, n
        do j = 1, ndf
          a(j,i) = b(j,i)
        enddo 
      enddo
      return
      end
c **********************************************************************
c 
