c *********************************************************************
c * CORRECT: corretor                                                 *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * u      - nao definido                                             *
c * x      - solucao do sistema  (K*du=R)                             *
c * numel  - numero de elementos                                      *
c * numm   - arranjo com os valores da equções da celula i            *
c * ndf    - graus de liberdade                                       *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c *  u     - valores atualizado para iteracao i+1                     *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine correct(u,x,num,numel,ndf)
      implicit none
      integer neq,ndf,i,j,numel,num(*)
      real*8 u(ndf,*),x(*)
      do i = 1, numel
        neq = num(i)
        do j = 1 , ndf  
          u(j,i) =  u(j,i)  +  x(neq)
        enddo      
      enddo
      return
      end
c *********************************************************************
c
c *********************************************************************
c * PREDICT: preditor                                                 *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * u      - valores da iteracao i                                    *
c * du     -                                                          *
c * neq    - numero de equacoes                                       *
c * ndf    - graus de liberdade                                       *
c * alfa   - para metro da tecnica utilizada                          *
c * dt     - passo de tempo                                           *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c *  u     - valores atualizado para iteracao i+1                     *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine predict(util,u,v,neq,ndf,alpha,dt)
      implicit none
      integer neq,ndf,i,j
      real*8 util(ndf,*),u(ndf,*),v(ndf,*),b,dt,alpha
      b = (1.0d0-alpha)*dt
      do i = 1, neq
        do j = 1 , ndf  
          util(j,i) = u(j,i) +  v(j,i)*b
          v(j,i)    = 0.0d0               
        enddo      
      enddo
      return
      end
