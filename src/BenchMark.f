c *********************************************************************
c * BENCHMARKDOT : seleciona o melhor produto escalar para o problema *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * openmp - flag do opnemp                                           *
c * neq    - diagonal principal de A                                  *
c * nThreads - numero de threads                                      *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine benchmarkDot(openmp,neq,nThreads)
      use Malloc
      implicit none
      include 'time.fi'
      include 'openmp.fi'
      integer neq,i,j,nThreads,mnt
      integer,parameter :: nsample=100000
      real*8 ,parameter :: ZERO=1.0d-15
      logical openmp
      real*8 timei
      real*8 d
      real*8 sumSquares
      integer*8 i_b1,i_b2 
c ... sequencial
      real*8 dot
c ... omp       
      real*8 dot_omp
      real*8 dot_ompL2,dot_ompL4,dot_ompL6,dot_ompL8
      real*8 dot_ompO2,dot_ompO4
      real*8 dot_ompO2L2
c ... alocando memoria
      i_b1  = alloc_8('bench1  ',             1,neq)
      i_b2  = alloc_8('bench2  ',             1,neq)
c .....................................................................
c
c ...
      call initDot(ia(i_b1),ia(i_b2),neq)
c .....................................................................
c
c ... sequencial
      timei = 0.0d0
      do i = 1, nsample
        timei = get_time() - timei
        d = dot(ia(i_b1),ia(i_b2),neq)
        timei = get_time() - timei
        if(dabs(d - sumSquares(neq)) .gt. ZERO) then
          print*,'benchmarkDot: dot nao passou no teste!!!.'
          stop
        endif
      enddo
      write(*,2000),'dot',timei,(d-sumSquares(neq))
c .....................................................................  
c 
c ... omp
      if(openmp) then
c ... dot_omp
c$      mnt = omp_get_max_threads()
        do j = 2, mnt,2
          timei = 0.0d0
          do i = 1, nsample
            timei = get_time() - timei
c$omp parallel num_threads(j)
            d = dot_omp(ia(i_b1),ia(i_b2),neq)
c$omp end parallel
            timei = get_time() - timei
            if(dabs(d - sumSquares(neq)) .gt. ZERO) then
              print*,'benchmarkDot: dot_omp nao passou no teste!!!.'
              stop
            endif
          enddo
          write(*,3000),'dot_omp',j,timei,(d-sumSquares(neq))
        enddo
c .....................................................................  
c
c ... dot_ompL2
        do j = 2, mnt,2
          timei = 0.0d0
          do i = 1, nsample
            timei = get_time() - timei
c$omp parallel num_threads(j)
            d = dot_ompL2(ia(i_b1),ia(i_b2),neq)
c$omp end parallel
            timei = get_time() - timei
            if(dabs(d - sumSquares(neq)) .gt. ZERO) then
              print*,'benchmarkDot: dot_ompL2 nao passou no teste!!!.'
              stop
            endif
          enddo
          write(*,3000),'dot_ompL2',j,timei,(d-sumSquares(neq))
        enddo
c .....................................................................  
c
c ... dot_ompL4
        do j = 2, mnt,2
          timei = 0.0d0
          do i = 1, nsample
            timei = get_time() - timei
c$omp parallel num_threads(j)
            d = dot_ompL4(ia(i_b1),ia(i_b2),neq)
c$omp end parallel
            timei = get_time() - timei
            if(dabs(d - sumSquares(neq)) .gt. ZERO) then
              print*,'benchmarkDot: dot_ompL4 nao passou no teste!!!.'
              stop
            endif
          enddo
          write(*,3000),'dot_ompL4',j,timei,(d-sumSquares(neq))
       enddo 
c .....................................................................  
c
c ... dot_ompL6
        do j = 2, mnt,2
          timei = 0.0d0
          do i = 1, nsample
            timei = get_time() - timei
c$omp parallel num_threads(j)
            d = dot_ompL6(ia(i_b1),ia(i_b2),neq)
c$omp end parallel
            timei = get_time() - timei
            if(dabs(d - sumSquares(neq)) .gt. ZERO) then
              print*,'benchmarkDot: dot_ompL6 nao passou no teste!!!.'
              stop
            endif
          enddo
          write(*,3000),'dot_ompL6',j,timei,(d-sumSquares(neq))
       enddo 
c .....................................................................  
c
c ... dot_ompL8
        do j = 2, mnt,2
          timei = 0.0d0
          do i = 1, nsample
            timei = get_time() - timei
c$omp parallel num_threads(j)
            d = dot_ompL8(ia(i_b1),ia(i_b2),neq)
c$omp end parallel
            timei = get_time() - timei
            if(dabs(d - sumSquares(neq)) .gt. ZERO) then
              print*,'benchmarkDot: dot_ompL8 nao passou no teste!!!.'
              stop
            endif
          enddo
          write(*,3000),'dot_ompL8',j,timei,(d-sumSquares(neq))
       enddo 
c .....................................................................  
c
c ... dot_ompO2
        do j = 2, mnt,2
          timei = 0.0d0
          do i = 1, nsample
            timei = get_time() - timei
c$omp parallel num_threads(j)
            d = dot_ompO2(ia(i_b1),ia(i_b2),neq)
c$omp end parallel
            timei = get_time() - timei
            if(dabs(d - sumSquares(neq)) .gt. ZERO) then
              print*,'benchmarkDot: dot_ompO2 nao passou no teste!!!.'
              stop
            endif
          enddo
          write(*,3000),'dot_ompO2',j,timei,(d-sumSquares(neq))
       enddo 
c .....................................................................  
c
c ... dot_ompO4
        do j = 2, mnt,2
          timei = 0.0d0
          do i = 1, nsample
            timei = get_time() - timei
c$omp parallel num_threads(j)
            d = dot_ompO4(ia(i_b1),ia(i_b2),neq)
c$omp end parallel
            timei = get_time() - timei
            if(dabs(d - sumSquares(neq)) .gt. ZERO) then
              print*,'benchmarkDot: dot_ompO4 nao passou no teste!!!.'
              stop
            endif
          enddo
          write(*,3000),'dot_ompO4',j,timei,(d-sumSquares(neq))
       enddo 
c .....................................................................  
c
c ... dot_ompO2L2
        do j = 2, mnt,2
          timei = 0.0d0
          do i = 1, nsample
            timei = get_time() - timei
c$omp parallel num_threads(j)
            d = dot_ompO2L2(ia(i_b1),ia(i_b2),neq)
c$omp end parallel
            timei = get_time() - timei
            if(dabs(d - sumSquares(neq)) .gt. ZERO) then
              print*,'benchmarkDot: dot_ompO2L2 nao passou no teste!!!.'
              stop
            endif
          enddo
          write(*,3000),'dot_ompO2L2',j,timei,(d-sumSquares(neq))
        enddo
c .....................................................................  
c 
c ... 
      endif 
      dottime = 0.0d0
c .....................................................................  
c 
c ... liberacao da  memoria
      i_b1  = dealloc('bench1  ')
      i_b2  = dealloc('bench2  ')
c .....................................................................
      return
 2000 format(a12,1x,'time :',1x,f20.3,1x,'erro :',d22.15)
 3000 format(a12,1x,'nthreads :',i2,1x,': time :'
     .      ,f20.3,1x,'erro :',d22.15)
      end
c *********************************************************************
c
c *********************************************************************
c * BENCHMARKDOT : seleciona o melhor produto escalar para o problema *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * openmp   - flag do opnemp                                         *
c * unsym    - matrix nao simetrica                                   *
c * iia      - apontador do csr                                       *
c * ja       - apontador do csr                                       *
c * ja       - apontador do csr                                       *
c * al       - coeficientes da parte inferior                         *
c * ad       - coeficientes da diagonal                               *
c * au       - coeficientes da parte superior                         *
c * neq      - numero de equacoes                                     *
c * nad      - numero nao nulos fora da diagonal principal            *
c * nThreads - numero de threads                                      *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine benchmarkCsr(openmp,unsym,iia,ja,al,ad,au,neq,nad
     .                       ,nThreads)
      use Malloc
      implicit none
      include 'time.fi'
      integer neq,nad,i,nThreads
      integer,parameter :: nsample=10000
      integer iia(*),ja(*)
      real*8 al(*),ad(*),au(*)
      logical openmp,unsym
      real*8 timei
      real*8 dot,d
      integer*8 i_x,i_y,i_b 
c ... alocando memoria
      i_x   = alloc_8('benchx  ',             1,neq)
      i_y   = alloc_8('benchy  ',             1,neq)
c .....................................................................
c
c ...
      call initMatVec(al,ad,au,ia(i_x),neq,nad,unsym)
c .....................................................................
c
c ... sequencial
      timei = 0.0d0
      if(unsym) then
        do i = 1, nsample
          timei = get_time() - timei
          call matvec_csrd(ad,al,au,ia(i_x),ia(i_y),iia,ja,neq)
          timei = get_time() - timei
        enddo
        d = dsqrt(dot(ia(i_y),ia(i_y),neq))
        write(*,2000),'matvec_csrd',timei,d
      else
        do i = 1, nsample
          timei = get_time() - timei
          call matvec_csrd_sym(ad,au,al,ia(i_x),ia(i_y),iia,ja,neq)
          timei = get_time() - timei
        enddo
        d = dsqrt(dot(ia(i_y),ia(i_y),neq))
        write(*,2000),'matvec_csrd_sym',timei,d
      endif
c .....................................................................  
c 
c ... omp
      if(openmp) then
c ... buffer
        i_b  = alloc_8('buffer  ',nThreads,neq)
        call partition_matrix(iia,ja,neq)
c .....................................................................
c
c ... matvec_csrd_omp
        timei = 0.0d0
        if(unsym) then
          do i = 1, nsample
            timei = get_time() - timei
c$omp parallel num_threads(nThreads)
            call matvec_csrd_omp(ad,al,au,ia(i_x),ia(i_y)
     .                          ,iia,ja,neq)
c$omp end parallel
            timei = get_time() - timei
          enddo
          d = dsqrt(dot(ia(i_y),ia(i_y),neq))
          write(*,2000),'matvec_csrd_omp',timei,d
        else
          do i = 1, nsample
            timei = get_time() - timei
c$omp parallel num_threads(nThreads)
            call matvec_csrd_sym_omp(ad,au,al,ia(i_x),ia(i_y)
     .                              ,iia,ja,neq,ia(i_b))
c$omp end parallel
            timei = get_time() - timei
          enddo
          d = dsqrt(dot(ia(i_y),ia(i_y),neq))
          write(*,2000),'matvec_csrd_sym_omp',timei,d
        endif
c .....................................................................
c
c ... matvec_csrd_ilu2_omp
        timei = 0.0d0
        if(unsym) then
          do i = 1, nsample
            timei = get_time() - timei
c$omp parallel num_threads(nThreads)
            call matvec_csrd_ilu2_omp(ad,al,au,ia(i_x),ia(i_y)
     .                            ,iia,ja,neq)
c$omp end parallel
            timei = get_time() - timei
          enddo
          d = dsqrt(dot(ia(i_y),ia(i_y),neq))
          write(*,2000),'matvec_csrd_ilu2_omp',timei,d
        else
          do i = 1, nsample
            timei = get_time() - timei
c$omp parallel num_threads(nThreads)
            call matvec_csrd_sym_ilu2_omp(ad,au,al,ia(i_x),ia(i_y)
     .                                   ,iia,ja,neq,ia(i_b))
c$omp end parallel
            timei = get_time() - timei
          enddo
          d = dsqrt(dot(ia(i_y),ia(i_y),neq))
          write(*,2000),'matvec_csrd_sym_ilu2_omp',timei,d
        endif
c .....................................................................
c
c ... matvec_csrd_ilu4_omp
        timei = 0.0d0
        if(unsym) then
          do i = 1, nsample
            timei = get_time() - timei
c$omp parallel num_threads(nThreads)
            call matvec_csrd_ilu4_omp(ad,al,au,ia(i_x),ia(i_y)
     .                            ,iia,ja,neq)
c$omp end parallel
            timei = get_time() - timei
          enddo
          d = dsqrt(dot(ia(i_y),ia(i_y),neq))
          write(*,2000),'matvec_csrd_ilu4_omp',timei,d
        else
          do i = 1, nsample
            timei = get_time() - timei
c$omp parallel num_threads(nThreads)
            call matvec_csrd_sym_ilu4_omp(ad,au,al,ia(i_x),ia(i_y)
     .                                   ,iia,ja,neq,ia(i_b))
c$omp end parallel
            timei = get_time() - timei
          enddo
          d = dsqrt(dot(ia(i_y),ia(i_y),neq))
          write(*,2000),'matvec_csrd_sym_ilu4_omp',timei,d
        endif
c .....................................................................
c
c ... matvec_csrd_ilo2_ilu2_omp
        timei = 0.0d0
        if(unsym) then
          do i = 1, nsample
            timei = get_time() - timei
c$omp parallel num_threads(nThreads)
            call matvec_csrd_ilo2_ilu2_omp(ad,al,au,ia(i_x),ia(i_y)
     .                            ,iia,ja,neq)
c$omp end parallel
            timei = get_time() - timei
          enddo
          d = dsqrt(dot(ia(i_y),ia(i_y),neq))
          write(*,2000),'matvec_csrd_ilo2_ilu2_omp',timei,d
        else
          do i = 1, nsample
            timei = get_time() - timei
c$omp parallel num_threads(nThreads)
            call matvec_csrd_sym_ilo2_ilu2_omp(ad,au,al,ia(i_x),ia(i_y)
     .                                        ,iia,ja,neq,ia(i_b))
c$omp end parallel
            timei = get_time() - timei
          enddo
          d = dsqrt(dot(ia(i_y),ia(i_y),neq))
          write(*,2000),'matvec_csrd_sym_ilo2_ilu2_omp',timei,d
        endif
c .....................................................................
c
c ... matvec_csrd_ilo2_ilu4_omp
        timei = 0.0d0
        if(unsym) then
        else
          do i = 1, nsample
            timei = get_time() - timei
c$omp parallel num_threads(nThreads)
            call matvec_csrd_sym_ilo2_ilu4_omp(ad,au,al,ia(i_x),ia(i_y)
     .                                        ,iia,ja,neq,ia(i_b))
c$omp end parallel
            timei = get_time() - timei
          enddo
          d = dsqrt(dot(ia(i_y),ia(i_y),neq))
          write(*,2000),'matvec_csrd_sym_ilo2_ilu4_omp',timei,d
        endif
c .....................................................................
      endif
c .....................................................................  
c 
c ... 
      matvectime = 0.0d0
c .....................................................................  
c 
c ... liberacao da  memoria
      i_b   = dealloc('buffer  ')
      i_y   = dealloc('benchy  ')
      i_x   = dealloc('benchx  ')
c .....................................................................
      return
 2000 format(a32,1x,':',1x,f20.3,1x,'erro :',d22.15)
      end
c *********************************************************************
c
c *********************************************************************
      subroutine initDot(x,y,neq)
      integer neq,i
      real*8 x(*),y(*)
      do i = 1, neq
        x(i) = i
        y(i) = 2.d0*i
      enddo
      return
      end
c *********************************************************************
c
c *********************************************************************
      subroutine initMatVec(al,ad,au,x,neq,nad,unsym)
      integer neq,nad,i
      logical unsym
      real*8 x(*),al(*),ad(*),au(*)
      if(unsym) then
        ad(1:neq) =  1.0d0
        al(1:nad) = -0.1d0
        au(1:nad) = -0.1d0
      else
        ad(1:neq) =  1.0d0
        al(1:nad) = -0.1d0
      endif
      do i = 1, neq
        if( mod(i,2) .eq. 0 ) then
          x(i) =  i/1000.0d0
        else
          x(i) = -2.0d0*i/1000.0d0
        endif  
      enddo
      return
      end
c *********************************************************************
c
c *********************************************************************
      real*8 function  sumSquares(x)
      integer x
      sumSquares = x*(x+1.0d0)*(2.0d0*x+1.0d0)/6.0d0 
      sumSquares = 2.0d0*sumSquares
      return
      end
c *********************************************************************
