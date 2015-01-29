c *********************************************************************
c * SOLVER : Metodos de solucao para o sistema linear (Ax=b)          *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ad     - diagonal principal de A                                  *
c * au     - coeficiantes fora da diagonal principal de A             *
c * al     - coeficiantes fora da diagonal principal de A             *
c * x      - nao definido                                             *
c * x      - nao definido                                             *
c * b      - termo indepedente                                        *
c * iax    - ponteiro do arranjo do csr da matrix A                   *
c * ja     - arranjo csr da matriz A                                  *
c * neq    - numero de equacoes                                       *
c * nad    - numero de elementos nao nulos fora da diagonal principal *
c * unsym  - matriz armazenada na forma nao simenteria                *
c * tol    - tolerancia do solver                                     *
c * maxit  - numero maximo de iteracao do solver                      *
c * code   - 1 CG                                                     *     
c *          2 BiCGSTAB                                               *
c *          3 GMRES                                                  *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * x    - solucao do sistema                                         *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine solver(ad,au,al,x,b,iax,ja,neq,nad,unsym,tol,maxit
     .                 ,code,matriz,nlit)
      Use Malloc
      implicit none
      include 'openmp.fi'
      include 'time.fi'
      include 'error.fi'
c ... ponteiros
      logical unsym
      integer*8 i_r,i_z,i_m,i_t,i_v,i_p,i_y
      real*8 ad(*),al(*),au(*),x(*),b(*),tol
      integer iax(*),ja(*),maxit,neq,nad,k,code,matriz
      integer nlit
c ... nao simetricos Csr
      external matvec_csrd,matvec_csr,matvec_csrc
      external matvec_csrd_omp
      external matvec_csrd_ilu2_omp,matvec_csrd_ilu4_omp
      external matvec_csrd_ilo2_ilu2_omp
c ... simetricos Csr  
      external matvec_csrd_sym
      external matvec_csrd_sym_omp
      external matvec_csrd_sym_ilu2_omp
      external matvec_csrd_sym_ilu4_omp
      external matvec_csrd_sym_ilo2_ilu2_omp
      external matvec_csrd_sym_ilo2_ilu4_omp
      external matvec_csrc_sym
c ... Ellpac
      external matvec_ellpack
      external matvec_ellpack_omp
c ... dot
      external dot
      external dot_omp,dot_omp2
      external dot_ompL2,dot_ompL4,dot_ompL6,dot_ompL8
      external dot_ompO2,dot_ompO4,dot_ompO6
      external dot_ompO2L2
      k = 15
c ... buffer para Csr simetrico
      i_threads_y = 1 
      if(openmpSolver .and. (.not. unsym)) then  
        i_threads_y  = alloc_8('buffer_y',nThreadsSolver,neq)
        call partition_matrix(iax,ja,neq)
      endif
c ....................................................................
c
c ... PCG
      if( code .eq. 1) then
c       if(unsym) then
c         write(*,2000)'Solver','Solver.f'
c    .         ,'PCG nao trabalha com matriz nao simetrica'
c         stop
c       endif
        i_r         = alloc_8('r       ',  1,neq)
        i_z         = alloc_8('z       ',  1,neq)
        i_m         = alloc_8('m       ',  1,neq)
c ... precondicionador diagonal 
        precontime = get_time() - precontime         
        call aequalb(ia(i_m),ad,neq)
        precontime = get_time() - precontime 
c .....................................................................
c
c ... CSRD
        if(matriz .eq. 1) then
c ... OPENMP
          if(openmpSolver) then
c ...
            if(unsym) then
              call pcg_omp(neq,nad,iax,ja,ad,au,al,ia(i_m),b,x,ia(i_z)
c ... comum     
     .                   ,ia(i_r),tol,maxit,matvec_csrd_omp
c ... loop interno desenrolado 2     
c    .                  ,ia(i_r),tol,maxit,matvec_csrd_ilu2_omp
c ... loop interno desenrolado 4      
c    .                  ,ia(i_r),tol,maxit,matvec_csrd_ilu4_omp
c ... loop externo desenrolado 2 / loop interno desenrolado 2      
c    .                  ,ia(i_r),tol,maxit,matvec_csrd_ilo2_ilu2_omp
c ... loop externo desenrolado 2 / loop interno desenrolado 4      
     .                  ,dot_ompL6,ia(i_threads_y),nlit,.false.)
c ....................................................................
c
c ...
            else
              call pcg_omp(neq,nad,iax,ja,ad,au,al,ia(i_m),b,x,ia(i_z)
c ... comum     
     .                    ,ia(i_r),tol,maxit,matvec_csrd_sym_omp
c ... loop interno desenrolado 2     
c    .                    ,ia(i_r),tol,maxit,matvec_csrd_sym_ilu2_omp
c ... loop interno desenrolado 4      
c    .                    ,ia(i_r),tol,maxit,matvec_csrd_sym_ilu4_omp
c ... loop externo desenrolado 2 / loop interno desenrolado 2      
c    .                    ,ia(i_r),tol,maxit,matvec_csrd_sym_ilo2_ilu2_omp
c ... loop externo desenrolado 2 / loop interno desenrolado 4      
c    .                    ,ia(i_r),tol,maxit,matvec_csrd_sym_ilo2_ilu4_omp
     .                    ,dot_ompL4,ia(i_threads_y),nlit,.false.)
c .....................................................................
            endif
c .....................................................................
c
c ... sequencial          
          else
            if(unsym) then
              call pcg(neq,nad,iax,ja,ad,au,al,ia(i_m),b,x,ia(i_z)
     .                ,ia(i_r),tol,maxit,matvec_csrd    
     .                ,dot,nlit,.true.)
            else
              call pcg(neq,nad,iax,ja,ad,au,al,ia(i_m),b,x,ia(i_z)
     .                ,ia(i_r),tol,maxit,matvec_csrd_sym
     .                ,dot,nlit,.true.)
            endif
c .....................................................................
          endif
c .....................................................................
c
c ... CSR
        elseif(matriz .eq. 2) then
          call pcg(neq,nad,iax,ja,ad,au,al,ia(i_m),b,x,ia(i_z),ia(i_r)
     .            ,tol,maxit,matvec_csr,nlit,.true.)
c .....................................................................
c
c ... CSRC
        elseif(matriz .eq. 3) then
          call pcg(neq,nad,iax,ja,ad,au,al,ia(i_m),b,x,ia(i_z)
     .            ,ia(i_r),tol,maxit,matvec_csrc_sym,nlit,.true.)
c .....................................................................
c
c ... ELLPACK
        elseif(matriz .eq. 4) then
c ... OPENMP
          if(openmpSolver) then
              call pcg_omp(neq,nad,iax,ja,ad,au,al,ia(i_m),b,x,ia(i_z)
c ... comum     
     .                    ,ia(i_r),tol,maxit,matvec_ellpack_omp
     .                    ,dot_ompL4,ia(i_threads_y),nlit,.false.)
c .....................................................................
c
c ...
          else
            call pcg(neq,nad,iax,ja,ad,au,al,ia(i_m),b,x,ia(i_z)
     .              ,ia(i_r),tol,maxit,matvec_ellpack,dot,nlit,.false.)
          endif
c .....................................................................
        endif
c .....................................................................
        i_m         = dealloc('m       ')
        i_z         = dealloc('z       ')
        i_r         = dealloc('r       ')
c .....................................................................
c
c ... PBiCGSTAB
      elseif(code .eq. 2) then
c ...
        i_r         = alloc_8('r       ',  1,neq)
        i_z         = alloc_8('z       ',  1,neq)
        i_m         = alloc_8('m       ',  1,neq)
        i_t         = alloc_8('t       ',  1,neq)
        i_v         = alloc_8('vi      ',  1,neq)
        i_p         = alloc_8('pp      ',  1,neq)
c ... precondicionador diagonal
        precontime = get_time() - precontime   
        call aequalb(ia(i_m),ad,neq)
        precontime = get_time() - precontime 
c .....................................................................
c
c ...
        if(matriz .eq. 1) then
c ... CSRD - coeficiente da matriz nao simetrica
          if(unsym) then
c ... OPENMP
            if(openmpSolver) then                      
              call pbicgstab_omp(neq,nad,iax,ja,ad,au,al,ia(i_m),b,x
     .                          ,ia(i_t),ia(i_v),ia(i_r),ia(i_p),ia(i_z)
     .                          ,tol,maxit
c ... matvec comum
c    .                          ,matvec_csrd_omp
c ... loop interno desenrolado 2   
     .                          ,matvec_csrd_ilu2_omp
c ... loop interno desenrolado 4   
c    .                          ,matvec_csrd_ilu4_omp
c ... loop externo desenrolado 2 / loop interno desenrolado 2  
c     .                          ,matvec_csrd_ilo2_ilu2_omp
     .                          ,dot_ompL4,nlit,.true.)
c .....................................................................
c
c ... sequencial
            else 
              call pbicgstab(neq,nad,iax,ja,ad,au,al,ia(i_m),b,x
     .                      ,ia(i_t),ia(i_v),ia(i_r),ia(i_p),ia(i_z)
     .                      ,tol,maxit,matvec_csrd,dot,nlit,.true.)
            endif
c .....................................................................
c
c ... CSRD - coeficiente da matriz simetrica
          else
            call pbicgstab(neq,nad,iax,ja,ad,au,al,ia(i_m),b,x
     .                    ,ia(i_t),ia(i_v),ia(i_r),ia(i_p),ia(i_z)
     .                    ,tol,maxit,matvec_csrd_sym,dot,nlit,.true.)
          endif 
        elseif(matriz .eq. 2) then
          write(*,2000)'Solver','Solver.f'
     .         ,'PBiCGSTAB nao implementado para o CSR'
          stop  
c .....................................................................
c
c ...
        elseif(matriz .eq. 3) then
c ... CSRC - coeficiente da matriz nao simetrica
           if(unsym) then 
             call pbicgstab(neq,nad,iax,ja,ad,au,al,ia(i_m),b,x
     .                     ,ia(i_t),ia(i_v),ia(i_r),ia(i_p),ia(i_z)
     .                     ,tol,maxit,matvec_csrc,dot,nlit,.true.)
c ... CSRC - coeficiente da matriz simetrica
           else
             call pbicgstab(neq,nad,iax,ja,ad,au,al,ia(i_m),b,x
     .                     ,ia(i_t),ia(i_v),ia(i_r),ia(i_p),ia(i_z)
     .                     ,tol,maxit,matvec_csrc_sym,dot,nlit,.true.)
           endif
c .....................................................................
c
c ... ELLPACK
        elseif(matriz .eq. 4) then
c ... OPENMP
          if(openmpSolver) then
              call pbicgstab_omp(neq,nad,iax,ja,ad,au,al,ia(i_m),b,x
     .                          ,ia(i_t),ia(i_v),ia(i_r),ia(i_p),ia(i_z)
     .                          ,tol,maxit
     .                          ,matvec_ellpack_omp
     .                          ,dot_ompL4,nlit,.false.)
c .....................................................................
c
c ...
          else
            call pbicgstab(neq,nad,iax,ja,ad,au,al,ia(i_m),b,x
     .                    ,ia(i_t),ia(i_v),ia(i_r),ia(i_p),ia(i_z)
     .                    ,tol,maxit,matvec_ellpack,dot,nlit,.false.)
          endif
c .....................................................................
        endif
c .....................................................................
c
c ...
        i_p         = dealloc('pp      ')
        i_v         = dealloc('vi      ')
        i_t         = dealloc('t       ')  
        i_m         = dealloc('m       ')
        i_z         = dealloc('z       ')
        i_r         = dealloc('r       ')
c .....................................................................
c
c ... PGMRES
      elseif(code .eq. 3) then  
        i_r         = alloc_8('g       ',k+1,neq)
        i_z         = alloc_8('h       ',  k,neq)
        i_m         = alloc_8('m       ',  1,neq)
        i_y         = alloc_8('y       ',  1,neq)
        i_t         = alloc_8('c       ',  k,neq)
        i_v         = alloc_8('s       ',  k,neq)
        i_p         = alloc_8('ee      ',k+1,neq)
c ... precondicionador diagonal 
        precontime = get_time() - precontime  
        call aequalb(ia(i_m),ad,neq)
        precontime = get_time() - precontime 
c .....................................................................
c
c ...
        if(matriz .eq. 1) then
c ... CSRD - coeficiente da matriz nao simetrica
          if(unsym) then
            call gmres(neq,nad,iax,ja,ad,au,al
     .                ,ia(i_m),b,x,k,ia(i_r),ia(i_z)
     .                ,ia(i_y),ia(i_t),ia(i_v),ia(i_p),tol,maxit
     .                ,matvec_csrd,nlit)
c ... CSRD - coeficiente da matriz simetrica
          else
            call gmres(neq,nad,iax,ja,ad,au,al
     .                ,ia(i_m),b,x,k,ia(i_r),ia(i_z)
     .                ,ia(i_y),ia(i_t),ia(i_v),ia(i_p),tol,maxit
     .                ,matvec_csrd_sym,nlit)
          endif
        elseif(matriz .eq. 2) then
          write(*,2000)'Solver','Solver.f'
     .                ,'PGMRES nao implementado para o CSR'
          stop  
        elseif(matriz .eq. 3) then
c ... CSRC - coeficiente da matriz nao simetrica
          if(unsym) then
            call gmres(neq,nad,iax,ja,ad,au,al
     .                ,ia(i_m),b,x,k,ia(i_r),ia(i_z)
     .                ,ia(i_y),ia(i_t),ia(i_v),ia(i_p),tol,maxit
     .                ,matvec_csrc,nlit)
c ... CSRC - coeficiente da matriz simetrico
          else
            call gmres(neq,nad,iax,ja,ad,au,al
     .                ,ia(i_m),b,x,k,ia(i_r),ia(i_z)
     .                ,ia(i_y),ia(i_t),ia(i_v),ia(i_p),tol,maxit
     .                ,matvec_csrc_sym)
          endif
        endif
c .....................................................................
c
c ...
        i_p         = dealloc('ee      ')
        i_v         = dealloc('s       ')
        i_t         = dealloc('c       ')
        i_y         = dealloc('y       ')   
        i_m         = dealloc('m       ')
        i_z         = dealloc('h       ')
        i_r         = dealloc('g       ')
c .....................................................................
c
c ...
      else
        write(*,2000)'Solver','Solver.f','solver nao especificado!!!'
        stop
      endif
c .....................................................................
      if(openmpSolver .and. (.not. unsym)) then  
        i_threads_y  = dealloc('buffer_y')
      endif
      return
      end
c *********************************************************************
c
c *********************************************************************
c * GETRES : relacionando valores da solucao entre as equacoes        * 
c * e as celulas                                                      *
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c *   u    - valores por celular                                      *
c *   x    - valores por equacoes                                     *   
c * num    - renumeracaos do elementos para as equacoes               *
c * neq    - numero de equacoes                                       *
c * numel  - numero de celulas                                        *
c * flag   - true celulas -> equacoes                                 *
c *        - false equacoes -> celulas                                *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * u(atualiziado) - flad = false                                     *
c * x(atualiziado) - flad = true                                      *      
c *********************************************************************
      subroutine getRes(u,x,num,neq,numel,flag)
      implicit none
      integer neq,numel,i
      real*8 u(numel),x(neq)
      integer num(numel)
      logical flag
c ... resgatar os valores das celulas para as equacoes
      if(flag) then
        do i = 1, numel
          x(num(i)) = u(i)
        enddo
c .....................................................................
c
c ... resgatar os valores da equacoes para as celulas
      else
        do i = 1, numel
           u(i) = x(num(i))
        enddo
      endif
c .....................................................................
c
c ...
      return
      end
c *********************************************************************
