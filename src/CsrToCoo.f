c********************************************************************* 
c* CSRTOCOO : conveter do formato CSR para COO                       * 
c*-------------------------------------------------------------------* 
c* Parametros de entrada:                                            * 
c*-------------------------------------------------------------------* 
c* lin -> indefinido                                                 * 
c* col -> indefinido                                                 * 
c* val -> indefinido                                                 * 
c* ia  -> vetor CSR                                                  * 
c* ja  -> vetor CSR                                                  * 
c* au  -> matrix de coeficientes                                     * 
c* ad  -> matrix de coeficientes                                     * 
c* al  -> matrix de coeficientes                                     * 
c* neq -> numero de equacoes                                         * 
c* bin -> matriz binaria                                             * 
c* nad -> numero de termos nao nulos                                 * 
c*-------------------------------------------------------------------* 
c* Parametros de saida:                                              * 
c*-------------------------------------------------------------------* 
c* lin -> numero da linha                                            * 
c* col -> numero da coluna                                           * 
c* val -> valor                                                      * 
c*-------------------------------------------------------------------* 
c* OBS:                                                              * 
c*-------------------------------------------------------------------* 
c********************************************************************* 
      subroutine csrToCoo(lin  ,col,val
     .                   ,ia   ,ja
     .                   ,al   ,ad ,au
     .                   ,neq  ,nad
     .                   ,unsym,bin)
      implicit none
      integer lin(*),col(*)
      real*8  val(*)
      integer ia(*),ja(*)
      real*8 al(*),ad(*),au(*)   
      integer neq ,nad
      logical unsym,bin
      integer nl,nc,kk,n,ipoint
c
      kk = 0
      do nl = 1, neq
        kk     = kk + 1
        lin(kk) = nl
        col(kk) = nl
        if(bin) then
          val(kk) = 1.0
        else
          val(kk) = ad(nl)
        endif
        n      = ia(nl+1) - ia(nl)
        ipoint = ia(nl)
        do nc = ia(nl), ia(nl+1) - 1
          kk = kk + 1
          lin(kk)  = nl
          col(kk)  = ja(nc)
          if(bin) then
            val(kk) = 1.0
          else
            if( col(kk) .lt. nl) then
              val(kk) = al(nc)
            endif
            if( col(kk) .gt. nl .and. unsym) then
              val(kk) = au(nc)
            endif
          endif
        enddo
      enddo
c
      return
      end
c *********************************************************************
c
c *********************************************************************
c* WRITECOO : escrita do sistem linear no formato COO                * 
c*-------------------------------------------------------------------* 
c* Parametros de entrada:                                            * 
c*-------------------------------------------------------------------* 
c* lin     -> numero da linha COO                                    * 
c* col     -> numero da coluna COO                                   * 
c* val     -> valor COO                                              * 
c* f       -> vetor de forcas                                        * 
c* neq     -> numero de equacoes                                     * 
c* nnz     -> numero totoal de elementos nao nulos                   * 
c* nad     -> numero de termos nao nulos                             * 
c* prename ->                                                        * 
c* fileout -> numero do arquivo de saida                             * 
c* flag    -> coo no formato binario                                 * 
c*-------------------------------------------------------------------* 
c* Parametros de saida:                                              * 
c*-------------------------------------------------------------------* 
c*-------------------------------------------------------------------* 
c* OBS:                                                              * 
c*-------------------------------------------------------------------* 
c********************************************************************* 
      subroutine writeCoo(lin,col,val,f,neq,nnz,prename,fileOut,flag)
      implicit none
      integer lin(*),col(*),nnz,i,neq
      real*8  val(*),f(*)
      integer fileOut
      logical flag
      character*80 prename
      character*1024 str,nameOut
c
      nameOut = trim(prename) //'.mtx'
c ... matriz de coeficientes
      open(fileOut , action='write', file=nameOut)
      str = '%%MatrixMarket matrix coordinate real symmetric'
      write(fileOut,'(a)')trim(str)
      write(fileOut,'(3i9)')neq,neq,nnz
      do i = 1, nnz
        write(fileOut,'(2i9,e24.16)')lin(i),col(i),val(i)
      enddo
      close(fileOut)
c .....................................................................
c
c ...
      if( flag ) return
c .....................................................................
c
c ... vetor de forca
      nameOut = trim(prename) //'_b.mtx'
      open(fileOut , action='write', file=nameOut)
      str ='%%MatrixMarket matrix array real general'
      write(fileOut,'(a)')trim(str)
      write(fileOut,'(2i9)')neq,1
      do i = 1, neq
        write(fileOut,'(e24.16)') f(i)
      enddo
      close(fileOut)
c .....................................................................
      return
      end
c *********************************************************************
