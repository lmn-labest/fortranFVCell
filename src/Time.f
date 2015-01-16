c *********************************************************************     
      real*8 function get_time()
      implicit none
      include 'openmp.fi'
      real*8 time
      time = 0.d0
c$    time = omp_get_wtime()
      get_time = time
      return
      end
c *********************************************************************
c
c *********************************************************************
c * WRITE_LOG : Arquivo de log de execucao                            *
c *********************************************************************
      subroutine write_log(ntime,prename,neqU1,neqU2,neqPc,neqE
     .                    ,nadU1,nadU2,nadPc,nadE
     .                    ,bandaU1,bandaU2,bandaPc,bandaE
     .                    ,cfl,re,prandlt,grashof,massa0,massa,fluxoM
     .                    ,vol)
      Use Malloc
      implicit none
      include 'time.fi'
      include 'string.fi'
      character*80 prename,fileout,name
      integer ntime
      integer neqU1,neqU2,neqPc,neqE
      integer nadU1,nadU2,nadPc,nadE
      integer bandaU1,bandaU2,bandaPc,bandaE
      real*8 cfl,re,prandlt,grashof,vol,massa0,massa,fluxoM
      real *8 use_work_vector
      logical unsym
c ...
      fileout = name(prename,0,49) 
      open(ntime, file=fileout ,action= 'write')
c ... Tempo
      write(ntime,'(a)') '********************************************'
      write(ntime,'(a)') 'Tempos:'
      write(ntime,'(a)') '********************************************'
      write(ntime,10) 'READ      : ',readTime
      write(ntime,10) 'ADJANCENCY: ',vizTime
      write(ntime,10) 'REORD     : ',readTime
      write(ntime,10) 'DATASTRUCT: ',dsTime
      write(ntime,10) 'MATVEC    : ',matvecTime
      write(ntime,10) 'PRECON    : ',preconTime
      write(ntime,10) 'DOT       : ',dotTime
      write(ntime,10) 'VECTIME   : ',vecTime
      write(ntime,10) 'TOTAL     : ',totalTime
      write(ntime,'(a)') '********************************************'
      write(ntime,'(a)') 'SIMPLE:'
      write(ntime,'(a)') '********************************************'
      write(ntime,10) 'SIMPLE    : ',simpleTime
      write(ntime,10) 'PFORMU1   : ',elmU1Time
      write(ntime,10) 'PFORMU2   : ',elmU2Time
      write(ntime,10) 'PFORMPC   : ',elmPcTime
      write(ntime,10) 'PFORME    : ',elmETime
      write(ntime,10) 'SOLVERU1  : ',solvU1Time
      write(ntime,10) 'SOLVERU2  : ',solvU2Time
      write(ntime,10) 'SOLVERPC  : ',solvPCTime
      write(ntime,10) 'SOLVERE   : ',solvETime
      write(ntime,10) 'POSVEL    : ',posVelocityTime
      write(ntime,10) 'UPSIMPLE  : ',simpleUpdateTime
      write(ntime,10) 'GFORMP    : ',grSimpleTime
c .....................................................................
c
c ... Informacao do sistem linear
      write(ntime,'(a)') '********************************************'
      write(ntime,'(a)') 'Sistema de equacoes:'
      write(ntime,'(a)') '********************************************'
      write(ntime,11) 'NEQU1     : ',neqU1
      write(ntime,11) 'NEQU2     : ',neqU2 
      write(ntime,11) 'NEQPc     : ',neqPc
      write(ntime,11) 'NEQE      : ',neqE       
      write(ntime,11) 'NADU1     : ',nadU1
      write(ntime,11) 'NADU2     : ',nadU2
      write(ntime,11) 'NADPc     : ',nadPc
      write(ntime,11) 'NADE      : ',nadE 
c      write(ntime,11) 'BANDAU1   : ',bandaU1
c      write(ntime,11) 'BANDAU2   : ',bandaU2
c      write(ntime,11) 'BANDAPc   : ',bandaPc
c      write(ntime,11) 'BANDAE    : ',bandaE 
c .....................................................................
c
c ... Informacao do sistem linear
      write(ntime,'(a)') '********************************************'
      write(ntime,'(a)') 'Sistema de equacoes:'
      write(ntime,'(a)') '********************************************'
      write(ntime,10) 'CFL       : ',cfl  
      write(ntime,10) 'Reynalds  : ',re
      write(ntime,10) 'Prandlt   : ',prandlt
      write(ntime,9 ) 'Grashof   : ',grashof
      write(ntime,9 ) 'Rayleigh  : ',grashof*prandlt
      write(ntime,9 ) 'MassaI    : ',Massa0            
      write(ntime,9 ) 'MassaF    : ',Massa              
      write(ntime,9 ) 'FluxoMassa: ',FluxoM               
      write(ntime,10) 'Volume(m3): ',vol        
c .....................................................................
c
c ... Informacao do sistem linear
      write(ntime,'(a)') '********************************************'
      write(ntime,'(a)') 'Uso da memoria: '    
      write(ntime,'(a)') '********************************************'
      write(ntime,10) 'Vetor de trabalho MB: '
     .                ,use_work_vector('MB      ')
c .....................................................................      
      return
    9 format(a,es15.5)
   10 format(a,f16.6)
   11 format(a,i10)
   12 format(a,10a)
      end
c .....................................................................
c *********************************************************************       