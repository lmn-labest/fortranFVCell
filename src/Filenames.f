      character*80 function name(NomeArqDados,NumArq,code)
c **********************************************************************
c *                                                                    *
c *   NAME: nomes de aquivos                                           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    NumArq       - numero do arquivo                                *
c *    code         - codigo de instrucao                              *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    NomeArqDados - nome do arquivo                                  *
c *                                                                    *
c **********************************************************************
      implicit none      
c     include 'parallel.fi'
      character*80 NomeArqDados,NomeArqGerado
      character*20 StrExtensao
      integer      iPonto,TamanhoNome,NumArq,code
c ......................................................................
      if    (code .eq. 0) then
         write( StrExtensao, '( A )' ) '.geo'
      elseif(code .eq. 1) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_trans1_'//trim(StrExtensao)//'.vtk'   
      elseif(code .eq. 2) then
         write( StrExtensao, '( A )' ) '.res'        
      elseif(code .eq. 3) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_temp_'//trim(StrExtensao)//'.vtk'
      elseif(code .eq. 4) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_flux_'//trim(StrExtensao)//'.vtk'
      elseif(code .eq. 5) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_grad_'//trim(StrExtensao)//'.vtk'
      elseif(code .eq. 6) then
        StrExtensao='_nel.dat'
      elseif(code .eq. 7) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_veloc_'//trim(StrExtensao)//'.vtk'
      elseif(code .eq. 8) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_elev_'//trim(StrExtensao)//'.vtk' 
      elseif(code .eq. 10) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_log_nl_'//trim(StrExtensao)//'.txt'
      elseif(code .eq. 11) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_pressure_'//trim(StrExtensao)//'.vtk'
      elseif(code .eq. 12) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_gradP_'//trim(StrExtensao)//'.vtk'
      elseif(code .eq. 13) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_energy_'//trim(StrExtensao)//'.vtk'
      elseif(code .eq. 14) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_masseps_'//trim(StrExtensao)//'.vtk'
      elseif(code .eq. 15) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_div_'//trim(StrExtensao)//'.vtk' 
      elseif(code .eq. 16) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_ResMass_'//trim(StrExtensao)//'.vtk'       
      elseif(code .eq. 24) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_mshape_'//trim(StrExtensao)//'.vtk'
      elseif(code .eq. 25) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_temp_node_'//trim(StrExtensao)//'.txt'
      elseif(code .eq. 27) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_temp_cell_'//trim(StrExtensao)//'.txt'
      elseif(code .eq. 28) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_grad_cell_'//trim(StrExtensao)//'.txt'
      elseif(code .eq. 29) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_nlplot_'//trim(StrExtensao)//'.vtk'
      elseif(code .eq. 30) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_save_'//trim(StrExtensao)//'.dat'
      elseif(code .eq. 31) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_gradU1_'//trim(StrExtensao)//'.vtk'
      elseif(code .eq. 32) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_gradU2_'//trim(StrExtensao)//'.vtk'
      elseif(code .eq. 33) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_pgeo_'//trim(StrExtensao)//'.vtk'  
      elseif(code .eq. 34) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_gradT_'//trim(StrExtensao)//'.vtk'  
      elseif(code .eq. 35) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_gradE_'//trim(StrExtensao)//'.vtk'
      elseif(code .eq. 36) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_vort2D_'//trim(StrExtensao)//'.vtk'  
      elseif(code .eq. 49) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_log_exc_'//trim(StrExtensao)//'.txt'
      elseif(code .eq. 50) then
        if( NumArq .eq. 1) then
          StrExtensao='pcg'
        elseif( NumArq .eq. 2) then
          StrExtensao='bicgstab'
        elseif( NumArq .eq. 3) then 
          StrExtensao='gmres'
        else
          StrExtensao='geral'
        endif
        StrExtensao='_solvit_'//trim(StrExtensao)//'.txt'
      elseif(code .eq. 51) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_log_simple_'//trim(StrExtensao)//'.txt'
      elseif(code .eq. 52) then
        write(strextensao,'( i6 )') numarq
        write(strextensao,'( a  )') adjustl(strextensao)
        strextensao='coo_'//trim(strextensao)
      elseif(code .eq. 53) then
        strextensao='_bcoo'
      elseif(code .eq. 54) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_eddyVisc_'//trim(StrExtensao)//'.vtk'
      elseif(code .eq. 55) then
        StrExtensao='_Estatic.txt'
      elseif(code .eq. 56) then
        StrExtensao='_MeanVel.vtk'
      elseif(code .eq. 57) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_pDs_'//trim(StrExtensao)//'.vtk'
      elseif(code .eq. 58) then
        write(StrExtensao,'( I6 )') NumArq
        write(StrExtensao,'( A  )') adjustl(StrExtensao)
        StrExtensao='_pyPlus_'//trim(StrExtensao)//'.vtk'
      endif    
c ................................................................
c
c ......................................................................      
      TamanhoNome = INDEX( NomeArqDados, ' '  )
      iPonto = INDEX( NomeArqdados, '.' )
      if( iPonto .EQ. 0 ) then
          NomeArqGerado = NomeArqDados(1:TamanhoNome-1) // StrExtensao
      else
          NomeArqGerado = NomeArqDados(1:iPonto-1) // StrExtensao
      endif
      name = NomeArqGerado
c ......................................................................      
      return
      end
