c *********************************************************************
c * ENTHALPYFORTEMP : conversao entre temperatura e entalpia sensivel *
c * ------------------------------------------------------------------*
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * hs     - entalpia sensivel                                        *
c * t      - temperatura                                              *
c * e      - propriendade do material                                 *
c * w      - velocidades                                              *
c * numel  - numero de elementos                                      *
c * ndm    - numero de dimensoes                                      *
c * conv   - false  temperatura -> entalpia                            *
c *          true  entalpia    -> entalpia                            *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * t ou hs atualizado                                                *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine EnthalpyForTemp(hS,t,e,w,numel,ndm,conv,iws) 
      implicit none
      include 'idealGas.fi'
      real*8 hs(*),t(*),e(10,*),tRef,cP,td,w(ndm,*),tKelvin
      integer numel,ndm,i,iws
      logical conv
c     
      td   = 0.d0
      cP   = e(4,1)
      tRef = e(5,1)
      if(iws .eq. 1) td = tConv 
c ... Enthapy fot Temp
      if(conv) then
        do i = 1, numel
          t(i) = tRef 
     .         + hS(i)/cP
        enddo
c .....................................................................
c
c ... Temp for Enthapy      
      else
        do i = 1, numel
          hS(i) = cP*(t(i) -tRef) 
        enddo
      endif
      return
      end
c ***********************************************************************
c
c *********************************************************************
c * IDEALGAS: lei dos gases ideiais                                   *
c * ------------------------------------------------------------------*
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * ro     - massa especifica                                         *
c * temp   - temperatura                                              *
c * numel  - numero de celulas                                        *
c * iws    - 0 - graus kelvin                                         *
c *          1 - graus Celsius                                        *
c * iws1   - coluna do vetor massa especifica a ser atualizada        *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * ro(iws1,*) - atualizado                                           *
c * ----------------------------------------------------------------- *
c *********************************************************************      
      subroutine idealGas(ro,temp,numel,iws,iws1)
      implicit none
      include 'simple.fi'
      include 'idealGas.fi'
      real*8 ro(3,*),temp(*),t,m,a,P,fIdealGas,tKelvin
      integer numel,i,iws,iws1
c ...
      t = 0.0d0
      if( iws .eq. 1) t = tConv
      m  = dO2*mMoleculaO2 + dN2*mMoleculaN2 + dAr*mMoleculaAr 
c .....................................................................
c
c ...      
      P = (cPc0/cPc)*Pa
      do i = 1, numel
        tKelvin  =  t+temp(i)
        a         = fIdealGas(tKelvin,P,m,R)
        ro(iws1,i) = (1.d0-underRo)*ro(iws1,i) + underRo*a
      enddo
c .....................................................................
c
c ...      
      return
      end
c ********************************************************************
c 
c ********************************************************************
      subroutine water(ro,temp,numel,typeT)
      implicit none
      real*8 ro(3,*),temp(*)
      integer numel,typeT,i
c ...
      do i = 1, numel
        ro(3,i) = 1.d0 - 1.95449d-5 * dabs(temp(i) - 3.98d0)**(1.68d0)
      enddo
c .....................................................................
      return
      end
c **********************************************************************
c
c **********************************************************************
c     real*8 function cp(t,iws)
c     implicit none
c     real*8 Pa,mMoleculaO2,mMoleculaN2,mMoleculaAr
c     real*8 dN2,dO2,dAr,R
c     parameter (Pa          = 1.01325d5)!Painteger iws
c     parameter (R           = 8.31431d0)! kJ,mol,kelvin
c     parameter (mMoleculaO2 = 32.0d-3 )! kg
c     parameter (mMoleculaN2 = 28.0d-3 )! kg
c     parameter (mMoleculaAr = 39.95d-3)! kg
c     parameter (dO2         = 0.21d0) 
c     parameter (dN2         = 0.78d0)
c     parameter (dAr         = 0.01d0)
c     real*8 t,ti,m
c     integer iws
c ...   
c     ti = 0.d0
c     if( iws .eq. 1) ti = 273.15d0
c     m  = dO2*mMoleculaO2 + dN2*mMoleculaN2 + dAr*mMoleculaAr
c ......................................................................
c
c ...              
c     cp = m/(R*(t+ti))
c ......................................................................
c     return
c     end
c **********************************************************************
c
c **********************************************************************
c * IDEALGAS: calculo da massa especifica de um gs ideal               *
c -------------------------------------------------------------------- *
c * Parametro de entrada :                                             *
c -------------------------------------------------------------------- *
c * tKelvin  - temperatira do gas em kelvin                            *
c * pressure - pressão do gas                                          *
c * mMolar   - massa molar do gas                                      *
c * R        - constante dos gases ideais                              *
c -------------------------------------------------------------------- *
c * Parametro de saida:                                                *
c -------------------------------------------------------------------- *
c * fdealGas - massa especifica do gas                                 *
c **********************************************************************      
      real*8 function fIdealGas(tKelvin,pressure,mMolar,R)
      implicit none
      real*8  tKelvin,pressure,mMolar,R
c ...             
      fidealGas = (pressure*mMolar)/(R*tKelvin)
c ......................................................................
      return
      end
c **********************************************************************   