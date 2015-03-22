#!/usr/bin/python
import sys

def main(argv):
  
  progArg = (['[massa especifica] [viscosidade molecular] ' 
              '[velociade caracteristica] [escala caracteristica'])
#checando os agumentos  
  nArgs = len(argv)
  if nArgs < 5:
    sys.stderr.write("Usage: %s "%argv[0])
    for arg in progArg:
      print arg 
    return 1
#  
  massaEspecifica      = float(argv[1])
  viscosidadeMolecular = float(argv[2])
  velocidade           = float(argv[3])
  escala               = float(argv[4])
#
  epsilon               = (velocidade**3)/escala
  viscosidadeCinematica = viscosidadeMolecular/massaEspecifica
#
  re                   = velocidade*escala/viscosidadeCinematica
#escala energetica-inercial (lei)
  lEI  = (0.43*escala)/6.0


#escala de kolmogorov
  neta                  =  ((viscosidadeCinematica**3)/epsilon)**(1.0/4.0) 
  tau                   =  (viscosidadeCinematica/epsilon)**(1.0/2.0) 
  vel                   =   (viscosidadeCinematica*epsilon)**(1.0/4.0) 
  reynoldsKolmogorov    =   vel*neta/viscosidadeCinematica

#escala dissipativa-inercial (lDI)
  lDI  = 60.0*neta


#
  print 'escala ener-inercial               : '     ,lEI    ,'  (m)'
  print 'escala diss-inercial               : '     ,lDI    ,'  (m)'
  print 'escala espacial de kolmogorov      : '     ,neta   ,'  (m)'
  print 'escala temporal de kolmogorov      : '     ,tau    ,'  (s)'
  print 'escala de velocidade de kolmogorov : ',vel    ,'(m/s)'
  print 'taxa de dissipacao                 : ',epsilon,'  (m**2/s**3)'
  print 'Reynolds                           : ',re
  print 'Reynolds escala de kolmogorov      : ',reynoldsKolmogorov 


if __name__ == '__main__':
  sys.exit(main(sys.argv))

