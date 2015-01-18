#!/bin/sh
# $URL: svn+ssh://henrique@146.164.26.19/mnt/shared/repositorio/mefpar/branches/geoacustica/mef/Makefile $
# $Id: Makefile 2 2009-04-03 17:46:50Z helano $

#make
#------------------Variaveis Globais-----------------------
PATH_MOD="src/"
PATH_INCLUDE="include"    
PRENAME=mvf
FC=gfortran
OPENMP=yes
DEBUG=yes
#------------------gerando o nome do excutavel-------------
ifeq ($(FC),ifort)
  COMPILER_NAME=intel
endif

ifeq ($(FC),gfortran)
  COMPILER_NAME=gnu
endif

NAME+=$(PRENAME)_$(COMPILER_NAME)
#-------------------Fontes---------------------------------
src_mef = src/Adjacency.f\
src/Assbly.f\
src/Celllib.f\
src/cellibAdvec.f\
src/CellibDif.f\
src/CellibPressure.f\
src/CellibSimple.f\
src/CellibTrans.f\
src/Correct.f\
src/Csr.f\
src/Datastruct.f\
src/Filenames.f\
src/Graph.f\
src/IdealGas.f\
src/Iterativos.f\
src/Iterativos_omp.f\
src/Main.f\
src/Matvec.f\
src/Matrix_partition.f\
src/Matvec_omp.f\
src/Numeq.f\
src/Pdata.f\
src/Pform.f\
src/Reord.f\
src/Rdata_mvf.f\
src/Shape_matrix.f\
src/Simple.f\
src/Solver.f\
src/Time.f\
src/Transporte.f\
src/Viz.f\
src/Vtk.f\
src/Write_vtk.f


mod_mef=src/Malloc.f
#-------------------Flags necessarios--------------------------------
NFLAGS=-I$(PATH_INCLUDE)    
LDFLAGS=
#--------------------compiladores------------------------------------
# intel ifort
ifeq ($(FC),ifort)
  LDFLAGS += 
  OFLAGS  += -module $(PATH_MOD) -W1 
  ifeq ($(OPENMP),yes)
    OFLAGS  += -openmp
  endif
endif
# gnu gcc
ifeq ($(FC),gfortran)
  LDFLAGS += 
  OFLAGS  += -J$(PATH_MOD) -Wall
  ifeq ($(OPENMP),yes)
    OFLAGS  += -fopenmp
  endif
endif

#--------------------------------------------------------------------
#---------------------------Debug------------------------------------
ifeq ($(DEBUG),yes)
  OFLAGS += -g 
else
  OFLAGS += -O2  
endif
#--------------------------------------------------------------------
FFLAGS= $(NFLAGS) $(OFLAGS) 

.SUFFIXES: 
.SUFFIXES: .for .f .h .fi .o
mod_src  = $(mod_mef:%.f=%.o)
objs_src= $(mod_src) $(src_mef:%.f=%.o)

build:	$(objs_src) $(mod_src)	
	mkdir -p bin
	$(FC) $(FFLAGS) $(objs_src) -o bin/$(NAME) $(LDFLAGS)   

tags:
	ctags -R src/*.f include/*.fi

.PHONY: cleantags
cleantags:
	@rm -fv tags
	
.PHONY: clean
clean:  
	@rm -fv src/*.o
	@rm -fv src/*.mod
	@rm -fv bin/$(NAME)

.PHONY: cleanall
cleanall:  
	@rm -fv tags
	@rm -fv src/*.o
	@rm -fv src/*.o
	@rm -fv src/*.mod
	@rm -fv src/*.f90
	@rm -fv bin/$(NAME)

.PHONY: cleanmod
cleanmod:  
	@rm -fv ../src/*.mod

.PHONY: help
help:
	@echo "Autor :$(AUTHOR)                              "
	@echo "Makefile para prepar para sitemas linux.      "
	@echo -e "\E[7;32mOpcoes:\E[1;0m                      "
	@echo "build         - compila o prepar              "
	@echo "build_modules - gera os modulos               "
	@echo "tags          - gera os tags                  "
	@echo "cleantags     - limpa os tags                 "
	@echo "clean         - limpa os obj, bin e mod       "
	@echo "cleaall       - limpa tudo obj,bin,mod e tags "
	@echo "cleanmod      - limpa os modulos              "

# DO NOT DELETE

# DO NOT DELETE
