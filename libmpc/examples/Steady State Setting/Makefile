 #                                                                        
 #  Criador: Zoé Magalhães (zr.magal@gmail.com)
 #  Mestrando do PPMEC-Unb, matrícula 170172767
 #  Disciplina: Controle Preditivo 01/2018
 #  Professor: André Murilo
 #  
 #  Este software tem como objetivo disponibilizar
 #  uma biblioteca  para a implementação de controladores preditivos.
 #   
 #  MAKEFILE do codigo de exemplo de utilização de biblioteca que mostra
 #  como criar um descritor 
 # 
 #  Histórico de modificações
 #  --------------------------------------------------------------------------
 #  28/04/2018 - Zoé Magalhães
 #  - Criação do header contendo o descritor de um controlador mpc. 
 #    Os membros (todos publicos) deste descritor  são as matrizes
 #    do modelo em espaço de estados discrto, as matrizes da função 
 #    custo sem otimização, e os tamanhos do horizonte de predição,
 #    do vetor de estados e do vetor de controle. O seu único método
 #    é o construtor que recebe as matrizes do modelo, as ponderações
 #    da função custo e com base nestes valores calcula os valores
 #    dos demais membros.
 #
 #  --------------------------------------------------------------------------
 #/

LIBMPC_PATH = ../../

HEADER_DIR = $(LIBMPC_PATH)/headers
EIGEN_DIR = $(LIBMPC_PATH)/extern/eigen
CC=g++
INCLUDES=-I$(HEADER_DIR)  -I$(EIGEN_DIR)
CFLAGS = -Wall -g
ODIR=./

SOURCES := ss_setting.c
LIBS := $(LIBMPC_PATH)/obj/libmpc.o
all:
		mkdir -p $(ODIR)
		$(CC) $(CFLAGS) $(INCLUDES) -o $(ODIR)/ss_setting $(SOURCES) $(LIBS)


clean:
	 rm -r $(ODIR)/ss_setting

