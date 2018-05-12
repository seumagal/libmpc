#                                                                        
 #  Criador: Zoé Magalhães (zr.magal@gmail.com)
 #  Mestrando do PPMEC-Unb, matrícula 170172767
 #  Disciplina: Controle Preditivo 01/2018
 #  Professor: André Murilo
 #  
 #  Este software tem como objetivo disponibilizar
 #  uma biblioteca  para a implementação de controladores preditivos.
 #   
 #  MAKEFILE para compilação da biblioteca libmpc
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

HEADER_DIR = ./headers
EIGEN_DIR = ./extern/eigen
CC=g++
INCLUDES=-I$(HEADER_DIR) -I$(EIGEN_DIR)
CFLAGS = -Wall -c -g
ODIR=obj

SOURCES := sources/mpc.c
all:
		mkdir -p $(ODIR)
		$(CC) $(CFLAGS) $(INCLUDES) -o $(ODIR)/libmpc.o $(SOURCES)

clean:
	 rm -r obj

