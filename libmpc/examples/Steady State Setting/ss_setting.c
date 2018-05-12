/*                                                                        
 *  Criador: Zoé Magalhães (zr.magal@gmail.com)
 *  Mestrando do PPMEC-Unb, matrícula 170172767
 *  Disciplina: Controle Preditivo 01/2018
 *  Professor: André Murilo
 *  
 *  Este software tem como objetivo disponibilizar
 *  uma biblioteca  para a implementação de controladores preditivos.
 *   
 *  Este arquivo apresenta o código fonte de um examplo de utilização
 *  da biblioteca libmpc, mostrando como:
 *  - Iniciar um descritor de controle preditivo
 *  - Configurar o estado estacionario desejado
 *
 *
 *  Histórico de modificações
 *  --------------------------------------------------------------------------
 *  28/04/2018 - Zoé Magalhães
 *  - Criação do header contendo o descritor de um controlador mpc. 
 *    Os membros (todos publicos) deste descritor  são as matrizes
 *    do modelo em espaço de estados discrto, as matrizes da função 
 *    custo sem otimização, e os tamanhos do horizonte de predição,
 *    do vetor de estados e do vetor de controle. O seu único método
 *    é o construtor que recebe as matrizes do modelo, as ponderações
 *    da função custo e com base nestes valores calcula os valores
 *    dos demais membros.
 *  --------------------------------------------------------------------------
 *  28/04/2018 - Zoé Magalhães
 *  - Adicionadad a impressão dos valores das matrizes da equação de
 *    controle.
 *  --------------------------------------------------------------------------
 *  08/05/2018 - Zoé Magalhães
 *  - Adicionadad ao example a configuracao do estado estacionário
 *  --------------------------------------------------------------------------
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpc.h>
#include <iostream>

#define NUM_STATES      ( 4 ) // Quantidade de estados
#define NUM_COMMANDS    ( 1 ) // Quantidade de sinais de controole
#define HORIZONT_LENGTH ( 100 ) // Quantidade de amostras no horizonte de predição 

using Eigen::VectorXd;
int main()
{
    // Inícia as matrizes que são definidas pela aplicação e não pela lib:
    // - Matrizes que definem o modelo em espaço de estado discreto
    // - Matrizes de ponderação da função custo
    MatrixXd A(NUM_STATES,NUM_STATES);  
    MatrixXd B(NUM_STATES,NUM_COMMANDS);
    MatrixXd C(NUM_STATES,NUM_STATES);
    VectorXd weight_u(1);
    VectorXd weight_y(4);
    MatrixXd Q_U(1,1);
    MatrixXd Q_Y(4,4);
    VectorXd y_d(1);
    uint32_t y_d_mask[1];
    A <<   0.999854522763908,   0.000059405368279,   0.000007229707831,   0.000133743102751,
           0.000031436545097,   0.999937700765884,  -0.000007230109596,  -0.000133750536876,
          -0.002989223149497,   0.002862272105557,   0.999559196252062,  -0.008154473807345,
          -0.000000119581212,   0.000000114501670,   0.000079982371186,   0.999999673802299;


    B <<  -0.019640201237548, 0.019641292833212, -0.402513917329487, -0.000016102506172;

    B = 0.000001*B;

    C = MatrixXd::Identity(4,4);
    weight_u << 1;
    weight_y << 1,1,1,1;   
   
    Q_U = weight_u.asDiagonal();
    Q_Y = weight_y.asDiagonal();
   

    // Inicia o descritor
    MPCDescriptor descriptor( A, B, C, Q_U, Q_Y, HORIZONT_LENGTH );
    
    // Imprime o valor dos memmbros do descritor para apresentar o resultado
    std::cout << "\r\nDescritor do MPC m:\r\n A :\r\n" << descriptor.A << 
                 "\r\n B   : \r\n" << descriptor.B   << 
                 "\r\n C_R : \r\n" << descriptor.C_R << 
                 "\r\n Q_U : \r\n" << descriptor.Q_U << 
                 "\r\n Q_Y : \r\n" << descriptor.Q_Y <<   
                 "\r\n H   : \r\n" << descriptor.H   << 
                 "\r\n F1  : \r\n" << descriptor.F1  << 
                 "\r\n F2  : \r\n" << descriptor.F2  << 
                 "\r\n F3  : \r\n" << descriptor.F3  << 
                 "\r\n K   : \r\n" << descriptor.K  << 
                 "\r\n L   : \r\n" << descriptor.G  << 
                 "\r\n G   : \r\n" << descriptor.L  << 
                 "\r\n n : " << descriptor.n << 
                 " nu  : " << descriptor.nu << 
                 " nx : " << descriptor.nx <<
                 "\n y_d : \r\n" << y_d << std::endl;


    y_d(0) = 0.001;  
    y_d_mask[0] = 1;

    descriptor.MPCSetDesiredSteadyState( y_d, y_d_mask, 1 );

}
