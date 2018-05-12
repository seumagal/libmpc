/* ****************************************************************************                                                                       
 *  Criador: Zoé Magalhães (zr.magal@gmail.com)
 *  Mestrando do PPMEC-Unb, matrícula 170172767
 *  Disciplina: Controle Preditivo 01/2018
 *  Professor: André Murilo
 *  
 *  Este software tem como objetivo disponibilizar
 *  uma biblioteca para a implementação de controladores preditivos.
 *
 *  Histórico de modificações
 *  --------------------------------------------------------------------------
 *  - 28/04/2018 - Zoé Magalhães
 *  - Criação do header contendo o descritor de um controlador mpc. 
 *    Os membros (todos publicos) deste descritor  são as matrizes
 *    do modelo em espaço de estados discrto, as matrizes da função 
 *    custo sem otimização, e os tamanhos do horizonte de predição,
 *    do vetor de estados e do vetor de controle. O seu único método
 *    é o construtor que recebe as matrizes do modelo, as ponderações
 *    da função custo e com base nestes valores calcula os valores
 *    dos demais membros.
 *  --------------------------------------------------------------------------
 *   - 30/04/2018 - Zoé Magalhães
 *   - Adicionado ao descritor as matrizes da equção de controle
 *   - Adicionado ao descritor um método para calcular o próximo
 *     sinal de controle. Este método recebe como parâmetro o 
 *     estado atual e o vetor de controle de estado estacionário
 *     desejado.
 *  --------------------------------------------------------------------------
 *   - 12/05/2018 - Zoé Magalhães
 *   - Disponibilizando métodos para calcular o valor desejado do
 *     comando de estado estacionário com base no valor desejado
 *     da saída de estado estacionário.
 *  --------------------------------------------------------------------------
 */
#ifndef __mpc__
#define __mpc__
#include <stdlib.h>
#include <stdint.h>
#include <Eigen/Dense>
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::FullPivLU;
/*
 * @file mpc.h
 * Este header define a interface da biblioteca libmpc.h, nele deve
 * estar contido a definição dos tipos e funções necessários
 * para utilizar a lib.
 */

/*
 * @brief protótipo da callback cadastrada para ler o valor 
 * do estado.
 */

/*
 * Descritor de um controlador preditivo
 */
class MPCDescriptor 
{
    public:
        MatrixXd A;    /**< matriz A do modelo discreto em espaço de estados */
        MatrixXd B;    /**< matriz B do modelo discreto em espaço de estados */
        MatrixXd C_R;  /**< matriz de observação y = C_R*x */
        MatrixXd Q_U;  /**< matriz de ponderação do vetor de comando na funcão custo*/
        MatrixXd Q_Y;  /**< matriz de ponderação do vetor de regulavel na função custo*/
        MatrixXd H;    /**< matriz Hessiano da função custo */
        MatrixXd F1;   /**< matriz F1  da função custo */
        MatrixXd F2;   /**< matriz F2  da função custo */
        MatrixXd F3;   /**< matriz F3  da função custo */
        MatrixXd K;    /**< matrix K da equação de controle, ganho no estado atual */
        MatrixXd G;    /**< matrix G da equação de controle, ganho da trajetória desejada */
        MatrixXd L;    /**< matrix L da equação de controle, ganho do comando em estado estacionário */
        VectorXd u;    /**< ùltimo vetor de comando caulcado */
        VectorXd u_s;  /**< estado desejado para regime permanente */
        uint32_t n;    /**< número de amostras do horizonte de predição */
        uint32_t nu;   /**< comprimento do vetor de controle */
        uint32_t nx;   /**< comprimento do vetor de estados */
        VectorXd u_d;  /**< vetor de comando desejado para o estado estacionário */
        VectorXd x_d;  /**< vetor de estado odesejado para o estado estacionário */ 
        
        /*
         * @brief Construtor do descritor de um controlador preditivo. Inicializa
         * todos os membros do descritor.
         *
         * @param arg_A é a matriz A do modelo discreto em espaço de estados
         * @param arg_B é a matriz B do modelo discreto em espaço de estados
         * @param arg_C_R é a matriz de observação y = C_Rx
         * @param arg_Q_U é a matriz de ponderação do vetor de comando na função custo
         * @param arg_Q_Y é a matriz de ponderação do erro no vetor regulado.
         * @param arg_n é a quantidade de amostras no horizonte de predição
         * @param arg_write_callback é a função cadastrada para leitura dos estados
         */
         MPCDescriptor( MatrixXd arg_A,
                        MatrixXd arg_B,
                        MatrixXd arg_C_R,
                        MatrixXd arg_Q_U,
                        MatrixXd arg_Q_Y,
                        uint32_t arg_n );
        
         /*
          * @breif Calcula o comando e os estados esperado para regime permanente
          * @param arg_desired_steady_y é o vetor regulado desejado em regime
          * permanente.
          * 
          * @param arg_mask é o vetor que indica os elementos do vetor regulado
          * que possuem valor de estado estacionário definido. 
          * @return retorna #TRUE se foi encontrada uma solução para o estado
          * e comando desejado em regime pernanente. Retorna #FALSE se
          * não foi encoontrada uma solução.
          */
         bool MPCSetDesiredSteadyState( VectorXd arg_desired_steady_y );
         
         /*
          * @breif Calcula o comando e os estados esperado para regime permanente,
          * quando nem todas as saídas y_d tem um estado desejado.
          *
          * @param arg_desired_steady_y é o vetor regulado desejado em regime
          * permanente.
          * @param arg_mask é o vetor que indica os elementos do vetor regulado
          * que possuem valor de estado estacionário definido.
          * @param arg_sel é um buffer com o indíce dos estados com valor
          * de estado estacionário definido.
          * @param arg_sel_size indica a quantidade de indices no buffer arg_sel. 
          *
          * @return retorna #TRUE se foi encontrada uma solução para o estado
          * e comando desejado em regime pernanente. Retorna #FALSE se
          * não foi encoontrada uma solução.
          */
         bool MPCSetDesiredSteadyState( VectorXd arg_desired_steady_y,
                                        uint32_t * arg_sel,
                                        uint32_t arg_sel_size );



         /*
         * @brief Calcula o próximo vetor de controle u 
         * @param arg_x é o estado atual
         * @param arg_tracking é a trajetória desejada na janela de predição
         * @return retorna #true em caso de sucesso.
         *         retorna #false se os vetores estiverem tamanhos incompatíveis
         *         retorna #false se o descritor não estiver devidamente construido
         */
         bool MPCWriteNextCommand( VectorXd arg_x, MatrixXd arg_tracking ) ; 
};

#endif
