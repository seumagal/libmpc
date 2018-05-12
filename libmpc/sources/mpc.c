/*
 *  Criador: Zoé Magalhães (zr.magal@gmail.com)
 *  Mestrando do PPMEC-Unb, matrícula 170172767
 *  Disciplina: Controle Preditivo 01/2018
 *  Professor: André Murilo
 *  
 *  Este software tem como objetivo disponibilizar
 *  uma biblioteca em C para a implementação 
 *  de controladores preditivos.
 *
 *  Histórico de modificações
 *  --------------------------------------------------------------------------
 *  - 28/04/2018 - Zoé Magalhães
 *  - Criação do arquivo contendo a implementação do construtor do descritor
 *    e da função interna que fornece a matriz de selação.
 *  --------------------------------------------------------------------------
 *  - 30/04/2018 - Zoé Magalhães
 *  - Adicionado ao construtor do descritor o calculo das matrizes da
 *    equação de controle. Criado corpo vazio do método que calcula
 *    o próximo vetor de controle.
 *  --------------------------------------------------------------------------
 *  - 08/05/2018 - Zoé Magalhães
 *  - Adicionado o metodo que configura o vetor de estado e o vetor de
 *    comando desejado para o estado estacionário baseado no vetor 
 *    de saída desejado para o estado estacionário. A versão atual 
 *    resolve apenas o caso em que  a solução é único e o caso 
 *    trivial em que a Y=X e y_d=0;
 *  --------------------------------------------------------------------------
 *  - 12/05/2018 - Zoé Magalhães
 *  - Conclusão da implementação dos métodos para calcular o valor do
 *    comando de estado estacionário desejado com base na saída
 *    de estado estacionário desejada.
 *  --------------------------------------------------------------------------
 */
#include <mpc.h>
#include <iostream>

/*
 * @brief Retorna uma matrix de seleção Pi
 * @param arg_n é o comprimento do horizonte de predição
 * @param arg_sel_size é o tamanho do vetor selecionado
 * @param arg_sel é o indice do vetor selecionado
 * @param arg_PI é o endereço onde a matriz
 * será copiada a matriz PI de tamanho arg_sel_size x arg_sel_size*arg_n,
 * tal que se U = [u_1; u_2; u_3; ...; u_<arg_n>] e length(u_i) = arg_sel_size,
 * então PI*U = u_<arg_sel>.
 *
 */
static  void sel_matrix( uint32_t arg_n,
                             uint32_t arg_sel_size,
                             uint32_t arg_sel,
                             MatrixXd * arg_PI )
{
    if( arg_PI == NULL )
        return; 
   
    if( arg_sel > arg_n )
        return;

    arg_PI->setZero( arg_sel_size , arg_sel_size*arg_n );
    arg_PI->block( 0, (arg_sel-1)*arg_sel_size, arg_sel_size, arg_sel_size ) = MatrixXd::Identity( arg_sel_size, arg_sel_size);
} 


/*
 * @brief Construtor do descritor de um controlador MPC. Inicializa
 * todos os membros do descritor.
 */
MPCDescriptor::MPCDescriptor( MatrixXd arg_A,
                              MatrixXd arg_B,
                              MatrixXd arg_C_R,
                              MatrixXd arg_Q_U,
                              MatrixXd arg_Q_Y,
                              uint32_t arg_n )
{
    uint64_t i,j;
    MatrixXd PHI, PSI, APOW, PI, AUX;
   
 //VALIDA OS PARÂMETROS
    // Q_U e B multiplicam U, então ambos tem que ter o mesmo número de colunas
    if( arg_Q_U.cols() != arg_B.cols() )
        return;

    // Q_Y e A e C_R . então todos devem ter o mesmo número de colunas 
    if( arg_Q_Y.cols() != arg_A.cols() )
        return;
    if( arg_Q_Y.cols() != arg_C_R.cols() )
        return;

    // A, C_R, Q_U e Q_Y devem ser matrizes quadradas
    if( arg_A.rows() != arg_A.cols() )
        return;
    if( arg_C_R.rows() != arg_C_R.cols() )
        return;
    if( arg_Q_U.rows() != arg_Q_U.cols() )
        return;
    if( arg_Q_Y.rows() != arg_C_R.cols() )
        return;

    // B e A devem ter o mesmo número de linhas
    if( arg_B.rows() != arg_A.rows() )
        return;


//INICIA OS MEMBROS QUE NÃO PRECISAM SER CALCULADOS
    // Cadastra as matrizes recebidas
    A = arg_A;
    B = arg_B;
    C_R = arg_C_R;
    Q_U = arg_Q_U;
    Q_Y = arg_Q_Y;
    // Inicia os comprimentos dos vetores de:
    n = arg_n;           // Predição
    nu = arg_B.cols();   // controle
    nx = arg_A.cols();   // estado

//CALCULA AS MATRIZES DE OTIMIZAÇÃO    
    // Define o tamanho e inicia as matrizes
    F1.setZero( nu*n ,nx   );
    F2.setZero( nu*n ,nx*n );
    F3.setZero( nu*n ,nu   );
    H.setZero( nu*n, nu*n );

    PHI.setIdentity( nx   , nx );
    for( i = 1; i <= n ; i++ )
    {
        // Calcula as matrizes de predição
        PHI = PHI*A;
        PSI.setZero( nx , nu*n  );
        APOW.setIdentity( A.rows(), A.cols() );
        for( j = i; j>=1; j-- )
        {

            sel_matrix( n , nu, j, &PI );
            PSI = PSI + APOW*B*PI;
            if( j > 1 )
                APOW = APOW*A;
        }


        // Calculo parcial das matrizes da função custo 
        sel_matrix( n, nu, i, &PI );  
        H = H + PSI.transpose()*C_R.transpose()*Q_Y*C_R*PSI + PI.transpose()*Q_U*PI;
        F1 = F1 + PSI.transpose()*C_R.transpose()*Q_Y*C_R*PHI;
        F3 = F3 + PI.transpose()*Q_U;
        sel_matrix( n, nx, i, &PI );  
        F2 = F2 + PSI.transpose()*C_R.transpose()*Q_Y*PI;
    }

    // Concluì o cálculo das matrizes da função custo
    H  = 2*H;
    F1 = 2*F1;
    F2 = -2*F2;
    F3 = 2*F3;
    
    // Calcula as matrizes da equação de controle
    // u[k] = -Kx[k] + G*(~yd[k]) + L*ud[k]
    sel_matrix( n, nu, 1, &PI );
    AUX = PI*H.inverse();
    K = AUX*F1;
    G = -AUX*F2;
    L = -AUX*F3;
}

 
bool MPCDescriptor::MPCSetDesiredSteadyState( VectorXd arg_desired_steady_y )
{

    MatrixXd Mss;
    VectorXd aux1;
    VectorXd aux2;

    // Verifica se o vetor recebido é compatível com o sistema configurado
    if( arg_desired_steady_y.size() != nx )
        return false;

    // Verifica se o sistema pode possuir solução única
    if( nx != nu )
        return false;
    

    Mss.setZero(2*nx, nx+nu);

    Mss.block( 0,   0,  nx, nx ) = MatrixXd::Identity( nx, nx) - A;
    Mss.block( 0,  nx,  nx, nu ) = -B;
    Mss.block( nx,  0,  nx, nx ) = C_R;
   
    Eigen::FullPivLU<MatrixXd> lu(Mss);
    if( !lu.isInvertible() )
        return false;
    
    aux1.setZero(2*nx);
    aux1.segment(nx,nx) = arg_desired_steady_y; 
    aux2.setZero(nx+nu);

    u_d.setZero(nu);
    x_d.setZero(nx);

    aux2 = Mss.inverse()*aux1;
    x_d = aux2.segment(0,nx);
    u_d = aux2.segment(nx,nu);

    return true;

}

bool MPCDescriptor::MPCSetDesiredSteadyState( VectorXd arg_desired_steady_y,
                                              uint32_t *arg_sel, uint32_t arg_sel_size )
{

    MatrixXd Mss;
    MatrixXd Amask(nu,nu);
    MatrixXd Bmask(nu,nu);
    MatrixXd Cmask;
    VectorXd desired_y;
    VectorXd aux2;
    VectorXd aux1;
    uint32_t i,j;


    // Verifica se os vetores recebidos tem o mesmo tamanho
    if( arg_desired_steady_y.size() != arg_sel_size )
        return false;

    //Verifica se os vetores recebidos tem o mesmo tamanho do vetor de controle
    //para que se tenha um sistema quadrado
    if( arg_desired_steady_y.size() != nu )
        return false;


    for( i=0; i< nu; i++ )
    {
        for( j=0; j < nu; j++ )
        {

            if( arg_sel[i] >= nu )
                return false;

            Amask(i,j) =   A(arg_sel[i],arg_sel[j]); 
            Cmask(i,j) = C_R(arg_sel[i],arg_sel[j]);
        }

        Bmask.row(i) = B.row(arg_sel[i]);

    }

    Mss.setZero(2*nu, 2*nu);

    Mss.block(  0,  0, nu, nu ) = MatrixXd::Identity( nu, nu) - Amask;
    Mss.block(  0, nu, nu, nu ) = -Bmask;

    Eigen::FullPivLU<MatrixXd> lu(Mss);
    if( !lu.isInvertible() )
        return false;
    
    aux1.setZero(2*nu);
    aux1.segment(nu,nu) = arg_desired_steady_y; 
    aux2.setZero(2*nu);

    u_d.setZero(nu);

    aux2 = Mss.inverse()*aux1;

    u_d = aux2.segment(nu,nu);

    return true;

}

