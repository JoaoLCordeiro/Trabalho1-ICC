//LIBGERAL: Define as funcoes usadas por mais de um método e as structs usadas

#ifndef __LIBGERAL__
#define __LIBGERAL__

//função que lê a entrada
int 	le_entrada (t_entrada* entrada);

//função que encontra o índice do maior número da coluna para baixo da linha = coluna
int 	encontra_pivo (double** matriz, int coluna, int n);

//retrossubstituicao em um vetor de double com índice
void 	retrossubs_v_i_double (double** m_A, t_i_double* v_X, double* v_B, int n);

//o delta pode ser desordenado na hora do pivoteamento, essa funcao o ordena de volta
void 	reordena_v_i_double (t_i_double* v_i_double, int n);

//calcula a norma de um vetor de double
double 	norma (double*  v_valores, int n);

//calcula a norma de um vetor de double com indice
double 	norma_i (t_i_double*  v_valores, int n);

//calcula as funcoes derivadas da funcao "funcao"
void 	calcula_derivadas (void* funcao, void** v_deriv, void*** m_deriv, char** v_vars, int n);

//aloca um vetor e uma matriz de funcoes derivadas de tamanho n
void 	alloca_v_m_derivs (void*** v_deriv, void**** m_deriv, int n);

//aloca um vetor de valores de resultados relativos a cada iteração
void 	alloca_v_double (double** v_double, int n);

//aloca os vetores: delta, xi1 e xi2
void 	alloca_v_delta_X (t_i_double** v_delta, double** v_X, int  n_vars);

//inicia o vetor xi1 com os valroes iniciais e coloca índice nos valores de delta
void 	inicia_X_delta (double* v_X_i1, t_i_double* v_delta, t_entrada* entrada);

//funcao que dá free na entrada
void 	free_entrada (t_entrada* entrada);

//retorna o maior dos valores dos parâmetros
int 	retorna_max (int a, int b, int c);

//imprime os resultados do programa
void 	imprime_resultados (t_entrada* entrada, double* v_res_np, double* v_res_nm, double* v_res_ni, int num_it_np, int num_it_nm, int num_it_ni, t_tempos tempos_np, t_tempos tempos_nm, t_tempos tempos_ni);

#endif