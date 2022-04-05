//LIBNP: Define as funcoes usadas pelo método de Newton padrão

#ifndef __LIBNP__
#define __LIBNP__

//troca as linhas do sistema linear
void troca_linhas (double** m_A, double* v_B, t_i_double* v_X, int i_1, int i_2, int n);

//funcao que resolve o sistema linear dado
void resolve_sistema_linear (double** m_A, t_i_double* v_X, double* v_B, int n);

//soma o vetor x1 e o delta no vetor x2
void soma_x_delta_pro_x (t_i_double* v_delta, double* v_X, int n);

//aloca um vetor e uma matriz de valores de funcoes de tamanho n
void alloca_v_m_funcao_it (double** v_f_iteracao, double*** m_f_iteracao, int n_vars);

//calcula os valores das funcoes derivadas
void calcula_valores_deriv (void** v_deriv, void*** m_deriv, double* v_f_iteracao, double** m_f_iteracao, int n_vars, char** v_vars, double* v_X);

//dá free no vetor e matriz de funcoes derivadas
void free_v_m_derivs (void**  v_deriv, void*** m_deriv, int n);

//dá free no vetor e matriz de valores de funcoes derivadas
void free_v_m_funcao_it (double*  v_f_iteracao, double** m_f_iteracao, int n);

//dá free nos vetores X e delta
void free_v_delta_X (t_i_double* v_delta, double* v_X_i1);

//faz o newton padrao
double* newton_padrao (t_entrada* entrada, int* num_it, t_tempos* tempos);

#endif