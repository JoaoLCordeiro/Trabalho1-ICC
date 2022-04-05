//LIBNM: Define as funcoes usadas pelo método de Newton modificado

#ifndef __LIBNM__
#define __LIBNM__

//aloca espaço para o vetor de derivadas primeiras e as duas matrizes LU
void alloca_v_mLU (double** v_fun_it, double*** m_fun_U, double*** m_fun_L, int n_vars);

//calcula o valor das derivadas primeiras (gradiente)
void calcula_val_deriv_v (void** v_deriv, double* v_fun_it, int n_vars, char** v_vars, double* v_X);

//calcula o valor das derivadas segundas (hessiana)
void calcula_val_deriv_LU (void*** m_deriv, double** m_fun_U, int n_vars, char** v_vars, double* v_X);

//realiza uma troca de linhas nas matrizes LU e guarda info das trocas em delta
void troca_linhas_LU (double** m_U, double** m_L, t_i_double* v_delta, int i_1, int i_2, int n);

//faz fatoração LU com pivoteameto
void LU_pivot (double** m_U, double** m_L, t_i_double* v_delta, int n);

//após encontrar um vetor gradiente, essa funcao o ordena de acordo com as trocas que ocorreram na fatoracao LU
void troca_v_fun_it (double* v_fun_it, t_i_double* v_delta, int n);

//faz retorssubstituicao com a matriz L, começando de cima
void retrossubs_v_L (double** m_A, double* v_X, double* v_B, int n);

//soma os valores de delta em v_X de acordo com seus índices
void passa_delta_pra_X (t_i_double* v_delta, double* v_X, int n);

//dá free na matriz LU
void free_v_m_LU (double* v_fun_it, double** m_fun_U, double** m_fun_L, int n);

//faz o newton padrão
double* newton_modificado (t_entrada* entrada, int* num_it, t_tempos* tempos);

#endif