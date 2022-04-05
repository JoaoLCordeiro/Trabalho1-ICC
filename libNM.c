//LIBNM: Define as funcoes usadas pelo método de Newton modificado

#include <stdio.h>
#include <stdlib.h>
#include <matheval.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "libDefine.h"
#include "libGeral.h"
#include "libNP.h"
#include "libNM.h"
#include "libNI.h"

#define MAX 1000

//aloca espaço para o vetor de derivadas primeiras e as duas matrizes LU
void alloca_v_mLU (double** v_fun_it, double*** m_fun_U, double*** m_fun_L, int n_vars){
	*v_fun_it		= (double *)  calloc (n_vars,sizeof(double));	
	*m_fun_U		= (double **) calloc (n_vars,sizeof(double*));
	*m_fun_L		= (double **) calloc (n_vars,sizeof(double*));	
	for (int i = 0 ; i < n_vars ; i++){
		(*m_fun_U)[i]		= (double *)  calloc (n_vars,sizeof(double));
		(*m_fun_L)[i]		= (double *)  calloc (n_vars,sizeof(double));
	}
}

//calcula o valor das derivadas primeiras (gradiente)
void calcula_val_deriv_v (void** v_deriv, double* v_fun_it, int n_vars, char** v_vars, double* v_X){
	for (int i = 0 ; i < n_vars ; i++){
		v_fun_it[i]	= -evaluator_evaluate(v_deriv[i], n_vars, v_vars, v_X);
	}
}

//calcula o valor das derivadas segundas (hessiana)
void calcula_val_deriv_LU (void*** m_deriv, double** m_fun_U, int n_vars, char** v_vars, double* v_X){
	for (int i = 0 ; i < n_vars ; i++)
		for (int j = 0 ; j < n_vars ; j++)
			m_fun_U[i][j] = evaluator_evaluate(m_deriv[i][j], n_vars, v_vars, v_X);
}

//realiza uma troca de linhas nas matrizes LU e guarda info das trocas em delta
void troca_linhas_LU (double** m_U, double** m_L, t_i_double* v_delta, int i_1, int i_2, int n){
	double* auxp;
	int aux;

	auxp 	 = m_U[i_1];				//troca as linhas da matriz U
	m_U[i_1] = m_U[i_2];
	m_U[i_2] = auxp;

	auxp 	 = m_L[i_1];				//troca as linhas da matriz L
	m_L[i_1] = m_L[i_2];
	m_L[i_2] = auxp;

	aux 	 		= v_delta[i_1].i;	//guarda as trocas no delta
	v_delta[i_1].i 	= v_delta[i_2].i;
	v_delta[i_2].i 	= aux;
}

//faz fatoração LU com pivoteameto
void LU_pivot (double** m_U, double** m_L, t_i_double* v_delta, int n){
	for (int i = 0 ; i < n ; i++){						//i = linha encontrando o pivo
		int iPivo	= encontra_pivo(m_U, i, n);

		if (i != iPivo)
			troca_linhas_LU (m_U, m_L, v_delta, i, iPivo, n);

		for (int j = i+1 ; j < n ; j++){				//j = linha subtraindo a linha do pivo	
			m_L[i][j] = 0.0;							//parte superior da matriz L é 0

			double m 	= m_U[j][i]/m_U[i][i];
			m_U[j][i]	= 0.0;

			for (int k = i+1 ; k < n ; k++)				//k = coluna fazendo a subtracao por elemento
				m_U[j][k]	-= m_U[i][k]*m;
			
			m_L[j][i] = m;								//elementos abaixo da diagonal principal
		}
		m_L[i][i] = 1.0;								//diagonal principal da matriz L
	}
	//agora as matrizes L e U estão prontas para serem usadas, com as mudanças salvas no delta
}

//após encontrar um vetor gradiente, essa funcao o ordena de acordo com as trocas que ocorreram na fatoracao LU
void troca_v_fun_it (double* v_fun_it, t_i_double* v_delta, int n){
	for (int j = 0 ; j < n ; j++){
		if (j < v_delta[j].i){					//quando o indice atual for menor do indice do delta
			for (int k = j+1 ; k < n ; k++){	//procura o lugar que esta com o indice correto e troca
				if (j == v_delta[k].i){
					double aux 	= v_fun_it[j];
					v_fun_it[j]	= v_fun_it[k];	//quando o indice atual é maior, a troca já aconteceu
					v_fun_it[k]	= aux;
					break;
				}
			}
		}
	}
}

//faz retorssubstituicao com a matriz L, começando de cima
void retrossubs_v_L (double** m_A, double* v_X, double* v_B, int n){
	double* v_res = (double *) calloc (n, sizeof(double));

	for (int i = 0 ; i < n ; i++){					//i = linha calculando o valor de v_X
		v_res[i] = v_B[i];

		for (int j = i-1 ; j >= 0 ; j--)			//j = coluna passando o valor de v_A[i][j]*v_X[i]
			v_res[i] -= m_A[i][j] * v_X[j];
		
		v_res[i] /= m_A[i][i];
	}

	for (int i = 0 ; i < n ; i++)
		v_X[i] = v_res[i];

	free (v_res);
	//agora, v_X possui o resultado do sistema linear
}

//soma os valores de delta em v_X de acordo com seus índices
void passa_delta_pra_X (t_i_double* v_delta, double* v_X, int n){
	for (int i = 0 ; i < n ; i++){
		int j = 0;
		while (v_delta[j].i != i){
			j++;
		}

		//aqui, encontramos o elemento de delta com o vetor correspondente

		v_X[i] += v_delta[j].n;
	}
	//agora, v_X possui x_V + v_delta;
}

//dá free na matriz LU
void free_v_m_LU (double* v_fun_it, double** m_fun_U, double** m_fun_L, int n){
	free (v_fun_it);
	for (int i = 0 ; i < n ; i++){
		free(m_fun_U[i]);
		free(m_fun_L[i]);
	}
	free(m_fun_L);
	free(m_fun_U);
}

//faz o newton padrão
double* newton_modificado (t_entrada* entrada, int* num_it, t_tempos* tempos){
	tempos->tmp_total = timestamp();

	void* 	funcao 	= entrada->funcao;
	char** 	v_vars;
	int 	n_vars;
	evaluator_get_variables (funcao, &v_vars, &n_vars);					//pega as informações da funcao atual
	
	void**  v_deriv;
	void*** m_deriv;
	alloca_v_m_derivs (&v_deriv, &m_deriv, n_vars);						//alloca o vetor e a matriz de derivadas
	tempos->tmp_deriv = timestamp();
	calcula_derivadas (funcao, v_deriv, m_deriv, v_vars, n_vars);		//calcula as derivadas no vetor e na matriz
	tempos->tmp_deriv = timestamp() - tempos->tmp_deriv;

	double*  v_fun_it;													//guarda os f'(xi)
	double** m_fun_U;													//guarda os f''(xi), no começo inteira e depois fica apenas a U
	double** m_fun_L;													//guarda a matriz L
	alloca_v_mLU (&v_fun_it, &m_fun_U, &m_fun_L, n_vars);				//alloca o vetor e a matriz que guardam o f'(xi) e f''(xi)

	t_i_double* v_delta;												//vetor delta que será calculado
	double* 	v_X;													//vetor que guarda os valores de x pra iteracao atual
	alloca_v_delta_X (&v_delta, &v_X, n_vars);							//aloca espaco pro vetor delta, x iteracao atual e x prox iteracao
	inicia_X_delta (v_X, v_delta, entrada);								//inicia os vetores

	double* v_res_it;													//guarda os resultados referentes às iteracoes
	alloca_v_double (&v_res_it, entrada->iteracoes);						//aloca o v_res_it
	double* v_Y;														//guarda o resultado intermediário da fatoração LU
	alloca_v_double (&v_Y, entrada->n_var);								//aloca o v_Y	

	int cont_it 	  = 0;
	v_res_it[cont_it] = evaluator_evaluate(funcao, n_vars, v_vars, v_X);
	int hess_steps 	  = entrada->n_var;

	double aux_tempos;
	tempos->tmp_SL = 0.0;

	while (cont_it < entrada->iteracoes){

		aux_tempos = timestamp();
		calcula_val_deriv_v (v_deriv, v_fun_it, n_vars, v_vars, v_X);

		cont_it++;
		if (norma(v_fun_it, n_vars) < entrada->epslon)
			break;

		if (((cont_it-1) % hess_steps) == 0){									//calcula a LU apenas a cada hess steps
			calcula_val_deriv_LU 	(m_deriv, m_fun_U, n_vars, v_vars, v_X);
			reordena_v_i_double 	(v_delta, n_vars);							//reordena o delta da ultima fatoração LU para ele guardar a da próxima
			LU_pivot				(m_fun_U, m_fun_L, v_delta, n_vars);		//eliminacao de gauss com pivoteamento na fat LU, com o delta guardando as trocas
		}

		troca_v_fun_it (v_fun_it, v_delta, n_vars);								//troca os elementos do v_fun de acordo com as torcas que aconteceram na LU

		retrossubs_v_L        (m_fun_L, v_Y, v_fun_it, n_vars);					//Ly = b
		retrossubs_v_i_double (m_fun_U, v_delta, v_Y, n_vars);					//Ux = y

		passa_delta_pra_X (v_delta, v_X, n_vars);

		aux_tempos = timestamp() - aux_tempos;
		tempos->tmp_SL += aux_tempos;

		v_res_it[cont_it] = evaluator_evaluate(funcao, n_vars, v_vars, v_X);

		if (isnan(v_res_it[cont_it]) || isinf(v_res_it[cont_it]))				//quando a funcao encontra um infinito ou um not a number, para de tentar
			break;

		if (norma_i(v_delta, n_vars) < entrada->epslon)
			break;
	}

	free_v_m_derivs 	(v_deriv, m_deriv, n_vars);
	free_v_m_LU		 	(v_fun_it, m_fun_U, m_fun_L, n_vars);
	free_v_delta_X 		(v_delta, v_X);
	free				(v_Y);

	*num_it = cont_it;
	tempos->tmp_total = timestamp() - tempos->tmp_total;
	return (v_res_it);
}