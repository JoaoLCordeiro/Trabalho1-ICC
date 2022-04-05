//LIBNP: Define as funcoes usadas pelo método de Newton padrão

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

//troca as linhas do sistema linear
void troca_linhas (double** m_A, double* v_B, t_i_double* v_X, int i_1, int i_2, int n){
	double 	aux;
	double* auxp;

	auxp		= m_A[i_1];				//troca as linhas da matriz
	m_A[i_1]	= m_A[i_2];
	m_A[i_2]	= auxp;

	aux 		= v_B[i_1];				//troca os valores do vetor B
	v_B[i_1] 	= v_B[i_2];
	v_B[i_2]	= aux;

	aux 		= v_X[i_1].n;			//troca os valores do vetor X
	v_X[i_1].n 	= v_X[i_2].n;
	v_X[i_2].n	= aux;

	aux 		= v_X[i_1].i;			//troca os indices do vetor X
	v_X[i_1].i 	= v_X[i_2].i;
	v_X[i_2].i	= aux;
}

//funcao que resolve o sistema linear dado
void resolve_sistema_linear (double** m_A, t_i_double* v_X, double* v_B, int n){
	//faz eliminação de gauss com pivoteamento parcial

	for (int i = 0 ; i < n ; i++){						//i = linha encontrando o pivo
		int iPivo	= encontra_pivo(m_A, i, n);

		if (i != iPivo)
			troca_linhas (m_A, v_B, v_X, i, iPivo, n);

		for (int j = i+1 ; j < n ; j++){				//j = linha subtraindo a linha do pivo
			double m 	= m_A[j][i]/m_A[i][i];
			m_A[j][i]	= 0.0;

			for (int k = i+1 ; k < n ; k++)				//k = coluna fazendo a subtracao por elemento
				m_A[j][k]	-= m_A[i][k]*m;
			
			v_B[j] -= v_B[i]*m;
		}
	}

	//faz retrosubstituicao
	retrossubs_v_i_double (m_A, v_X, v_B, n);

	//agora, o sistema está resolvido em v_X, com seus respectivos indices, mesmo após as trocas
}

//soma o vetor x1 e o delta no vetor x2
void soma_x_delta_pro_x (t_i_double* v_delta, double* v_X, int n){
	for (int i = 0 ; i < n ; i++){
		v_X[i] += v_delta[i].n;
	}
}

//aloca um vetor e uma matriz de valores de funcoes de tamanho n
void alloca_v_m_funcao_it (double** v_f_iteracao, double*** m_f_iteracao, int n_vars){
	alloca_v_double (v_f_iteracao, n_vars);	
	*m_f_iteracao	= (double **) calloc (n_vars,sizeof(double*));	
	for (int i = 0 ; i < n_vars ; i++)
		alloca_v_double (&((*m_f_iteracao)[i]), n_vars);
}

//calcula os valores das funcoes derivadas
void calcula_valores_deriv (void** v_deriv, void*** m_deriv, double* v_f_iteracao, double** m_f_iteracao, int n_vars, char** v_vars, double* v_X){
	for (int i = 0 ; i < n_vars ; i++){												//guarda os valores das funcoes derivadas
		v_f_iteracao[i]	= -evaluator_evaluate(v_deriv[i], n_vars, v_vars, v_X);		//com os X atuais e deixa o valor com sinal trocado para operação
		for (int j = 0 ; j < n_vars ; j++){
			m_f_iteracao[i][j] = evaluator_evaluate(m_deriv[i][j], n_vars, v_vars, v_X);
		}
	}
}

//dá free no vetor e matriz de funcoes derivadas
void free_v_m_derivs (void**  v_deriv, void*** m_deriv, int n){
	for (int i = 0 ; i < n ; i++){
		evaluator_destroy(v_deriv[i]);
		for (int j = 0 ; j < n ; j++)
			evaluator_destroy(m_deriv[i][j]);
		free (m_deriv[i]);
	}

	free (v_deriv);
	free (m_deriv);
}

//dá free no vetor e matriz de valores de funcoes derivadas
void free_v_m_funcao_it (double*  v_f_iteracao, double** m_f_iteracao, int n){
	free (v_f_iteracao);
	for (int i = 0 ; i < n ; i++)
		free (m_f_iteracao[i]);
	free (m_f_iteracao);
}

//dá free nos vetores X e delta
void free_v_delta_X (t_i_double* v_delta, double* v_X_i1){
	free (v_delta);
	free (v_X_i1);
}

//faz o newton padrao
double* newton_padrao (t_entrada* entrada, int* num_it, t_tempos* tempos){
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

	double*  v_f_iteracao;												//guarda os f'(xi)
	double** m_f_iteracao;												//guarda os f''(xi)
	alloca_v_m_funcao_it (&v_f_iteracao, &m_f_iteracao, n_vars);		//alloca o vetor e a matriz que guardam o f'(xi) e f''(xi)

	t_i_double* v_delta;												//vetor delta que será calculado
	double* 	v_X;													//vetor que guarda os valores de x pra iteracao
	alloca_v_delta_X (&v_delta, &v_X, n_vars);							//aloca espaco pro vetor delta e x da iteracao
	inicia_X_delta (v_X, v_delta, entrada);								//inicia os vetores

	double* v_res_it;													//guarda os resultados referentes às iteracoes
	alloca_v_double (&v_res_it, entrada->iteracoes);					//aloca o v_res_it

	int 	cont_it	  = 0;												//contador de iteracoes
	v_res_it[cont_it] = evaluator_evaluate(funcao, n_vars, v_vars, v_X);

	double aux_tempos;
	tempos->tmp_SL = 0.0;

	while (cont_it < entrada->iteracoes){
	
		cont_it++;

		aux_tempos = timestamp();
		calcula_valores_deriv (v_deriv, m_deriv, v_f_iteracao, m_f_iteracao, n_vars, v_vars, v_X);

		if (norma(v_f_iteracao, n_vars) < entrada->epslon)
			break;

		resolve_sistema_linear 	(m_f_iteracao, v_delta, v_f_iteracao, n_vars);
		reordena_v_i_double 	(v_delta, n_vars);
		soma_x_delta_pro_x  	(v_delta, v_X, n_vars);					//agora, v_X possui os X da futura iteracao (i+1)
		aux_tempos = timestamp() - aux_tempos;
		tempos->tmp_SL += aux_tempos;

		v_res_it[cont_it] = evaluator_evaluate(funcao, n_vars, v_vars, v_X);
	
		if (isnan(v_res_it[cont_it]) || isinf(v_res_it[cont_it]))		//quando a funcao encontra um infinito ou um not a number, para de tentar
			break;

		if (norma_i(v_delta, n_vars) < entrada->epslon)
			break;

	}

	free_v_m_derivs 	(v_deriv, m_deriv, n_vars);
	free_v_m_funcao_it 	(v_f_iteracao, m_f_iteracao, n_vars);
	free_v_delta_X 		(v_delta, v_X);

	*num_it = cont_it;
	tempos->tmp_total = timestamp() - tempos->tmp_total;
	return (v_res_it);
}