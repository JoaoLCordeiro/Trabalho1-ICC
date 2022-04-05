//LIBNI: Define as funcoes usadas pelo método de Newton inexato

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

//funcao que inicia o vetor x com seus valores iniciais
void inicia_x (double* v_X, int n, t_entrada* entrada){
	for (int i = 0 ; i < n ; i++){
		v_X[i] = entrada->valores_ini[i];								//carrega v_X com os valores iniciais
	}
}

//copia o x em delta
void copia_X_delta (double* v_X, double* v_delta, int n){
	for (int i = 0 ; i < n ; i++)
		v_delta[i] = v_X[i];
}

//torna o delta zero
void zera_delta (double* v_delta, int n){
	for (int i = 0 ; i < n ; i++)
		v_delta[i] = 0.0;
}

//soma delta no valor atual de X
void soma_delta_X (double* v_X, double* v_delta, int n){
	for (int i = 0 ; i < n ; i++)
		v_X[i] += v_delta[i];
}

//faz o newton inexato
double* newton_inexato (t_entrada* entrada, int* num_it, t_tempos* tempos){
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

	double*		v_X;
	double* 	v_delta;												//vetor que guarda os valores de x pra iteracao i e i+1
	alloca_v_double	(&v_X, n_vars);
	alloca_v_double	(&v_delta,n_vars);

	inicia_x		(v_X, n_vars, entrada);

	double* v_res_it;													//guarda os resultados referentes às iteracoes
	alloca_v_double (&v_res_it, entrada->iteracoes);					//aloca o v_res_it

	int 	cont_it	  = 0;												//contador de iteracoes
	v_res_it[cont_it] = evaluator_evaluate(funcao, n_vars, v_vars, v_X);

	double aux_tempos;

	while (cont_it < entrada->iteracoes){
		aux_tempos = timestamp();

		cont_it++;
		calcula_valores_deriv (v_deriv, m_deriv, v_f_iteracao, m_f_iteracao, n_vars, v_vars, v_X);
		
		double norma_gs = entrada->epslon + 1.0;
		int it = 0;
		while ((norma_gs > entrada->epslon)&&(it < 50)){
			norma_gs = 0.0;

			zera_delta   (v_delta, n_vars);
			// começa o Gauss-Seidal
			for (int i = 0 ; i < entrada->n_var ; i++){
				double sum = 0.0;
				// Primeiro lado da matriz, até a diagonal principal
				for (int j = 0; j < i; j++)
					sum += m_f_iteracao[i][j] * v_delta[j];

				// Segundo lado da matriz, depois da diagonal principal
				for (int j = i+1; j < entrada->n_var; j++)
					sum += m_f_iteracao[i][j] * v_delta[j];

				// Cálculo para a diagonal principal
				v_delta[i] = (v_f_iteracao[i] - sum) / m_f_iteracao[i][i];

				if (fabs(v_delta[i] - v_X[i]) > norma_gs)
					norma_gs = fabs(v_delta[i] - v_X[i]);
			}

			it++;
		}
		//termina o Gauss-Seidal

		aux_tempos = timestamp() - aux_tempos;
		tempos->tmp_SL += aux_tempos;

		soma_delta_X (v_X, v_delta, n_vars);							//soma o delta no valor atual de X, formando o X da prox iteracao

		v_res_it[cont_it] = evaluator_evaluate(funcao, n_vars, v_vars, v_X);

		if (isnan(v_res_it[cont_it]) || isinf(v_res_it[cont_it]))		//quando a funcao encontra um infinito ou um not a number, para de tentar
			break;

		if (norma (v_f_iteracao, n_vars) < entrada->epslon)
			break;
		
	}

	free_v_m_derivs 	(v_deriv, m_deriv, n_vars);
	free_v_m_funcao_it 	(v_f_iteracao, m_f_iteracao, n_vars);	
	free				(v_X);
	free 				(v_delta);

	*num_it = cont_it;
	tempos->tmp_total = timestamp() - tempos->tmp_total;
	return (v_res_it);
}