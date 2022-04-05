//LIBGERAL: Define as funcoes usadas por mais de um método e as structs usadas

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

#define MAX 1000

//função que lê a entrada
int 	le_entrada (t_entrada* entrada){
	if ((scanf("%d\n", &(entrada->n_var))) == EOF)									//lê o número de incógnitas
		return (-1);		//o arquivo acabou, pode parar :)

	char* funcao_string = (char *) malloc ((MAX+1)*sizeof(char));					//malloca e lê a string da funcao
	fgets(funcao_string, MAX, stdin);
	strtok(funcao_string, "\n");													//tira o \n do string para acertar a sintaxe para o evaluator
	entrada->funcao	= evaluator_create(funcao_string);								//cria a funcao e poe em funcao
	free (funcao_string);

	if (entrada->funcao == NULL){
		fprintf(stderr, "ERRO: sintaxe de funcao errada!\n");
		return (-2);		//erro na evaluacao: erro de sintaxe na funcao
	}

	printf ("%s", funcao_string);

	entrada->valores_ini = (double *) calloc (entrada->n_var,sizeof(double));		//alloca e lê um vetor de valores iniciais para as incógnitas

	if (entrada->valores_ini == NULL){
		fprintf(stderr, "ERRO: não foi possível alocar vetor de valores de incógnitas iniciais\n");
		return (-3);
	}

	for (int i = 0; i < entrada->n_var ; i++){
		scanf("%lf ",&(entrada->valores_ini[i]));
	}

	scanf("%le\n", &(entrada->epslon));												//lê o número para a tolerância epslon

	scanf("%d\n", &(entrada->iteracoes));											//lê o número máximo de iterações
	return (0);
}

//função que encontra o índice do maior número da coluna para baixo da linha = coluna
int 	encontra_pivo (double** matriz, int coluna, int n){
	int max 	= fabs (matriz[coluna][coluna]);
	int max_i	= coluna;

	for (int i = coluna+1 ; i < n ; i++){
		if (fabs (matriz[i][coluna]) > max){				//se o número analisado for maior que o máximo atual
			max 	= fabs (matriz[i][coluna]);				//se torna o max atual, guarda o índice dele também
			max_i	= i;
		}
	}

	return (max_i);
}

//retrossubstituicao em um vetor de double com índice
void 	retrossubs_v_i_double (double** m_A, t_i_double* v_X, double* v_B, int n){
	for (int i = n-1 ; i >=0 ; i--){					//i = linha calculando o valor de v_X
		v_X[i].n = v_B[i];

		for (int j = i+1 ; j < n ; j++)					//j = coluna passando o valor de v_A[i][j]*v_X[i]
			v_X[i].n -= m_A[i][j] * v_X[j].n;
		
		v_X[i].n /= m_A[i][i];
	}
	//agora, v_X possui o resultado do sistema linear
}

//o delta pode ser desordenado na hora do pivoteamento, essa funcao o ordena de volta
void 	reordena_v_i_double (t_i_double* v_i_double, int n){
	//selection_sort onde os elementos são os índices dos elementos de v_delta (v_delta[j].i)
	int menor, menor_i, aux;

	for (int j = 0 ; j < n ; j++){					//ordena o v_delta em relação ao índice
		menor 	= v_i_double[j].i;						//dos elementos
		menor_i	= j;

		for (int k = j+1 ; k < n ; k++){
			if (v_i_double[k].i < menor){
				menor 	= v_i_double[k].i;
				menor_i = k;
			}
		}

		if (j != menor_i){
			aux						= v_i_double[j].n;			//troca os números
			v_i_double[j].n			= v_i_double[menor_i].n;
			v_i_double[menor_i].n	= aux;

			aux						= v_i_double[j].i;			//troca os índices
			v_i_double[j].i			= v_i_double[menor_i].i;
			v_i_double[menor_i].i	= aux;
		}
	}
}

//calcula a norma de um vetor de double
double 	norma (double*  v_valores, int n){
	double total = 0.0;
	for (int i = 0 ; i < n ; i++)
		total += v_valores[i]*v_valores[i];
	
	total = sqrt(total);
	return (total);
}

//calcula a norma de um vetor de double com indice
double 	norma_i (t_i_double*  v_valores, int n){
	double total = 0.0;
	for (int i = 0 ; i < n ; i++)
		total += v_valores[i].n*v_valores[i].n;
	
	total = sqrt(total);
	return (total);
}

//calcula as funcoes derivadas da funcao "funcao"
void 	calcula_derivadas (void* funcao, void** v_deriv, void*** m_deriv, char** v_vars, int n){
	void* d_funcao;
	void* dd_funcao;

	for (int i = 0 ; i < n ; i++){										//enche o vetor v_deriv de derivadas
		d_funcao 	= evaluator_derivative(funcao, v_vars[i]);			//parciais primeiras
		v_deriv[i] 	= d_funcao;											//gradiente
	}

	for (int i = 0 ; i < n ; i++){										//enche a matriz m_deriv de derivadas
		d_funcao = v_deriv[i];											//parciais segundas
		for (int j = 0 ; j < n ; j++){
			dd_funcao 		= evaluator_derivative(d_funcao, v_vars[j]);
			m_deriv[i][j] 	= dd_funcao;								//hessiana
		}
	}
}

//aloca um vetor e uma matriz de funcoes derivadas de tamanho n
void 	alloca_v_m_derivs (void*** v_deriv, void**** m_deriv, int n){
	*v_deriv	= (void **)  calloc(n,sizeof(void *));		//malloca um vetor de funcoes derivadas primeiras
	*m_deriv	= (void ***) calloc(n,sizeof(void **));		//malloca uma matriz de funcoes derivadas segundas
	for (int i = 0 ; i < n ; i++)
		(*m_deriv)[i]	= (void **)  calloc(n,sizeof(void *));
}

//aloca um vetor de valores de resultados relativos a cada iteração
void 	alloca_v_double (double** v_double, int n){
	*v_double	= (double *) calloc (n,sizeof(double));
}

//aloca os vetores: delta e v_X
void 	alloca_v_delta_X (t_i_double** v_delta, double** v_X, int  n_vars){
	*v_delta	= (t_i_double *) calloc (n_vars,sizeof(t_i_double));
	alloca_v_double (v_X, n_vars);
}

//inicia o vetor xi1 com os valroes iniciais e coloca índice nos valores de delta
void 	inicia_X_delta (double* v_X_i1, t_i_double* v_delta, t_entrada* entrada){
	for (int i = 0 ; i < entrada->n_var ; i++){
		v_X_i1[i]		= entrada->valores_ini[i];							//carrega x_i1 com os valores iniciais
		v_delta[i].i	= i;												//coloca índices no vetor delta
	}
}

//funcao que dá free na entrada
void 	free_entrada (t_entrada* entrada){
	evaluator_destroy(entrada->funcao);
	free (entrada->valores_ini);
}

//retorna o maior dos valores dos parâmetros
int 	retorna_max (int a, int b, int c){
	if (a > b){
		if (a > c)
			return a;
	}
	else{
		if (b > c)
			return b;
	}
	return c;
}

//imprime os resultados do programa
void 	imprime_resultados (t_entrada* entrada, double* v_res_np, double* v_res_nm, double* v_res_ni, int num_it_np, int num_it_nm, int num_it_ni, t_tempos tempos_np, t_tempos tempos_nm, t_tempos tempos_ni){
	int max = retorna_max (num_it_np, num_it_nm, num_it_ni);

	printf("%d\n", entrada->n_var);
	printf("%s\n", evaluator_get_string(entrada->funcao));
	printf("#\n");
	printf("|Iteracao	||Newton Padrão	||Newton Modificado	||Newton Inexato	|\n");
	for (int i = 0 ; i < max ; i++){
		printf("|%d		|",i+1);
		if (i < num_it_np)
			printf("|%le	|",v_res_np[i]);
		else
			printf("|		|");

		if (i < num_it_nm)
			printf("|%le		|",v_res_nm[i]);
		else
			printf("|			|");

		if (i < num_it_ni)
			printf("|%le		|",v_res_ni[i]);
		else
			printf("|			|");

		printf("\n");
	}

	printf("|Tempo Total	||%le	||%le		||%le		|\n", tempos_np.tmp_total, tempos_nm.tmp_total, tempos_ni.tmp_total);
	printf("|Tempo Derivadas||%le	||%le		||%le		|\n", tempos_np.tmp_deriv, tempos_nm.tmp_deriv, tempos_ni.tmp_deriv);
	printf("|Tempo SL	||%le	||%le		||%le		|\n", 	  tempos_np.tmp_SL   , tempos_nm.tmp_SL   , tempos_ni.tmp_SL);
	printf("#\n");
	printf("\n");
}