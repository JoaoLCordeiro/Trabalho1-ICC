#include <stdio.h>
#include <stdlib.h>
#include <matheval.h>
#include <assert.h>
#include <string.h>

#define MAX 1000

struct t_entrada{			//tipo entrada: armazena os dados da entrada atual
	int		n_var;			//número de variáveis
	void*	funcao;			//string da funcao (TALVEZ ALTERAR AQUI)
	double*	valores_ini;	//vetor de valores iniciais para as incógnitas
	double	epslon;			//tolerância epslon
	int		iteracoes;		//número máximo de iterações
} typedef t_entrada;

struct t_i_double{			//um tipo double com índice, para o vetor delta
	int 	i;
	double 	n;
} typedef t_i_double;

int le_entrada(t_entrada* entrada){					//função que lê a entrada
	if ((scanf("%d\n", &(entrada->n_var))) == EOF)									//lê o número de incógnitas
		return (-1);		//erro na leitura: o arquivo acabou

	char* funcao_string = (char *) malloc ((MAX+1)*sizeof(char));					//malloca e lê a string da funcao
	fgets(funcao_string, MAX, stdin);
	strtok(funcao_string, "\n");													//tira o \n do string para acertar a sintaxe para o evaluator
	entrada->funcao	= evaluator_create(funcao_string);								//cria a funcao e poe em funcao
	free (funcao_string);

	if (entrada->funcao == NULL)
		return (-2);		//erro na evaluacao: erro de sintaxe na funcao

	printf ("%s", funcao_string);

	entrada->valores_ini = (double *) malloc (entrada->n_var * sizeof(double));		//malloca e lê um vetor de valores iniciais para as incógnitas
	for (int i = 0; i < entrada->n_var ; i++){
		scanf("%lf ",&(entrada->valores_ini[i]));
	}

	entrada->epslon;
	scanf("%le\n", &(entrada->epslon));												//lê o número para a tolerância epslon

	entrada->iteracoes;
	scanf("%d\n", &(entrada->iteracoes));											//lê o número máximo de iterações
	return (0);
}

void imprime_v_m (void** v_deriv, void*** m_deriv, int n_vars){
	for (int i = 0 ; i < n_vars ; i++){
		printf ("  f'[%d](x^(i))) = %s\n", i, evaluator_get_string (v_deriv[i]));
	}

	for (int i = 0 ; i < n_vars ; i++){
		for (int j = 0 ; j < n_vars ; j++){
			printf ("  f''[%d][%d](x^(i))) = %s\t", i, j, evaluator_get_string (m_deriv[i][j]));
		}
		printf("\n");
	}
}

int encontra_pivo (double** matriz, int coluna){		//função que encontra o índice do maior número da
	int max 	= matriz[coluna][coluna];				//coluna para baixo da linha = coluna
	int max_i	= coluna;

	for (int i = coluna+1 ; i < coluna ; i++){
		if (matriz[i][coluna] > max){					//se o número analisado for maior que o máximo atual
			max 	= matriz[i][coluna];				//se torna o max atual, guarda o índice dele também
			max_i	= i;
		}
	}

	return (max_i);
}

void troca_linhas (double** m_A, double* v_B, t_i_double* v_X, int i_1, int i_2, int n){
	double aux;
	for (int j = 0 ; j < n ; j++){
		aux			= m_A[i_1][j];		//troca as linhas da matriz
		m_A[i_1][j]	= m_A[i_2][j];
		m_A[i_2][j] = aux;
	}

	aux 		= v_B[i_1];				//troca os valores do vetor B
	v_B[i_1] 	= v_B[i_2];
	v_B[i_2]	= aux;

	aux 		= v_X[i_1].n;			//troca os valores do vetor_x
	v_X[i_1].n 	= v_X[i_2].n;
	v_X[i_2].n	= aux;

	aux 		= v_X[i_1].i;			//troca os indices do vetor_x
	v_X[i_1].i 	= v_X[i_2].i;
	v_X[i_2].i	= aux;
}

void resolve_sistema_linear (double** m_A, t_i_double* v_X, double* v_B, int n){
	//faz eliminação de gauss com pivoteamento parcial

	for (int i = 0 ; i < n ; i++){						//i = linha encontrando o pivo
		int iPivo	= encontra_pivo(m_A, i);

		if (i == iPivo)
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

	for (int i = n-1 ; i >=0 ; i--){					//i = linha calculando o valor de v_X
		v_X[i].n = v_B[i];

		for (int j = j+1 ; j < n ; j++)					//j = coluna passando o valor de v_A[i][j]*v_X[i]
			v_X[i].n -= m_A[i][j] * v_X[j].n;
		
		v_X[i].n /= m_A[i][i];
	}

	//agora, o sistema está resolvido em v_X, com seus respectivos indices, mesmo após as trocas
}

void reordena_v_delta (t_i_double* v_delta, int n){
	//selection_sort onde os elementos são os índices dos elementos de v_delta (v_delta[j].i)
	int menor, menor_i, aux;

	for (int j = 0 ; j < n ; j++){					//ordena o v_delta em relação ao índice
		menor 	= v_delta[j].i;						//dos elementos
		menor_i	= j;

		for (int k = j ; k < n ; k++){
			if (v_delta[j].i < menor){
				menor 	= v_delta[j].i;
				menor_i = k;
			}
		}

		aux					= v_delta[j].n;			//troca os números
		v_delta[j].n		= v_delta[menor_i].n;
		v_delta[menor_i].n	= aux;

		aux					= v_delta[j].i;			//troca os índices
		v_delta[j].i		= v_delta[menor_i].i;
		v_delta[menor_i].i	= aux;
	}
}

void newton_padrao (t_entrada* entrada){
	void* 	funcao 	= entrada->funcao;
	char** 	v_vars;
	int 	n_vars;
	evaluator_get_variables (funcao, &v_vars, &n_vars);					//pega as informações da funcao atual

	void**  v_deriv		= (void **)  malloc(n_vars*sizeof(void *));		//malloca um vetor de funcoes derivadas primeiras
	void*** m_deriv		= (void ***) malloc(n_vars*sizeof(void **));	//malloca uma matriz de funcoes derivadas segundas
	for (int i = 0 ; i < n_vars ; i++)
		m_deriv[i]		= (void **)  malloc(n_vars*sizeof(void *));

	void* d_funcao;
	void* dd_funcao;

	for (int i = 0 ; i < n_vars ; i++){									//enche o vetor v_deriv de derivadas
		d_funcao 	= evaluator_derivative(funcao, v_vars[i]);			//parciais primeiras
		v_deriv[i] 	= d_funcao;											//gradiente
	}

	for (int i = 0 ; i < n_vars ; i++){									//enche a matriz m_deriv de derivadas
		d_funcao = v_deriv[i];											//parciais segundas
		for (int j = 0 ; j < n_vars ; j++){
			dd_funcao 		= evaluator_derivative(d_funcao, v_vars[j]);
			m_deriv[i][j] 	= dd_funcao;								//hessiana
		}
	}

	double*  v_f_iteracao	= (double *)  malloc (n_vars*sizeof(double));	//guarda valores atuais do gradiente para os valroes de x atuais
	double** m_f_iteracao	= (double **) malloc (n_vars*sizeof(double*));	//guarda valores atuais da hessiana para os valroes de x atuais
	for (int i = 0 ; i < n_vars ; i++)
		m_f_iteracao[i]		= (double *)  malloc (n_vars*sizeof(double));

	t_i_double* v_delta	= (t_i_double *) malloc (n_vars*sizeof(t_i_double));//vetor delta que será calculado
	double* v_X_i1		= (double *) malloc (n_vars*sizeof(double));		//vetor que guarda os valores de x pra iteracao atual
	for (int i = 0 ; i < n_vars ; i++){
		v_X_i1[i]	= entrada->valores_ini[i];
		v_delta[i].i	= i;
	}
	double* v_X_i2		= (double *) malloc (n_vars*sizeof(double));		//vetor que guarda os valores de x pra proxima iteracao

	while (/*CONDIÇÃO DE PARADA*/ 1){
		for (int i = 0 ; i < n_vars ; i++){												//gaurda os valores das funcoes derivadas
			v_f_iteracao[i]	= evaluator_evaluate(v_deriv[i], n_vars, v_vars, v_X_i1);//com os X atuais
			for (int j = 0 ; j < n_vars ; j++){
				m_f_iteracao[i][j] = evaluator_evaluate(m_deriv[i][j], n_vars, v_vars, v_X_i1);
			}
		}

		resolve_sistema_linear (m_f_iteracao, v_delta, v_f_iteracao, n_vars);
		reordena_v_delta (v_delta, n_vars);
	}

	//imprime_v_m (v_deriv, m_deriv, n_vars); funcao de debug: imprime as funcoes derivadas encontradas
	//dar free no final nesses vetores/matrizes
}

void imprime_entrada (t_entrada *entrada){				//funcao só pra debuggar :)
	printf ("%d\n", entrada->n_var);
	printf ("%s", entrada->funcao);
	for (int i = 0; i < entrada->n_var ; i++)
		printf ("%lf ", entrada->valores_ini[i]);
	printf ("\n");
	printf ("%le\n", entrada->epslon);
	printf ("%d\n", entrada->iteracoes);
}

int main(){
	t_entrada* entrada_atual 	= (t_entrada *) malloc (sizeof(t_entrada));

	int cont = 0;

	while ((le_entrada(entrada_atual) == 0)&&(cont < 10)){
		
		//imprime_entrada (entrada_atual);

		newton_padrao (entrada_atual);
		cont++;
		//free_entrada(); dar free na entrada para a próxima
	}

	return 1;
}