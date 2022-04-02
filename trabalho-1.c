#include <stdio.h>
#include <stdlib.h>
#include <matheval.h>
#include <assert.h>
#include <string.h>
#include <math.h>

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

//função que lê a entrada
int le_entrada(t_entrada* entrada){
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

	entrada->epslon;
	scanf("%le\n", &(entrada->epslon));												//lê o número para a tolerância epslon

	entrada->iteracoes;
	scanf("%d\n", &(entrada->iteracoes));											//lê o número máximo de iterações
	return (0);
}

/*----------------COMECA AS FUNCOES DO NEWTON PADRAO---------------*/

//função que encontra o índice do maior número da coluna para baixo da linha = coluna
int encontra_pivo (double** matriz, int coluna, int n){
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

//troca as linhas do sistema linear
void troca_linhas (double** m_A, double* v_B, t_i_double* v_X, int i_1, int i_2, int n){
	double aux;
	for (int j = 0 ; j < n ; j++){
		aux			= m_A[i_1][j];		//troca as linhas da matriz
		m_A[i_1][j]	= m_A[i_2][j];		//ALTERAR PARA TROCAR APENAS O PONTEIRO
		m_A[i_2][j] = aux;
	}

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

	for (int i = n-1 ; i >=0 ; i--){					//i = linha calculando o valor de v_X
		v_X[i].n = v_B[i];

		for (int j = i+1 ; j < n ; j++)					//j = coluna passando o valor de v_A[i][j]*v_X[i]
			v_X[i].n -= m_A[i][j] * v_X[j].n;
		
		v_X[i].n /= m_A[i][i];
	}
	//agora, o sistema está resolvido em v_X, com seus respectivos indices, mesmo após as trocas
}

//o delta pode ser desordenado na hora do pivoteamento, essa funcao o ordena de volta
void reordena_v_delta (t_i_double* v_delta, int n){
	//selection_sort onde os elementos são os índices dos elementos de v_delta (v_delta[j].i)
	int menor, menor_i, aux;

	for (int j = 0 ; j < n ; j++){					//ordena o v_delta em relação ao índice
		menor 	= v_delta[j].i;						//dos elementos
		menor_i	= j;

		for (int k = j+1 ; k < n ; k++){
			if (v_delta[k].i < menor){
				menor 	= v_delta[k].i;
				menor_i = k;
			}
		}

		if (j != menor_i){
			aux					= v_delta[j].n;			//troca os números
			v_delta[j].n		= v_delta[menor_i].n;
			v_delta[menor_i].n	= aux;

			aux					= v_delta[j].i;			//troca os índices
			v_delta[j].i		= v_delta[menor_i].i;
			v_delta[menor_i].i	= aux;
		}
	}
}

//calcula a norma de um vetor de valores
double norma (double*  v_valores, int n){
	double total = 0.0;
	for (int i = 0 ; i < n ; i++)
		total += v_valores[i]*v_valores[i];
	
	total = sqrt(total);
	return (total);
}

//copia o vetor vx2 no vx1
void copia_vx2_vx1 (double* v_X_i2, double* v_X_i1, int n){
	for (int i = 0 ; i < n ; i++)
		v_X_i1[i] = v_X_i2[i];
}

//soma o vetor x1 e o delta no vetor x2
void soma_x1_delta_pro_x2 (t_i_double* v_delta, double* v_X_i1, double* v_X_i2, int n){
	for (int i = 0 ; i < n ; i++){
		v_X_i2[i] = v_delta[i].n + v_X_i1[i];
	}
}

//imprime um vetor, não usada no final APAGAR ESSA FUNCAO AQUI
void imprime_vetor (double* vetor, int n){
	for (int i = 0 ; i < n ; i++)
		printf("%le	", vetor[i]);
	printf("\n");
}

//calcula as funcoes derivadas da funcao "funcao"
void calcula_derivadas (void* funcao, void** v_deriv, void*** m_deriv, char** v_vars, int n){
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
void alloca_v_m_derivs (void*** v_deriv, void**** m_deriv, int n){
	*v_deriv	= (void **)  calloc(n,sizeof(void *));		//malloca um vetor de funcoes derivadas primeiras
	*m_deriv	= (void ***) calloc(n,sizeof(void **));		//malloca uma matriz de funcoes derivadas segundas
	for (int i = 0 ; i < n ; i++)
		(*m_deriv)[i]	= (void **)  calloc(n,sizeof(void *));
}

//aloca um vetor e uma matriz de valores de funcoes de tamanho n
void alloca_v_m_funcao_it (double** v_f_iteracao, double*** m_f_iteracao, int n_vars){
	*v_f_iteracao	= (double *)  calloc (n_vars,sizeof(double));	
	*m_f_iteracao	= (double **) calloc (n_vars,sizeof(double*));	
	for (int i = 0 ; i < n_vars ; i++)
		(*m_f_iteracao)[i]		= (double *)  calloc (n_vars,sizeof(double));
}

//aloca os vetores: delta, xi1 e xi2
void alloca_v_delta_X (t_i_double** v_delta, double** v_X_i1, double** v_X_i2,int  n_vars){
	*v_delta	= (t_i_double *) calloc (n_vars,sizeof(t_i_double));
	*v_X_i1		= (double *) 	 calloc (n_vars,sizeof(double));
	*v_X_i2		= (double *) 	 calloc (n_vars,sizeof(double));
}

//inicia o vetor xi1 com os valroes iniciais e coloca índice nos valores de delta
void inicia_Xi1_delta (double* v_X_i1, t_i_double* v_delta, t_entrada* entrada){
	for (int i = 0 ; i < entrada->n_var ; i++){
		v_X_i1[i]		= entrada->valores_ini[i];							//carrega x_i1 com os valores iniciais
		v_delta[i].i	= i;												//coloca índices no vetor delta
	}
}

//aloca um vetor de valores de resultados relativos a cada iteração
void aloca_v_res_it (double** v_res_it, int it){
	*v_res_it	= (double *) calloc (it,sizeof(double));
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

//faz o newton padrao
double* newton_padrao (t_entrada* entrada, int* num_it){
	void* 	funcao 	= entrada->funcao;
	char** 	v_vars;
	int 	n_vars;
	evaluator_get_variables (funcao, &v_vars, &n_vars);					//pega as informações da funcao atual

	void**  v_deriv;
	void*** m_deriv;
	alloca_v_m_derivs (&v_deriv, &m_deriv, n_vars);						//alloca o vetor e a matriz de derivadas
	calcula_derivadas (funcao, v_deriv, m_deriv, v_vars, n_vars);		//calcula as derivadas no vetor e na matriz

	double*  v_f_iteracao;												//guarda os f'(xi)
	double** m_f_iteracao;												//guarda os f''(xi)
	alloca_v_m_funcao_it (&v_f_iteracao, &m_f_iteracao, n_vars);		//alloca o vetor e a matriz que guardam o f'(xi) e f''(xi)

	t_i_double* v_delta;												//vetor delta que será calculado
	double* 	v_X_i1;													//vetor que guarda os valores de x pra iteracao atual
	double* 	v_X_i2;													//vetor que guarda os valores de x pra proxima iteracao
	alloca_v_delta_X (&v_delta, &v_X_i1, &v_X_i2, n_vars);				//aloca espaco pro vetor delta, x iteracao atual e x prox iteracao
	inicia_Xi1_delta (v_X_i1, v_delta, entrada);						//inicia os vetores

	double* v_res_it;													//guarda os resultados referentes às iteracoes
	aloca_v_res_it (&v_res_it, entrada->iteracoes);						//aloca o v_res_it

	int 	cont_it	= 0;												//contador de iteracoes

	while (cont_it < entrada->iteracoes){
	
		cont_it++;
		calcula_valores_deriv (v_deriv, m_deriv, v_f_iteracao, m_f_iteracao, n_vars, v_vars, v_X_i1);

		if (norma(v_f_iteracao, n_vars) < entrada->epslon)
			break;

		resolve_sistema_linear 	(m_f_iteracao, v_delta, v_f_iteracao, n_vars);
		reordena_v_delta 		(v_delta, n_vars);
		soma_x1_delta_pro_x2	(v_delta, v_X_i1, v_X_i2, n_vars);

		copia_vx2_vx1(v_X_i2, v_X_i1, n_vars);				//copia i2 pra i1, logo, i1 agr possui o vetor da iteracao futura

		v_res_it[cont_it] = evaluator_evaluate(funcao, n_vars, v_vars, v_X_i1);
		if (norma(v_X_i1, n_vars) < entrada->epslon)
			break;

	}

	*num_it = cont_it;
	return (v_res_it);

	//imprime_v_m (v_deriv, m_deriv, n_vars); funcao de debug: imprime as funcoes derivadas encontradas
	//dar free no final nesses vetores/matrizes
	
}

/*----------------ACABA AS FUNCOES DO NEWTON PADRÃO----------------*/

//funcao que dá free na entrada
void free_entrada (t_entrada* entrada){
	free (entrada->valores_ini);
}

void imprime_resultados (double* v_res, int n, int n_f){
	printf("Funcao número %d\n", n_f);
	printf("|Iteracao	|Resultado		|\n");
	for (int i = 0 ; i < n ; i++){
		printf("|%d		|%le		|\n", i, v_res[i]);
	}
	printf("\n");
}

/*void imprime_entrada (t_entrada *entrada){				//funcao só pra debuggar :)
	printf ("%d\n", entrada->n_var);
	printf ("%s", entrada->funcao);
	for (int i = 0; i < entrada->n_var ; i++)
		printf ("%lf ", entrada->valores_ini[i]);
	printf ("\n");
	printf ("%le\n", entrada->epslon);
	printf ("%d\n", entrada->iteracoes);
}*/

int main(){
	t_entrada* 	entrada_atual 	= (t_entrada *) calloc (1,sizeof(t_entrada));

	int cont_f = 0;
	int num_it;

	while (le_entrada(entrada_atual) == 0){
		double*		v_resultados_np	= (double *) calloc (entrada_atual->n_var,sizeof(double));	//vetor de resultados do newton padrão

		v_resultados_np	= newton_padrao (entrada_atual, &num_it);
		printf("iteracoes: %d\n", num_it);
		imprime_resultados (v_resultados_np, num_it, cont_f);

		free (v_resultados_np);
		free_entrada(entrada_atual); 	//dar free na entrada para a próxima

		cont_f++;
	}

	return 1;
}