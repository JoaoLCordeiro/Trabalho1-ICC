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

/*---------------------COMECA AS FUNCOES GERAIS--------------------*/

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

//retrossubstituicao em um vetor de double com índice
void retrossubs_v_i_double (double** m_A, t_i_double* v_X, double* v_B, int n){
	for (int i = n-1 ; i >=0 ; i--){					//i = linha calculando o valor de v_X
		v_X[i].n = v_B[i];

		for (int j = i+1 ; j < n ; j++)					//j = coluna passando o valor de v_A[i][j]*v_X[i]
			v_X[i].n -= m_A[i][j] * v_X[j].n;
		
		v_X[i].n /= m_A[i][i];
	}
	//agora, v_X possui o resultado do sistema linear
}

//o delta pode ser desordenado na hora do pivoteamento, essa funcao o ordena de volta
void reordena_v_i_double (t_i_double* v_i_double, int n){
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
double norma (double*  v_valores, int n){
	double total = 0.0;
	for (int i = 0 ; i < n ; i++)
		total += v_valores[i]*v_valores[i];
	
	total = sqrt(total);
	return (total);
}

//calcula a norma de um vetor de double com indice
double norma_i (t_i_double*  v_valores, int n){
	double total = 0.0;
	for (int i = 0 ; i < n ; i++)
		total += v_valores[i].n*v_valores[i].n;
	
	total = sqrt(total);
	return (total);
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

//aloca um vetor de valores de resultados relativos a cada iteração
void aloca_v_double (double** v_double, int n){
	*v_double	= (double *) calloc (n,sizeof(double));
}

//aloca os vetores: delta, xi1 e xi2
void alloca_v_delta_X (t_i_double** v_delta, double** v_X, int  n_vars){
	*v_delta	= (t_i_double *) calloc (n_vars,sizeof(t_i_double));
	aloca_v_double (v_X, n_vars);
}

//inicia o vetor xi1 com os valroes iniciais e coloca índice nos valores de delta
void inicia_X_delta (double* v_X_i1, t_i_double* v_delta, t_entrada* entrada){
	for (int i = 0 ; i < entrada->n_var ; i++){
		v_X_i1[i]		= entrada->valores_ini[i];							//carrega x_i1 com os valores iniciais
		v_delta[i].i	= i;												//coloca índices no vetor delta
	}
}

/*---------------------ACABA AS FUNCOES GERAIS---------------------*/

/*----------------COMECA AS FUNCOES DO NEWTON PADRAO---------------*/

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
	aloca_v_double (v_f_iteracao, n_vars);	
	*m_f_iteracao	= (double **) calloc (n_vars,sizeof(double*));	
	for (int i = 0 ; i < n_vars ; i++)
		aloca_v_double (&((*m_f_iteracao)[i]), n_vars);
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

void imprime_matriz (void*** m_deriv, int n, char** v_vars, double* v_X){
	printf ("MATRIZ NEWTON PADRÃO:\n");
	for (int i = 0 ; i < n ; i++){
		for (int j = 0 ; j < n ; j++)
			printf("%2.le	",evaluator_evaluate(m_deriv[i][j], n, v_vars, v_X));
		printf("\n");
	}
	printf("\n");
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
	double* 	v_X;													//vetor que guarda os valores de x pra iteracao
	alloca_v_delta_X (&v_delta, &v_X, n_vars);							//aloca espaco pro vetor delta e x da iteracao
	inicia_X_delta (v_X, v_delta, entrada);								//inicia os vetores

	double* v_res_it;													//guarda os resultados referentes às iteracoes
	aloca_v_double (&v_res_it, entrada->iteracoes);						//aloca o v_res_it

	int 	cont_it	  = 0;												//contador de iteracoes
	v_res_it[cont_it] = evaluator_evaluate(funcao, n_vars, v_vars, v_X);

	while (cont_it < entrada->iteracoes){
	
		cont_it++;
		calcula_valores_deriv (v_deriv, m_deriv, v_f_iteracao, m_f_iteracao, n_vars, v_vars, v_X);

		if (norma(v_f_iteracao, n_vars) < entrada->epslon)
			break;

		resolve_sistema_linear 	(m_f_iteracao, v_delta, v_f_iteracao, n_vars);
		reordena_v_i_double 	(v_delta, n_vars);
		soma_x_delta_pro_x  	(v_delta, v_X, n_vars);					//agora, v_X possui os X da futura iteracao (i+1)

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
	return (v_res_it);
}

/*----------------ACABA AS FUNCOES DO NEWTON PADRÃO----------------*/

/*--------------COMECA AS FUNCOES DO NEWTON MODIFICADO-------------*/

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

void free_v_m_LU (double* v_fun_it, double** m_fun_U, double** m_fun_L, int n){
	free (v_fun_it);
	for (int i = 0 ; i < n ; i++){
		free(m_fun_U[i]);
		free(m_fun_L[i]);
	}
	free(m_fun_L);
	free(m_fun_U);
}

void imprime_LU (void*** m_deriv, double** m_fun_L, double** m_fun_U, int n, char** v_vars, double* v_X){
	printf ("MATRIZ NEWTON MODIFICADO:\n");
	for (int i = 0 ; i < n ; i++){
		for (int j = 0 ; j < n ; j++)
			printf("%2.le	",evaluator_evaluate(m_deriv[i][j], n, v_vars, v_X));
		printf("\n");
	}
	printf("\n");

	/*printf ("L:\n");
	for (int i = 0 ; i < n ; i++){
		for (int j = 0 ; j < n ; j++)
			printf("%2.le	",m_fun_L[i][j]);
		printf("\n");
	}
	printf("\n");

	printf ("U:\n");
	for (int i = 0 ; i < n ; i++){
		for (int j = 0 ; j < n ; j++)
			printf("%2.le	",m_fun_U[i][j]);
		printf("\n");
	}
	printf("\n");*/
}

double* newton_modificado (t_entrada* entrada, int* num_it){
	void* 	funcao 	= entrada->funcao;
	char** 	v_vars;
	int 	n_vars;
	evaluator_get_variables (funcao, &v_vars, &n_vars);					//pega as informações da funcao atual
	
	void**  v_deriv;
	void*** m_deriv;
	alloca_v_m_derivs (&v_deriv, &m_deriv, n_vars);						//alloca o vetor e a matriz de derivadas
	calcula_derivadas (funcao, v_deriv, m_deriv, v_vars, n_vars);		//calcula as derivadas no vetor e na matriz

	double*  v_fun_it;													//guarda os f'(xi)
	double** m_fun_U;													//guarda os f''(xi), no começo inteira e depois fica apenas a U
	double** m_fun_L;													//guarda a matriz L
	alloca_v_mLU (&v_fun_it, &m_fun_U, &m_fun_L, n_vars);				//alloca o vetor e a matriz que guardam o f'(xi) e f''(xi)

	t_i_double* v_delta;												//vetor delta que será calculado
	double* 	v_X;													//vetor que guarda os valores de x pra iteracao atual
	alloca_v_delta_X (&v_delta, &v_X, n_vars);							//aloca espaco pro vetor delta, x iteracao atual e x prox iteracao
	inicia_X_delta (v_X, v_delta, entrada);								//inicia os vetores

	double* v_res_it;													//guarda os resultados referentes às iteracoes
	aloca_v_double (&v_res_it, entrada->iteracoes);						//aloca o v_res_it
	double* v_Y;														//guarda o resultado intermediário da fatoração LU
	aloca_v_double (&v_Y, entrada->n_var);								//aloca o v_Y	

	int cont_it 	  = 0;
	v_res_it[cont_it] = evaluator_evaluate(funcao, n_vars, v_vars, v_X);
	int hess_steps 	  = entrada->n_var;

	while (cont_it < entrada->iteracoes){

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
	return (v_res_it);
}

/*--------------ACABA AS FUNCOES DO NEWTON MODIFICADO--------------*/

//funcao que dá free na entrada
void free_entrada (t_entrada* entrada){
	evaluator_destroy(entrada->funcao);
	free (entrada->valores_ini);
}

int retorna_max (int a, int b){
	if (a > b)
		return a;
	else
		return b;
}

void imprime_resultados (double* v_res_np, double* v_res_nm, int num_it_np, int num_it_nm, int n_f){
	int max = retorna_max (num_it_np, num_it_nm);

	printf("Funcao número %d\n", n_f);
	printf("|Iteracao	||Newton Padrão	||Newton Modificado	|\n");
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

		printf("\n");
	}
	printf("\n");
}

int main(){
	t_entrada* 	entrada_atual 	= (t_entrada *) calloc (1,sizeof(t_entrada));

	int cont_f = 0;
	int num_it_np, num_it_nm;

	while (le_entrada(entrada_atual) == 0){
		double*		v_res_np;					//vetor de resultados do newton padrão
		double*		v_res_nm;					//vetor de resultados do newton modificado

		//printf("\nFUNCAO NUMERO %d:\n\n", cont_f);

		v_res_np	= newton_padrao 	(entrada_atual, &num_it_np);
		v_res_nm 	= newton_modificado (entrada_atual, &num_it_nm);
		
		imprime_resultados (v_res_np, v_res_nm, num_it_np, num_it_nm, cont_f);

		free (v_res_np);
		free (v_res_nm);
		free_entrada(entrada_atual); 			//dar free na entrada para a próxima

		cont_f++;
	}
	
	free (entrada_atual);
	return 1;
}