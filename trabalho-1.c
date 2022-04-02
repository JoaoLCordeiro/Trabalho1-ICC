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

int encontra_pivo (double** matriz, int coluna, int n){	//função que encontra o índice do maior número da
	int max 	= matriz[coluna][coluna];				//coluna para baixo da linha = coluna
	int max_i	= coluna;

	for (int i = coluna+1 ; i < n ; i++){
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
		m_A[i_1][j]	= m_A[i_2][j];		//ALTERAR PARA TROCAR APENAS O PONTEIRO
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

double norma (double*  v_valores, int n){
	/*double total = 0.0;
	for (int i = 0 ; i < n ; i++)
		total += v_valores[i]*v_valores[i];
	
	total = sqrt(total);
	printf("Norma:\n%le\n",total);
	return (total);*/
	double max, dif_atual;

	if (n == 1){
		printf("Norma:\n%le\n",fabs(v_valores[0]));
		return (fabs(v_valores[0]));
	}
	else{
		max	= v_valores[1] - v_valores[0];
		max = fabs (max);

		for (int i = 2 ; i < n ; i++){
			dif_atual	= v_valores[i] - v_valores[i-1];
			dif_atual	= fabs(dif_atual);

			if (dif_atual > max)
				max		= dif_atual;
		}
		printf("Norma:\n%le\n",max);
	}
	return (max);
}

void copia_vx2_vx1 (double* v_X_i2, double* v_X_i1, int n){
	for (int i = 0 ; i < n ; i++)
		v_X_i1[i] = v_X_i2[i];
}

void soma_x1_delta_pro_x2 (t_i_double* v_delta, double* v_X_i1, double* v_X_i2, int n){
	for (int i = 0 ; i < n ; i++){
		v_X_i2[i] = v_delta[i].n + v_X_i1[i];
	}
}

void imprime_vetor (double* vetor, int n){
	for (int i = 0 ; i < n ; i++)
		printf("%le	", vetor[i]);
	printf("\n");
}

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

double* newton_padrao (t_entrada* entrada){
	void* 	funcao 	= entrada->funcao;
	char** 	v_vars;
	int 	n_vars;
	evaluator_get_variables (funcao, &v_vars, &n_vars);					//pega as informações da funcao atual

	void**  v_deriv		= (void **)  malloc(n_vars*sizeof(void *));		//malloca um vetor de funcoes derivadas primeiras
	void*** m_deriv		= (void ***) malloc(n_vars*sizeof(void **));	//malloca uma matriz de funcoes derivadas segundas
	for (int i = 0 ; i < n_vars ; i++)
		m_deriv[i]		= (void **)  malloc(n_vars*sizeof(void *));

	calcula_derivadas (funcao, v_deriv, m_deriv, v_vars, n_vars);

	//DAR FREE NO D_FUNCAO E DD_FUNCAO

	double*  v_f_iteracao	= (double *)  malloc (n_vars*sizeof(double));	//guarda valores atuais do gradiente para os valroes de x atuais
	double** m_f_iteracao	= (double **) malloc (n_vars*sizeof(double*));	//guarda valores atuais da hessiana para os valroes de x atuais
	for (int i = 0 ; i < n_vars ; i++)
		m_f_iteracao[i]		= (double *)  malloc (n_vars*sizeof(double));

	t_i_double* v_delta	= (t_i_double *) malloc (n_vars*sizeof(t_i_double));//vetor delta que será calculado
	double* v_X_i1		= (double *) malloc (n_vars*sizeof(double));		//vetor que guarda os valores de x pra iteracao atual
	double* v_X_i2		= (double *) malloc (n_vars*sizeof(double));		//vetor que guarda os valores de x pra proxima iteracao

	for (int i = 0 ; i < n_vars ; i++){
		v_X_i1[i]		= entrada->valores_ini[i];							//carrega x_i1 com os valores iniciais
		v_delta[i].i	= i;												//coloca índices no vetor delta
	}

	int 	cont_it	= 0;		//contador de iteracoes

	while (cont_it < entrada->iteracoes){
	
		for (int i = 0 ; i < n_vars ; i++){												//gaurda os valores das funcoes derivadas
			v_f_iteracao[i]	= evaluator_evaluate(v_deriv[i], n_vars, v_vars, v_X_i1);	//com os X atuais
			v_f_iteracao[i] = -v_f_iteracao[i];											//deixa o valor com sinal trocado para operação
			for (int j = 0 ; j < n_vars ; j++){
				m_f_iteracao[i][j] = evaluator_evaluate(m_deriv[i][j], n_vars, v_vars, v_X_i1);
			}
		}

		printf("Iteracao numero %d:\n\n", cont_it);
		printf("vetor evaluate:\n");
		imprime_vetor (v_f_iteracao, n_vars);

		printf("matriz evaluate:\n");
		for (int i = 0 ; i < n_vars ; i++)
			imprime_vetor (m_f_iteracao[i], n_vars);

		if (norma(v_f_iteracao, n_vars) < entrada->epslon)
			break;

		resolve_sistema_linear 	(m_f_iteracao, v_delta, v_f_iteracao, n_vars);
		reordena_v_delta 		(v_delta, n_vars);
		soma_x1_delta_pro_x2	(v_delta, v_X_i1, v_X_i2, n_vars);

		copia_vx2_vx1(v_X_i2, v_X_i1, n_vars);				//copia i2 pra i1, logo, i1 agr possui o vetor da iteracao futura

		printf("Vetor x_i1:\n");
		imprime_vetor			(v_X_i1, n_vars);

		if (norma(v_X_i1, n_vars) < entrada->epslon)
			break;

		printf("\n\n");

		cont_it++;
	}

	return (v_X_i1);

	//imprime_v_m (v_deriv, m_deriv, n_vars); funcao de debug: imprime as funcoes derivadas encontradas
	//dar free no final nesses vetores/matrizes
	
}

void imprime_resultados (double* v_res, int n){
	printf("|Variável	|Resultado		|\n");
	for (int i = 0 ; i < n ; i++){
		printf("|%d		|%le		|\n", i, v_res[i]);
	}
	printf("\n");
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
	t_entrada* 	entrada_atual 	= (t_entrada *) malloc (sizeof(t_entrada));

	int cont = 0;

	while ((le_entrada(entrada_atual) == 0)&&(cont < 1)){
		//imprime_entrada (entrada_atual);
		double*		v_resultados	= (double *) malloc (entrada_atual->n_var*sizeof(double));

		v_resultados	= newton_padrao (entrada_atual);
		imprime_resultados (v_resultados, entrada_atual->n_var);

		cont++;
		free (v_resultados);
		//free_entrada(); dar free na entrada para a próxima
	}

	return 1;
}