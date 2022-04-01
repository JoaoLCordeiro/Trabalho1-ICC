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

void newton_padrao (t_entrada* entrada){
	void* funcao 	= entrada->funcao;
	char** v_vars;
	int n_vars;
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
		v_f_iteracao[i]		= (double *)  malloc (n_vars*sizeof(double));

	double* v_delta		= (double *) malloc (n_vars*sizeof(double));		//vetor delta que será calculado
	double* v_X_atual	= (double *) malloc (n_vars*sizeof(double));		//vetor que guarda os valores de x atuais
	for (int i = 0 ; i < n_vars ; i++)
		v_X_atual[i]	= entrada->valores_ini[i];

	while (/*CONDIÇÃO DE PARADA*/ 1){
		for (int i = 0 ; i < n_vars ; i++){
			v_f_iteracao[i]	= evaluator_evaluate(v_deriv[i], n_vars, v_vars, v_X_atual);
			for (int j = 0 ; j < n_vars ; j++){
				m_f_iteracao[i][j] = evaluator_evaluate(m_deriv[i][j], n_vars, v_vars, v_X_atual);
			}
		}
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