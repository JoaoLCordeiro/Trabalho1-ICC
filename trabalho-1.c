#include <stdio.h>
#include <stdlib.h>

#define MAX 1000

struct t_entrada{			//tipo entrada: armazena os dados da entrada atual
	int		n_var;			//número de variáveis
	char*	funcao;			//string da funcao (TALVEZ ALTERAR AQUI)
	double*	valores_ini;	//vetor de valores iniciais para as incógnitas
	double	epslon;			//tolerância epslon
	int		iteracoes;		//número máximo de iterações
} typedef t_entrada;

int le_entrada(t_entrada* entrada){					//função que lê a entrada
	scanf("%d\n", &(entrada->n_var));												//lê o número de incógnitas

	entrada->funcao = (char *) malloc ((MAX+1)*sizeof(char));						//malloca e lê a string da funcao
	fgets(entrada->funcao, MAX, stdin);

	entrada->valores_ini = (double *) malloc (entrada->n_var * sizeof(double));		//malloca e lê um vetor de valores iniciais para as incógnitas
	for (int i = 0; i < entrada->n_var ; i++){
		scanf("%lf ",&(entrada->valores_ini[i]));
	}

	entrada->epslon;
	scanf("%lf\n", &(entrada->epslon));												//lê o número para a tolerância epslon

	entrada->iteracoes;
	scanf("%d\n", &(entrada->iteracoes));											//lê o número máximo de iterações
}

void imprime_entrada (t_entrada *entrada){				//funcao só pra debuggar :)
	printf ("%d\n", entrada->n_var);
	printf ("%s", entrada->funcao);
	for (int i = 0; i < entrada->n_var ; i++)
		printf ("%lf ", entrada->valores_ini[i]);
	printf ("\n");
	printf ("%lf\n", entrada->epslon);
	printf ("%d\n", entrada->iteracoes);
}

int main(){
	t_entrada* entrada_atual = (t_entrada *) malloc (sizeof(t_entrada));

	if (le_entrada(entrada_atual) != 0)
		//tratar erros na leitura de funcoes
	
	imprime_entrada (entrada_atual);

	//if (free_funcoes() != 0)
		//tratar erros nos frees das variaveis das funcoes
	return 1;
}