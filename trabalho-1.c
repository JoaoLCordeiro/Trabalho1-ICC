#include <stdio.h>
#include <stdlib.h>

#define MAX 1000

struct t_entrada{
    int     n_var;
    char*   funcao;
    double* valores_ini;
    double  epslon;
    int     iteracoes;
} typedef t_entrada;

int le_funcao(t_entrada* entrada){
    scanf("%d\n", &(entrada->n_var));

    entrada->funcao = (char *) malloc ((MAX+1)*sizeof(char));
    fgets(entrada->funcao, MAX, stdin);

    entrada->valores_ini = (double *) malloc (entrada->n_var * sizeof(double));
    for (int i = 0; i < entrada->n_var ; i++){
        scanf("%lf ",&(entrada->valores_ini[i]));
    }

    entrada->epslon;
    scanf("%lf\n", &(entrada->epslon));

    entrada->iteracoes;
    scanf("%d\n", &(entrada->iteracoes));
}

void imprime_entrada (t_entrada *entrada){
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

    if (le_funcao(entrada_atual) != 0)
        //tratar erros na leitura de funcoes
    
    imprime_entrada (entrada_atual);

    //if (free_funcoes() != 0)
        //tratar erros nos frees das variaveis das funcoes
    return 1;
}