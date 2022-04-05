//LIBDEFINE: Define as structs usadas

#ifndef __LIBDEF__
#define __LIBDEF__

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

struct t_tempos{			//tipo que guarda os três tempos de cada funcao
	double tmp_total;
	double tmp_deriv;
	double tmp_SL;
} typedef t_tempos;

#endif