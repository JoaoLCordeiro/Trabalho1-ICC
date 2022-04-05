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

int main(){
	t_entrada* 	entrada_atual 	= (t_entrada *) calloc (1,sizeof(t_entrada));

	int 		cont_f = 0;
	int 		num_it_np, num_it_nm, num_it_ni;
	t_tempos	tempos_np, tempos_nm, tempos_ni;

	while (le_entrada(entrada_atual) == 0){
		double*		v_res_np;					//vetor de resultados do newton padrão
		double*		v_res_nm;					//vetor de resultados do newton modificado
		double* 	v_res_ni;					//vetor de resultados do newton inexato

		v_res_np	= newton_padrao 	(entrada_atual, &num_it_np, &tempos_np);
		v_res_nm 	= newton_modificado (entrada_atual, &num_it_nm, &tempos_nm);
		v_res_ni	= newton_inexato	(entrada_atual, &num_it_ni, &tempos_ni);

		imprime_resultados (entrada_atual, v_res_np, v_res_nm, v_res_ni, num_it_np, num_it_nm, num_it_ni, tempos_np, tempos_nm, tempos_ni);

		free (v_res_np);
		free (v_res_nm);
		free (v_res_ni);
		free_entrada(entrada_atual); 			//dar free na entrada para a próxima

		cont_f++;
	}
	
	free (entrada_atual);
	return 1;
}