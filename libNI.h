//LIBNI: Define as funcoes usadas pelo método de Newton inexato

#ifndef __LIBNI__
#define __LIBNI__

//funcao que inicia o vetor x com seus valores iniciais
void inicia_x (double* v_X, int n, t_entrada* entrada);

//passa as valores de X pro delta
void copia_X_delta (double* v_X, double* v_delta, int n);

//soma delta no valor atual de X
void soma_delta_X (double* v_X, double* v_delta, int n);

//faz o newton inexato
double* newton_inexato (t_entrada* entrada, int* num_it, t_tempos* tempos);

#endif