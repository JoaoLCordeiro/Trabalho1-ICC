TRABALHO 1 DE ICCC - IMPLEMENTAÃO DOS MÉTODOS DE NEWTON

Alunos:
Gabriel Razzolini Pires De Paula
GRR20197155

João Lucas Cordeiro
GRR20190427

Sobre a organização:

libDefine: 	possui a definição das structs usadas no programa.
libGeral:	implementa as funções usadas na main e por mais de uma implementação do Método de Newton
libNP:		implementa as funções usadas apenas no Newton Padrão
libNM:		implementa as funções usadas apenas no Newton Modificado
libNI:		implementa as funções usadas apenas no Newton Inexato
utils:		implementa as funções que contam o tempo do programa
main:		é o arquivo principal do programa

Sobre a execução:

Na main, allocamos espaço para um t_entrada (t_entrada) e iniciamos um contador de iterações e um t_tempos (tipo tempos) para cada método. O
tipo entrada armazena o número de variáveis, a função, os valores de x iniciais, o erro epslon e o número máximo de iterações para cada
entrada. O tipo tempo armazena o tempo total, das derivadas e do sistema linear de cada entrada, usamos um para cada método.