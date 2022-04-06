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

	Na main, allocamos espaço para um t_entrada (t_entrada) e iniciamos um contador de iterações e um t_tempos (tipo tempos) para cada
método. O tipo entrada armazena o número de variáveis, a função, os valores de x iniciais, o erro epslon e o número máximo de iterações para cada entrada. O tipo tempo armazena o tempo total, das derivadas e do sistema linear de cada entrada, então usamos um para cada método.
	Dentro do laço que passa pelas entradas, declaramos vetores que armazenarão os resultados de cada iteração para cada método. Então
chamamos cada método implementado para a entrada atual. Depois disso, temos os resultados e os tempos de cada método, então os imprimimos, 
e finalizamos o laço dando free no que alocamos.
	Acabando o laço, damos free na variável que guardava a entrada atual e finalizamos a execução.

	No Método de Newton de padrão, começamos adquirindo o número de variáveis e seus nomes, então allocamos espaço para guardar as funções
derivadas primeiras em um vetor (vetor gradiente) e derivadas segundas em uma matriz (matriz hessiana). Também allocamos espaço para um vetor
e uma matriz de mesmo tamanho das anteriores, elas irão guardar os resultados das derivadas a cada iteração. Allocamos então espaço para um
vetor para guardar os valores atuais de X, iniciando-os com os valores iniciais, e um vetor de double com índice (t_i_double) onde
guardaremos o delta. O uso do double com índice vem da necessidade de gravar em algum lugar as trocas de linhas que ocorrem mais pra frente. 
Terminamos allocando espaço para o vetor que guardará os resultados de cada iteração e iniciando um auxiliador para contagem de tempo e um
contador de iterações.
	Começamos então um laço, aumentamos o número de iterações e calculamos os valores das derivadas dados os X atuais. Verificamos se a norma
do vetor gradiente já é menor que o erro epslon, caso sim, quebramos o laço. Logo depois disso, resolvemos o sistema linear matriz hessiana
resolvida "vezes" delta "igual" vetor gradiente. Depois disso, reordenamos o delta já que ele foi trocado de lugar no pivoteamento que
acontece na resolução do sistema linear, ficando na ordem original. Então somamos o delta no X atual e temos o X da próxima iteração e
armazenamos o resultado atual no vetor de resultados. Verificamos se esse resultado é um NAN ou um INF, se sim, quebramos o laço. Outra
verificação é a norma do vetor delta, se for menos que o erro epslon, quebramos o laço. O laço também quebra quando passamos o máximo de
iterações.
	Quando terminamos o laço, damos free nas variáveis necessárias, alteramos o valor do contador de iterações do método e retornamos o
endereço do vetor de resultados, encerrando a execução deste método.

	Já no Método de Newton Modificado, o começo com allocações e iniciações é bem parecido, a diferença é a allocação de um vetor que irá
armazenar o resultado intermediário da fatoração LU, o vetor Y, e a de uma matriz a mais, já que agora temos a L e a U. Também temos uma
variável para guardar o "hess_steps", para sabermos quando precisamos re-calcular os valores da matriz hessiana.
	No laço, a verificação de quebras de laço são as mesmas. Começamos calculando os valores do vetor gradiente. Quando estamos em uma
iteração módulo de hess_steps, calculamos os valores da matriz hessiana, reordenamos o delta e então fazemos a fatoração LU com
pivoteamento. Depois disso, trocamos as posições dos elementos dos valores do vetor gradiente de acordo com as trocas que aconteceram na
fatoração LU. Então fazemos retrossubstituição resolvendo Ly = b e depois Ux = y. Então, somamos delta nos X atuais, sem re-ordenar o
delta para ele manter as informações de troca de linhas até a próxima hess_steps. Terminando o laço, armazenamos o resultado atual no
vetor de resultados.
	O fim também é parecido, damos free nas variáveis necessárias, alteramos o número de iterações do método e retornamos o vetor
de resultados.

	Terminando com o Método de Newton Inexato, também temos um início com parecido com o padrão, a diferença está na não iniciação
do vetor delta, já que ele é zerado no começo das iterações do Gauss-Seidel.
	O escopo do laço também possui semelhanças: os critérios de quebra de laço só possuem uma diferença, não verificamos mais a norma do
vetor gradiente. Também carregamos o vetor resultados com o resultado atual no fim do laço. Voltando ao começo, zeramos o delta, como
comentado anteriormente. Depois, usamos o método de Gauss-Seidel para calcular o delta atual. Após calculado, somamos no X atual, 
conseguindo o X da iteração futura.
	O fim mantém a semelhança, com os free nas variáveis, alteração no número de iterações do método e o retorno do vetor de resultados.