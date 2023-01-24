# Método Simplex

Resolução de um problema de programação linear de maximização ou minimização através do método Simplex, em linguagem Python 🐍.

Dado uma função objetivo e suas restrições o programa monta o quadro e executa o algoritmo simples ou duas fases, dependendo da entrada.

<div align="right", font-size: 150%><em>Feito por <b>Lucas Guerreiro Pellegrini;</b></em></div>
<div align="right", font-size: 85%><em>Para a disciplina de Computação Científica e Otimização.</em></div>


## Saídas do algoritmo para um caso de maximização através do método simples simples:

 <div align="center", font-size: 150%><em>Exemplo de Execução</em></div>
 
 <pre><code>
lista_restricoes = [[4, 9, 400], [10, 6, 600]]
lista_desigualdades = ['<=', '<=']
z = [6, 8]

p = Problema(z, lista_restricoes, lista_desigualdades, True)
# A diretiva True indica que é maximização

Problema sem normalizar:
MAX Z = +6X0 +8X1
Sujeito a:
+4X0 +9X1  <= 400
+10X0 +6X1  <= 600
Com todos os Xi, fi, ei e ai >= 0

=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

Problema normalizado
MAX Z = +6X0 +8X1
Sujeito a:
+4X0 +9X1 +f0 = 400
+10X0 +6X1 +f1 = 600
Com todos os Xi, fi, ei e ai >= 0

=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

Quadro Simplex
cT = [6, 8, 0, 0]
xT = ['x1', 'x2', 'f0', 'f1']
A = [[4, 9, 1, 0], [10, 6, 0, 1]]
B = [400, 600]
x0 = ['f0', 'f1']
c0 = [0. 0.]
li = [-6.0, -8.0, 0, 0]

=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

Solucao
x2 = 24.24242;
x1 = 45.45455;

=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 </code></pre>
