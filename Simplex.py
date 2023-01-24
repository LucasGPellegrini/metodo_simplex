import numpy as np


# Programa para resolver Simplex através do Método Simples ou Duas Fases
# Por Lucas Guerreiro Pellegrini

# classe para representar o problema
# função objetivo Z na forma de lista de coeficientes
# lista de inequações
# lista de desigualdades para cada inequação
# flag dizendo se o problema é de maximização (default ou True) ou minimização (False).
# Uma lista para cada tipo de variável extra (folga, excesso e artificial), cada uma de tamanho n, onde n é
#   a quantidade de inequações de restrição.
class Problema:
    def __init__(self, z, lista_rest, lista_desig, is_max=True):
        self.fc_obj = z
        self.restricoes = lista_rest
        self.desigualdades = lista_desig
        self.folga = np.zeros(len(lista_rest))
        self.excesso = np.zeros(len(lista_rest))
        self.artificial = np.zeros(len(lista_rest))
        self.is_normal = False
        self.tipo = is_max

    # Função para imprimir na tela os dados do problema
    # Usada para os testes de normalização (TESTES NO FINAL DO ARQUIVO)
    def printa_problema(self):
        z = "MAX Z = " if self.tipo else "MIN Z = "
        rest_list = []

        for i, coef in enumerate(self.fc_obj):
            if coef >= 0:
                z += "+"
            z += str(coef)
            z += "X" + str(i) + " "

        sinal_m = " -" if self.tipo else " +"
        for j, a in enumerate(self.artificial):
            if a == 1:
                z += sinal_m
                z += "Ma" + str(j)

        # Restrições:
        for k, rest in enumerate(self.restricoes):
            temp_r = ""
            igualdade = ""
            for idx, cf in enumerate(rest):
                if idx != len(rest) - 1:
                    if cf >= 0:
                        temp_r += "+"
                    temp_r += str(cf)
                    temp_r += "X" + str(idx) + " "
                else:
                    igualdade = str(cf)

            if self.is_normal == True:
                if self.desigualdades[k] == '<' or self.desigualdades[k] == '<=':
                    temp_r += "+f" + str(k)
                elif self.desigualdades[k] == '>' or self.desigualdades[k] == '>=':
                    temp_r += "-e" + str(k) + " +a" + str(k)
                elif self.desigualdades[k] == '=':
                        temp_r += "+a" + str(k)
            else:
                temp_r += " " + self.desigualdades[k] + " "

            if self.is_normal:
                temp_r += " = "
            temp_r += igualdade

            rest_list.append(temp_r)

        print(z)
        print("Sujeito a:")
        for res in rest_list:
            print(res)
        print("Com todos os Xi, fi, ei e ai >= 0")

    # Função para normalizar o problema,
    # É a primeira função executada ao tentar resolver o problema.
    # Testes de normalização no final do arquivo.
    def normalizar(self):
        for i, desig in enumerate(self.desigualdades):
            if desig == '<' or desig == '<=':
                self.folga[i] = 1
            elif desig == '>' or desig == '>=':
                self.excesso[i] = 1
                self.artificial[i] = 1
            elif desig == '=':
                self.artificial[i] = 1
            else:
                print("Atenção, a entrada está inconsistente, verifique os dados")
                return

        self.is_normal = True

    # Função generalista para resolver o problema de Programação Linear.
    # Primeiro normaliza o problema;
    # Depois resolve o problema através do método simples ou duas fases;
    #   dependendo se existem ou não variáveis artificiais.
    #   Nesta etapa chamam-se duas funções, uma que monta o quadro simplex,
    #       e outra que resolve.
    def resolver(self):
        self.normalizar()
        if any(self.artificial):
            cT, xT, A, B, x0, c0, li, liM = self.monta_simplex_duas_fases()
            return self.resolve_simplex_duas_fases(cT, xT, A, B, x0, c0, li, liM)
        else:
            cT, xT, A, B, x0, c0, lin_inf = self.monta_simplex_simples()
            return self.resolve_simplex_simples(cT, xT, A, B, x0, c0, lin_inf)

    # Função para montar o quadro simplex Simples.
    # Retorna todos os dados de um quadro simplex, possibilitando resolver o problema.
    def monta_simplex_simples(self):
        cT = self.fc_obj.copy()
        x0 = []
        c0 = np.zeros(len(self.restricoes))
        B = []
        xT = []
        A = []

        k = 1
        for i in self.restricoes:
            B.append(i[-1])
            linha_A = i.copy()
            linha_A.pop()
            A.append(linha_A)

        for _ in self.restricoes[0][:-1]:
            xT.append("x" + str(k))
            k += 1

        for ind, val in enumerate(self.folga):
            if val != 0:
                cT.append(0)
                xT.append("f"+str(ind))
                x0.append("f" + str(ind))
                for i, item in enumerate(A):
                    item.append(1) if i == ind else item.append(0)

        linha_inferior = []

        for ind, val in enumerate(cT):
            x = 0
            if ind < len(A):
                for equacao in A:
                    x += c0[ind]*equacao[ind]
            x -= val
            linha_inferior.append(x)

        return cT, xT, A, B, x0, c0, linha_inferior

    # Função para resolver o quadro simplex Simples.
    # Recebe todos os dados da função que monta o quadro simplex;
    # Resolve através do algoritmo descrito no material disponível no moodle.
    def resolve_simplex_simples(self, cT, xT, A, B, x0, c0, lin_inf):
        # Localizar Coluna de Trabalho:
        coluna_trab = min(range(len(lin_inf)), key=lin_inf.__getitem__)

        # Localizar Pivô:
        potenciais_pivos = []
        for ind, val in enumerate(B):
            if A[ind][coluna_trab] <= 0:
                potenciais_pivos.append(99999)
            else:
                potenciais_pivos.append(val / A[ind][coluna_trab])

        pivo = min(range(len(potenciais_pivos)), key=potenciais_pivos.__getitem__)

        # Trocar variaveis:
        var_troca = xT[coluna_trab]
        xT[coluna_trab] = x0[pivo]
        x0[pivo] = var_troca

        # Igualar o Pivo = 1
        B[pivo] = B[pivo] / A[pivo][coluna_trab]
        A[pivo] = [x / A[pivo][coluna_trab] for x in A[pivo]]

        # Zerar Coluna de Trabalho:
        # ---   Linhas de A:    ---
        for ind, val in enumerate(A):
            if ind != pivo:
                A[ind] = [x - (val[coluna_trab]*A[pivo][idx]) for idx, x in enumerate(val)]
                B[ind] = B[ind] - val[coluna_trab]*B[pivo]
        # ---   Linha Inferior: ---
        lin_inf = [x - (lin_inf[coluna_trab] * A[pivo][idx]) for idx, x in enumerate(lin_inf)]

        resposta = [(x, B[ind]) for ind, x in enumerate(x0)]

        # Chamada recursiva caso ainda existam elementos negativos na linha inferior.
        return self.resolve_simplex_simples(cT, xT, A, B, x0, c0, lin_inf) if any(x < 0 for x in lin_inf) else resposta

    # Função para montar o quadro simplex Duas Fases.
    # Retorna todos os dados de um quadro simplex duas fases, possibilitando resolver o problema.
    def monta_simplex_duas_fases(self):
        cT = self.fc_obj.copy()
        x0 = []
        c0 = ["M" for _ in self.restricoes]
        B = []
        xT = []
        A = []

        k = 1
        for i in self.restricoes:
            B.append(i[-1])
            linha_A = i.copy()
            linha_A.pop()
            A.append(linha_A)

        for _ in self.restricoes[0][:-1]:
            xT.append("x" + str(k))
            k += 1

        for ind, val in enumerate(self.folga):
            if val != 0:
                cT.append(0)
                xT.append("f" + str(ind))
                x0.append("f" + str(ind))
                for i, item in enumerate(A):
                    item.append(1) if i == ind else item.append(0)

        for ind, val in enumerate(self.excesso):
            if val != 0:
                cT.append(0)
                xT.append("e" + str(ind))
                for i, item in enumerate(A):
                    item.append(-1) if i == ind else item.append(0)

        for ind, val in enumerate(self.artificial):
            if val != 0:
                cT.append(0)
                xT.append("a" + str(ind))
                x0.append("a" + str(ind))
                for i, item in enumerate(A):
                    item.append(1) if i == ind else item.append(0)

        linha_inferior = cT.copy()
        linha_inferior_M = []

        for ind, var in enumerate(xT):
            valor_m = 0
            if var[0] != 'a':
                for elem in A:
                    valor_m += elem[ind]
            linha_inferior_M.append(-valor_m)

        return cT, xT, A, B, x0, c0, linha_inferior, linha_inferior_M

    # Função para resolver o quadro simplex Duas Fases.
    # Recebe todos os dados da função que monta o quadro simplex;
    # Resolve através do algoritmo descrito no material disponível no moodle.
    def resolve_simplex_duas_fases(self, cT, xT, A, B, x0, c0, li, liM, fase2=False):
        # ------------------------------------------ Caso Segunda Fase: ------------------------------------------
        if fase2:
            coluna_trab = -1
            # Localizar Coluna de Trabalho:
            for ind, val in enumerate(liM):
                if val == 0:
                    if li[ind] < 0:
                        coluna_trab = ind
                        break
            if coluna_trab == -1:
                return [(x, B[ind]) for ind, x in enumerate(x0)]

            # Localizar Pivô:
            potenciais_pivos = []
            for ind, val in enumerate(B):
                if A[ind][coluna_trab] <= 0:
                    potenciais_pivos.append(99999)
                else:
                    potenciais_pivos.append(val / A[ind][coluna_trab])

            pivo = min(range(len(potenciais_pivos)), key=potenciais_pivos.__getitem__)

            # Trocar variaveis:
            var_troca = xT[coluna_trab]
            xT[coluna_trab] = x0[pivo]
            x0[pivo] = var_troca

            # Igualar o Pivo = 1
            B[pivo] = B[pivo] / A[pivo][coluna_trab]
            A[pivo] = [x / A[pivo][coluna_trab] for x in A[pivo]]

            # Zerar Coluna de Trabalho:
            # ---   Linhas de A:    ---
            for ind, val in enumerate(A):
                if ind != pivo:
                    A[ind] = [x - (val[coluna_trab] * A[pivo][idx]) for idx, x in enumerate(val)]
                    B[ind] = B[ind] - val[coluna_trab] * B[pivo]
            # ---   Linha Inferior: ---
            li = [x - (li[coluna_trab] * A[pivo][idx]) for idx, x in enumerate(li)]

            resposta = [(x, B[ind]) for ind, x in enumerate(x0)]

            return self.resolve_simplex_duas_fases(cT, xT, A, B, x0, c0, li, liM, True)
        # ------------------------------------------ Fim do caso Segunda Fase ------------------------------------------


        # Localizar Coluna de Trabalho:
        coluna_trab = min(range(len(liM)), key=liM.__getitem__)

        # Localizar Pivô:
        potenciais_pivos = []
        for ind, val in enumerate(B):
            if A[ind][coluna_trab] <= 0:
                potenciais_pivos.append(99999)
            else:
                potenciais_pivos.append(val / A[ind][coluna_trab])

        pivo = min(range(len(potenciais_pivos)), key=potenciais_pivos.__getitem__)

        # Trocar variaveis:
        var_troca = xT[coluna_trab]
        if x0[pivo][0] == 'a':
            artf = xT.index(x0[pivo])
            cT.pop(artf)
            xT.pop(artf)
            li.pop(artf)
            liM.pop(artf)
            for linha in A:
                linha.pop(artf)
        xT[coluna_trab] = x0[pivo]
        x0[pivo] = var_troca

        # Igualar o Pivo = 1
        B[pivo] = B[pivo] / A[pivo][coluna_trab]
        A[pivo] = [x / A[pivo][coluna_trab] for x in A[pivo]]

        # Zerar Coluna de Trabalho:
        # ---   Linhas de A:    ---
        for ind, val in enumerate(A):
            if ind != pivo:
                A[ind] = [x - (val[coluna_trab] * A[pivo][idx]) for idx, x in enumerate(val)]
                B[ind] = B[ind] - val[coluna_trab] * B[pivo]
        # ---   Linhas Inferiores: ---
        li = [x - (li[coluna_trab] * A[pivo][idx]) for idx, x in enumerate(li)]
        liM = [x - (liM[coluna_trab] * A[pivo][idx]) for idx, x in enumerate(liM)]

        # Chamada recursiva da função.
        #   Chama para a primeira fase, caso existam elementos negativos na linha inferior de M;
        #   Caso contrário, chama para a segunda fase, através da diretiva True.
        return self.resolve_simplex_duas_fases(cT, xT, A, B, x0, c0, li, liM) if any(x < 0 for x in liM) \
            else self.resolve_simplex_duas_fases(cT, xT, A, B, x0, c0, li, liM, True)



def main():
    lista_restricoes = [[4, 9, 400], [10, 6, 600]]
    lista_desigualdades = ['<=', '<=']
    z = [6, 8]
    p = Problema(z, lista_restricoes, lista_desigualdades, True)

    print("Solucao")
    solucao = p.resolver()
    for var, val in solucao:
        print(f'{var} = {val:.5f};')
    print("\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")


if __name__ == '__main__':
    main()

# -------------------------------------  TESTES ------------------------------------- #
"""
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

"""

# --------------------------------  TESTE DUAS FASES -------------------------------- #
"""
lista_restricoes = [[0.1, 0, 0.4], [0, 0.1, 0.6], [0.1, 0.2, 2], [0.2, 0.1, 1.7]]
lista_desigualdades = ['>=', '>=','>=','>=']
z = [80, 32]

p = Problema(z, lista_restricoes, lista_desigualdades, False)
# A diretiva False indica que é minimização

Problema sem normalizar: 
MAX Z = +80X0 +32X1 
Sujeito a:
+0.1X0 +0X1  >= 0.4
+0X0 +0.1X1  >= 0.6
+0.1X0 +0.2X1  >= 2
+0.2X0 +0.1X1  >= 1.7
Com todos os Xi, fi, ei e ai >= 0

=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

Problema normalizado
MAX Z = +80X0 +32X1  -Ma0 -Ma1 -Ma2 -Ma3
Sujeito a:
+0.1X0 +0X1 -e0 +a0 = 0.4
+0X0 +0.1X1 -e1 +a1 = 0.6
+0.1X0 +0.2X1 -e2 +a2 = 2
+0.2X0 +0.1X1 -e3 +a3 = 1.7
Com todos os Xi, fi, ei e ai >= 0

=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

Quadro Simplex
cT = [80, 32, 0, 0, 0, 0, 0, 0, 0, 0]
xT = ['x1', 'x2', 'e0', 'e1', 'e2', 'e3', 'a0', 'a1', 'a2', 'a3']
A = [[0.1, 0, -1, 0, 0, 0, 1, 0, 0, 0], [0, 0.1, 0, -1, 0, 0, 0, 1, 0, 0], [0.1, 0.2, 0, 0, -1, 0, 0, 0, 1, 0], [0.2, 0.1, 0, 0, 0, -1, 0, 0, 0, 1]]
B = [0.4, 0.6, 2, 1.7]
x0 = ['a0', 'a1', 'a2', 'a3']
c0 = ['M', 'M', 'M', 'M']
li = [80, 32, 0, 0, 0, 0, 0, 0, 0, 0]
liM = [-0.4, -0.4, 1, 1, 1, 1, 0, 0, 0, 0]

=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

Solucao
x1 = 8.00000;
x2 = 6.00000;
e3 = 0.50000;
e0 = 0.40000;

=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
"""

# --------------------------------------  FIM! -------------------------------------- #
