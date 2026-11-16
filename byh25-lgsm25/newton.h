#ifndef NEWTON_H
#define NEWTON_H
 
/*
 * Método de Newton para o SNLB de Broyden.
 *
 * Parâmetros:
 *   n, x0, eps, max : dimensão, valor inicial (igual p/ todas as variáveis),
 *                     tolerância e número máximo de iterações
 *   nIter           : (saída) número de iterações efetivamente realizadas
 *   tTotal          : (saída) tempo total do método (ms)
 *   tJac            : (saída) soma dos tempos de cálculo das Jacobianas (ms)
 *   tSL             : (saída) soma dos tempos de resolução dos sistemas lineares (ms)
 *
 * Retorno: histórico das iterações, matriz (nIter+1) x n em row-major
 *          (linha k = X^(k)). O chamador deve liberar com free().
 *          Retorna NULL em caso de erro de alocação.
 */
double *newton(int n, double x0, double eps, int max, int *nIter,
               double *tTotal, double *tJac, double *tSL);
 
#endif