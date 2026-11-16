#include <stdlib.h>
#include <math.h>
#include "newton.h"
#include "funcao.h"
#include "jacobiana.h"
#include "solver.h"
#include "utils.h"
 
/* normaMax: norma do máximo (infinito) de um vetor. */
static double normaMax(const double *v, int n)
{
    double m = 0.0;
    for (int i = 0; i < n; ++i) {
        double a = fabs(v[i]);
        if (a > m) m = a;
    }
    return m;
}
 
/*
 * newton: Método de Newton  J(X^(i)) Delta^(i) = -F(X^(i)),  X^(i+1) = X^(i) + Delta^(i).
 * Critérios de parada: ||F(X^(i))|| < eps, ||Delta^(i)|| < eps ou i = max.
 * O histórico de iterações é guardado em memória (cresce dinamicamente) para
 * que a impressão fique FORA das medições de tempo e dos marcadores LIKWID.
 */
double *newton(int n, double x0, double eps, int max, int *nIter,
               double *tTotal, double *tJac, double *tSL)
{
    TriDiag J;
    double *F = malloc(n * sizeof(double));
    double *b = malloc(n * sizeof(double));
    double *delta = malloc(n * sizeof(double));
    int cap = (max < 64 ? max : 64) + 1;          /* linhas do histórico */
    double *hist = malloc((size_t)cap * n * sizeof(double));
 
    if (!F || !b || !delta || !hist || alocaTriDiag(&J, n) != 0) {
        free(F); free(b); free(delta); free(hist);
        return NULL;
    }
 
    for (int j = 0; j < n; ++j) hist[j] = x0;
 
    int k = 0;
    double tj = 0.0, ts = 0.0, t, t0;
 
    LIKWID_MARKER_START("NEWTON");
    t0 = timestamp();
 
    for (int i = 0; i < max; ++i) {
        /* garante espaço para a linha i+1 */
        if (i + 2 > cap) {
            cap *= 2;
            double *tmp = realloc(hist, (size_t)cap * n * sizeof(double));
            if (!tmp) { k = -1; break; }
            hist = tmp;
        }
        double *x = hist + (size_t)i * n;
 
        calculaF(x, F, n);
        if (normaMax(F, n) < eps) break;
 
        /* Jacobiana */
        t = timestamp();
        LIKWID_MARKER_START("JACOBIANA");
        calculaJacobiana(x, &J);
        LIKWID_MARKER_STOP("JACOBIANA");
        tj += timestamp() - t;
 
        for (int j = 0; j < n; ++j) b[j] = -F[j];
 
        /* Sistema linear */
        t = timestamp();
        LIKWID_MARKER_START("SISTEMA_LINEAR");
        resolveTriDiag(&J, b, delta);
        LIKWID_MARKER_STOP("SISTEMA_LINEAR");
        ts += timestamp() - t;
 
        double *xn = hist + (size_t)(i + 1) * n;
        for (int j = 0; j < n; ++j) xn[j] = x[j] + delta[j];
        k = i + 1;
 
        if (normaMax(delta, n) < eps) break;
    }
 
    *tTotal = timestamp() - t0;
    LIKWID_MARKER_STOP("NEWTON");
 
    *tJac = tj;
    *tSL = ts;
    *nIter = k;
 
    free(F); free(b); free(delta);
    liberaTriDiag(&J);
 
    if (k < 0) { free(hist); return NULL; }
    return hist;
}
 