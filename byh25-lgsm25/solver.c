#include "solver.h"
 
/*
 * resolveTriDiag: Eliminação de Gauss para matrizes tridiagonais.
 *  - Triangularização: elimina a subdiagonal (apenas d[i+1] e b[i+1] mudam).
 *  - Retrossubstituição: só a superdiagonal c está acima da diagonal.
 */
void resolveTriDiag(TriDiag *J, double *b, double *x)
{
    int n = J->n;
    double *a = J->a, *d = J->d, *c = J->c;
 
    /* Triangularização */
    for (int i = 0; i < n - 1; ++i) {
        double m = a[i] / d[i];
        a[i] = 0.0;
        d[i + 1] -= c[i] * m;
        b[i + 1] -= b[i] * m;
    }
 
    /* Retrossubstituição */
    x[n - 1] = b[n - 1] / d[n - 1];
    for (int i = n - 2; i >= 0; --i)
        x[i] = (b[i] - c[i] * x[i + 1]) / d[i];
}