#include <float.h>
#include <math.h>
#include <stddef.h>
#include "solver.h"

/*
 * resolveTriDiag: resolve o sistema tridiagonal pelo algoritmo de Thomas
 * (eliminação de Gauss especializada), sem pivoteamento.
 * Antes de cada divisão, verifica se o pivô é suficientemente diferente
 * de zero em relação à escala dos coeficientes da matriz.
 */
int resolveTriDiag(TriDiag *J, double *b, double *x)
{
    const int n = J->n;
    double *a = J->a;
    double *d = J->d;
    double *c = J->c;
    double escala = 0.0;

    for (int i = 0; i < n; ++i) {
        if (!isfinite(d[i]) || !isfinite(b[i])) return -1;
        if (fabs(d[i]) > escala) escala = fabs(d[i]);
    }
    for (int i = 0; i < n - 1; ++i) {
        if (!isfinite(a[i]) || !isfinite(c[i])) return -1;
        if (fabs(a[i]) > escala) escala = fabs(a[i]);
        if (fabs(c[i]) > escala) escala = fabs(c[i]);
    }

    /* Matriz nula ou pivô indistinguível de zero em precisão double. */
    const double limite = 16.0 * DBL_EPSILON * escala;

    /* Eliminação para triangularizar a matriz. */
    for (int i = 0; i < n - 1; ++i) {
        if (fabs(d[i]) <= limite) return -1;
        const double m = a[i] / d[i];
        d[i + 1] -= c[i] * m;
        b[i + 1] -= b[i] * m;
        if (!isfinite(d[i + 1]) || !isfinite(b[i + 1])) return -1;
    }

    /* Retrossubstituição. */
    if (fabs(d[n - 1]) <= limite) return -1;
    x[n - 1] = b[n - 1] / d[n - 1];
    if (!isfinite(x[n - 1])) return -1;

    for (int i = n - 2; i >= 0; --i) {
        if (fabs(d[i]) <= limite) return -1;
        x[i] = (b[i] - c[i] * x[i + 1]) / d[i];
        if (!isfinite(x[i])) return -1;
    }
    return 0;
}
