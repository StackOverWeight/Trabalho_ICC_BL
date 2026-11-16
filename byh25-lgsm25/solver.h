#ifndef SOLVER_H
#define SOLVER_H
 
#include "jacobiana.h"
 
/*
 * Resolve o sistema linear tridiagonal J x = b pelo Método da Eliminação de
 * Gauss (sem pivoteamento) seguido de retrossubstituição.
 *
 * ATENÇÃO: J->d e b são modificados (triangularização feita in-place).
 */
void resolveTriDiag(TriDiag *J, double *b, double *x);
 
#endif