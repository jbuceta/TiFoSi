#include <iostream>
using namespace std;
#include "main.h"

<especies_cpp_1>

double especies::calcula_r(celula *c1, celula *c2)
{
  double tx, ty;
  
  tx=(c1->x)-(c2->x);
  ty=(c1->y)-(c2->y);
  
  tx=tx*tx;
  ty=ty*ty;
  
  return sqrt(tx+ty);
  
}

double especies::function_hill(double x, double k, int n)
{
  int i;
  double resultat;
  
  resultat=(x/k);
  for(i=1; i<n; i++)
  {
    resultat=resultat*(x/k);
  }
  
  resultat=1./(1.+resultat);
  
  return resultat;
}

double especies::function_hill_inverse(double x, double k, int n)
{
  int i;
  double resultat;
  
  resultat=(x/k);
  for(i=1; i<n; i++)
  {
    resultat=resultat*(x/k);
  }
  
  resultat=resultat/(1.+resultat);
  
  return resultat;
}

double especies::function_hill_f(double x, double k, double n)
{
  double resultat;
  
  resultat=pow((x/k), n);
  
  resultat=1./(1.+resultat);
  
  return resultat;
}

double especies::function_hill_f_inverse(double x, double k, double n)
{
  double resultat;
  
  resultat=pow((x/k), n);
  
  resultat=resultat/(1.+resultat);
  
  return resultat;
}

double especies::white_noise()
{
  return (box_muller(0., 1., &c->p->llavor)/SQRTDELTAT);
}

double especies::f_step(double x, double l)
{
  return (x>l?1.:0.);
}

double especies::f_step_inversa(double x, double l)
{
  return (x>l?0.:1.);
}


<especies_cpp_difusio>

<especies_cpp_senyal>

<especies_cpp_senyal_raw>
