#include <iostream>
using namespace std;
#include "main.h"

void especies::calcula_dinamica_especies()
{
}

void especies::actualitza_especies()
{
}

void especies::inicia_especies(celula *ctemp)
{
  c = ctemp;
  if(c->ctype == 2)
  {
  }
  if(c->ctype == 1)
  {
  }
}

void especies::divideix_celula(celula *cfilla)
{
}

void especies::escriu_especies(ofstream &arxiu)
{
  arxiu << 0 << " ";
}

std::string especies::especies_en_ordre()
{
  return "";
}

void especies::llegeix_especies(ifstream &arxiu)
{
  int temp;
  arxiu >> temp;
}

void especies::inicia_constants(int stage)
{
  if(stage == 1)
  {
  }
}


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

double especies::funcio_hill_f_inversa(double x, double k, double n)
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




