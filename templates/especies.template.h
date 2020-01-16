#pragma once

#include "auxiliar.h"
#include "main.h"

#define ESCALA 0.1

class especies
{
public:
  celula *c;
  
  
  
  <especies_h_constants>
  
public:
  void calcula_dinamica_especies();
  void actualitza_especies();
  void inicia_especies(celula *ctemp);
  void inicia_constants(int stage);
  void divideix_celula(celula *cfilla);
  void escriu_especies(ofstream &arxiu);
  std::string especies_en_ordre();
  void llegeix_especies(ifstream &arxiu);
  
  double calcula_r(celula *c1, celula *c2);
  double function_hill(double x, double k, int n);
  double function_hill_inverse(double x, double k, int n);
  double funcio_hill_f(double x, double k, double n);
  double funcio_hill_f_inversa(double x, double k, double n);
  double white_noise();
  double f_step(double x, double l);
  double f_step_inversa(double x, double l);
  
  
  <especies_h_funcions>
};
