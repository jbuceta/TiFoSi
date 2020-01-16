#pragma once

#include "auxiliar.h"
#include "main.h"


#define DELTAT 0.001
#define SQRTDELTAT 0.0316227766017

/*Els passos totals de la dinamica es el producte de les constants PASSOS_BUCLE_PARCIAL*PASSOS_BUCLE_PRINCIPAL*/
#define PASSOS_STAGE1 200
#define PASSOS_BUCLE_PARCIAL1 5000

#define N_CELULES 200000
#define N_VERTEXS 800000
#define N_ARESTES 800000

#define NVERTEXS_ARXIU_ESTAT_INICIAL 20

#define MAXIM_NX_INICIAL 24
#define MAXIM_NY_INICIAL 8

#define MAXIM_ELEMENTS_DESTRUIR 10
#define N_TIPUS_CELULA 3
#define N_MAXIM_REGIONS 24

class poblacio
{
public:
  unsigned int idx_interior;
  unsigned int idx_exterior;
  double timepos;
  
  int n_matriu_c;
  int n_matriu_a;
  int n_matriu_v;
  
  celula *matriu_c;
  vertex *matriu_v;
  aresta *matriu_a;
  
  celula celula_buida;
  vertex vertex_buit;
  aresta aresta_buida;
  
  int stage;
  
  int c_track_idx;
  
  int celules_per_destruir[MAXIM_ELEMENTS_DESTRUIR];
  int n_celules_destruir;
  int vertexs_per_destruir[MAXIM_ELEMENTS_DESTRUIR];
  int n_vertexs_destruir;
  int arestes_per_destruir[MAXIM_ELEMENTS_DESTRUIR];
  int n_arestes_destruir;
  
  //double k_kappa;                       //Constant del potencial. Es igual a kappa/2 en l'article.
  double k_gamma_aresta[N_TIPUS_CELULA][N_TIPUS_CELULA];
  double k_force_x[N_TIPUS_CELULA];
  double k_force_y[N_TIPUS_CELULA];
  double k_kappa[N_TIPUS_CELULA];
  double k_gamma[N_TIPUS_CELULA];
  double k_lambda_h[N_TIPUS_CELULA][N_TIPUS_CELULA];
  double k_lambda_v[N_TIPUS_CELULA][N_TIPUS_CELULA];

  double pas_cicle[N_TIPUS_CELULA];
  int nfases[N_TIPUS_CELULA];
  double dispersio_cicle[N_TIPUS_CELULA];
  double area0[N_TIPUS_CELULA][MAXIM_FASES];
  double area0_final[N_TIPUS_CELULA][MAXIM_FASES];
  double reldiv[N_TIPUS_CELULA][MAXIM_FASES];
  double proporcio_cicle[N_TIPUS_CELULA][MAXIM_FASES];
  double divisiondispersion[N_TIPUS_CELULA];
  double divisiondispersionlimit[N_TIPUS_CELULA];
  double divisionshift[N_TIPUS_CELULA];
  
  
  
  
  
  regio posicio_tipus[N_MAXIM_REGIONS];
  int mapa_tipus_celules[MAXIM_NX_INICIAL][MAXIM_NY_INICIAL];
  
  long llavor;
  
  
  
  
  ofstream forces;
  ofstream energia;
  ofstream celules;
  ofstream vertexs;
  ofstream dades;
  ofstream dades_aux;
  ofstream direccio_divisio;
  ofstream flog;
  ofstream ftime;
  ofstream constantspotencial;
  
  clock_t start_time;
  clock_t partial_start_time;
  
  /*********************************************************
   * 
   * Iniciem la definició de funcions.
   *
   *********************************************************/
  
  void bucle_dinamica_principal();
  void bucle_dinamica_parcial_stage_1();
  void actualitza_constants_stage_1();
  void escriu_constants_stage_1();
  void defineix_constants_stage(int s);
  void guarda_dades();
  void guarda_dades_final(std::string stageIdx);
  
  void crea_poblacio();
  void destrueix_poblacio();
  
  void crea_celules_inicials(int ni, int nj);
  void crea_celula_buida_vertex_buit_i_aresta_buida();
  vertex *crea_vertex(celula *c, vertex *vanterior, double cx, double cy);
  void defineix_lambda(double valorh, double valorv, int tipus1, int tipus2);
  
  void fes_vertexs_random();
  void fes_fase_random();
  
  void destrueix_elements();
  
  void crea_mapa_tipus_celules(int nx, int ny);
  
  double function_hill(double x, double k, int n);
  double function_hill_inverse(double x, double k, int n);
  double function_hill_f(double x, double k, double n);
  double funcio_hill_f_inversa(double x, double k, double n);
  double white_noise();
  double f_step(double x, double l);
  double f_step_inversa(double x, double l);
  
  
  void crea_celules_inicials_des_de_arxiu();
};
