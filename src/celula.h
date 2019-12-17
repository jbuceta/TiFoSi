#pragma once

#include "auxiliar.h"
#include "main.h"

#define MAXIM_ARESTES_C 20                      //Numero maxim d'arestes que pot tenir una celula
#define MAXIM_VERTEXS_C 20                      //Numero maxim de vertexs que pot tenir una celula

class celula
{
public:
  poblacio *p;                            //Punter a la poblacio a la cual pertany;
  int id;                                 //Identificador de la celula
  
  int tipus;                              //Tipus de celula
  
  string c_track_id;
  
  double k_gamma;
  double k_kappa;
  double k_force_x;
  double k_force_y;
  
  double duracio_cicle[MAXIM_FASES];
  double posicio_cicle;
  double pas_cicle;
  double area_growth;
  int fase;
  int allow_change_phase;
  
  double area;
  double area0;
  double kappa_area_area0;                //Guarda el producte k_kappa*(area - area0)
  double perimeter;
  
  double area_centre;                     //Variable interna que s'utilitza per calcular el centre.
  //Correspon a l'area amb signe.
  
  double x;                               //Posicio x i y del centre de masses
  double y;
  
  int nvertexs;                           //Numero de vertexs de la celula
  vertex *v[MAXIM_VERTEXS_C];             //Llista de vertexs en ordre horari
  
  int narestes;                           //Numero d'arestes de la celula
  aresta *a[MAXIM_ARESTES_C];
  
  int ncelules;                           //Numero de celules veines
  celula *c[MAXIM_ARESTES_C];             //Llista de celules veines en l'ordre de les arestes
  
  especies proteines;
  
  //**********************************************************************
  //**********************************************************************
  //Funcions
  //**********************************************************************
  //**********************************************************************
  
  void inicia_celula(poblacio *pt, int i, int tipus_celula, double a0, int nc=0, int nv=0, int na=0);
  void introdueix_vertex(vertex *vi, vertex *v1, vertex *v2);     //Introdueix un vertex en ordre horari,
  //entre els vertexs v1 i v2. Si algun dels dos es
  //el vertex buit. L'afegeix al final de la llista
  //en ordre horari.
  void elimina_vertex(vertex *ve);                        //Elimina el vertex ve;
  //      void elimina_aresta(aresta *ae);                        //Elimina l'aresta ae;
  //      void afegeix_aresta(aresta *aa);                        //Afegeix l'aresta a;
  void troba_arestes();                                   //Torna a trobar les arestes,
  //despres d'un canvi de configuracio. Els
  //vertexs i les seves arestes han d'estar
  //ben definits.
  
  vertex *torna_vertex_anterior(vertex *vc);              //Tornen els vertexs anterior o posterior
  vertex *torna_vertex_posterior(vertex *vc);             //al vertex vc en ordre horari.
  
  void escriu_celula(ofstream &arxiu);
  void escriu_informacio_celula(ofstream &arxiu);
  
  void troba_celules_veines();                            //Troba les celules veines a aquesta celula
  //basant-se en la informacio de les arestes
  
  void calcula_kappa_area_area0();
  void calcula_perimetre();
  void calcula_centre();
  void calcula_area();
  void calcula_area_dinamica();                           //Calcula l'area i mira si la celula s'ha de dividir
  void calcula_rellotge();                                //Avança el rellotge de la celula.
  void calcula_propietats();                              //Calcula algunes propietats de la celula.
  
  void divideix_celula();
  void aux_divideix_celula(celula *cn, vertex *vn1, vertex *vn2, aresta *anc, aresta *an1, aresta *an2);  //Divideix la celula
  //Els punters apunten a la celula, vertexs
  //i aresta nous, respectivament.
  punt calcula_direccio_divisio();
  
  void cambia_referencia_celula(celula *dv, celula *dn);
  
  
  void recalcula_fase();
  void actualitza_constants();
  
        void actualitza_speed_1();
};
