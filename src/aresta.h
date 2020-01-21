#pragma once

#include "auxiliar.h"
#include "main.h"

#define MAXIM_VERTEXS_A 2               //Numero maxim de vertexs que pot tenir una aresta
#define MAXIM_CELULES_A 2               //Numero maxim de celules que pot tenir una aresta

#define TOLERANCE (HEXAGONAL_CELL_EDGE/10.)
#define PASSOS_ABANS_CANVI 1000

class aresta
{
public:
  poblacio *p;                            //Punter a la poblacio a la cual pertany;
  int id;                         //Identificador de l'aresta.
  
  double x[2];                    //Posicio x[1], y[1], x[2] i y[2] que defineixen l'aresta.
  double y[2];
  
  
  
  double lambdah;                 //Constant quant l'aresta es horitzontal
  double lambdav;                 //Constant quant l'aresta es vertical
  double templh2lv;               //(lambdah-2.*lambdav)
  double templv2lh;               //(lambdav-2.*lambdah)
  double k_gamma_aresta;
  
  celula *c[MAXIM_CELULES_A];     //Llista de celules adjacents a l'aresta.
  
  vertex *v[MAXIM_VERTEXS_A];     //Llista de vertexs adjacents a l'aresta.
  celula *cp[MAXIM_CELULES_A];    //Celules que toquen al vertex v[i] i que l'aresta no separa entre si.
  int canvi;                      //Variable que fa de rellotge per evitar que es vagin alternant 
  //processos t1 seguits.
  
  double l;                       //Longitud de l'aresta.
  
  //**********************************************************************
  //**********************************************************************
  //Funcions
  //**********************************************************************
  //**********************************************************************
  
  void escriu_aresta(ofstream &arxiu);
  void escriu_informacio_aresta(ofstream &arxiu);
  
  void calcula_longitud();        //Funcio que calcula l.
  void calcula_longitud_dinamica();//Funcio que calcula l. Si l es menor que un valor, executa un proces t1.
  
  void actualitza_coordenades();  //Funcio que actualitza x1, y1, x2 i y2
  punt comprova_si_intersecta_amb_recta(double alfa, double beta);        //Funcio que comprova si l'aresta
  //intersecta amb la recta que te
  //pendent alfa i ordenada en
  //l'origen beta. Si es aixi, el camp
  //actiu de punt es !=1.
  
  void crea_aresta(poblacio *pt, vertex *v1, vertex *v2, celula *c1, celula *c2);
  void troba_celules_puntes_aresta();
  void calcula_lambda();
  
  void canvia_conexio(vertex *vo, vertex *vv, vertex *vn);        //Cambia la conexio de l'aresta que dels
  //vertexs vo a vv per una conexio de
  //vo a vn. Suposa que els veins de l'aresta
  //no canvien.
  void canvia_celules_veines(celula *cv, celula *cn); //Cambia les celules veines. Mante ce estatica (per comprovar), i canvia cv per cn.
  
  void proces_t1();
  void proces_t2();
  void proces_t3(aresta *acentral);
  
  void cambia_referencia_aresta(aresta *dv, aresta *dn);          //Cambia la direccio de l'aresta i arregla les referencies en les celules i els vertexs.
  matriu troba_moment_inercia(double x, double y);
  void copia_des_de_aresta(aresta *aoriginal);
};
