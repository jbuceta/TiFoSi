#pragma once

#include "auxiliar.h"
#include "main.h"

#define MAXIM_ARESTES_V 3               //Numero maxim d'arestes que pot tenir un vertex
#define MAXIM_CELULES_V 3               //Numero maxim de celules que pot tenir un vertex

class vertex
{
public:
  poblacio *p;                            //Punter a la poblacio a la cual pertany;
  int id;                                 //Identificador del vertex
  
  double x;                               //Posicio x i y del vertex
  double y;
  
  int ncelules;                           //*Numero de celules adjacents al vertex
  celula *c[MAXIM_CELULES_V];             //*Llista de celules adjacents al vertex
  vertex *vc1[MAXIM_CELULES_V];           //*Punters als vertexs veins a aquest, en ordre horari
  vertex *vc2[MAXIM_CELULES_V];           //i que pertanyen a la célula "i"
  double base[MAXIM_CELULES_V];           //*Longitud de la base del triangle definit pels dos vertexs
  //adjacents i la celula [i];
  
  int narestes;                           //Numero d'arestes adjacents al vertex
  aresta *a[MAXIM_ARESTES_V];             //Llista d'arestes ("i")
  vertex *v[MAXIM_ARESTES_V];             //Llista dels vertexs on van a parar l'aresta "i" de la declaracio 
  //anterior.
  
  double forca_x;
  double forca_y;
  
  double f0_x;                            //Termes per aplicar una força constant al vertex.
  double f0_y;
  
  double fx[4];
  double fy[4];
  
  double energiaf0;
  double energia[4];
  
  //**********************************************************************
  //**********************************************************************
  //Funcions
  //**********************************************************************
  //**********************************************************************
  
  void escriu_vertex();
  void escriu_informacio_vertex(ofstream &arxiu);
  void escriu_informacio_forces(ofstream &arxiu);
  
  void inicia_vertex(poblacio *pt, int i, double cx, double cy);
  
  void calcula_base(int i);               //Calcula la base per la celula "i"
  void calcula_totes_les_bases();         //Calcula totes les bases
  
  void troba_vertexs_veins();             //Busca els vertexs vc1 i vc2
  void canvia_celules_veines(celula *cv, celula *cn);     //canvia la celula cv per la cn. En cas que la cn 
  //sigui una celula buida, l'elimina.
  void elimina_aresta(aresta *aeliminar);         //Elimina l'aresta aeliminar.
  void canvia_conexio_amb_vertex(vertex *vv, aresta *av, vertex *vn, aresta *an); //Canvia una conexio del vertex. Del vertex vv amb l'aresta av al vertex vn amb l'aresta an.
  
  vertex *troba_vertex_celula_celula(celula *c1, celula *c2);     //Troba el vertex que esta unit amb l'actual
  //a traves d'una aresta que te les celules
  //c1 i c2 com a veines.
  
  void calcula_forca();
  void calcula_forces_per_separat();
  void desplaca_vertex();
  void calcula_f0();
  
  void calcula_energia();
  void escriu_informacio_energia(ofstream &arxiu);
  
  void cambia_referencia_vertex(vertex *dv, vertex *dn);
  
  void copia_des_de_vertex(vertex *voriginal);
  aresta* troba_aresta_vertex_oposat(celula *cref, vertex *vref);
};
