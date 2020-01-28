#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>

class especies;
class poblacio;
class celula;
class vertex;
class aresta;

class regio
{
public:
  int tipus;
  int x;
  int y;
  int amplada;
  int alcada;
  
  int prova_si_esta_dins(int px, int py)
  {
    int resultat;
    
    if(x==-1||y==-1)
    {
      resultat=0;
    }
    else
    {
      if((px>=x)&&(px<(x+amplada)))
      {
        if((py>=y)&&(py<(y+alcada)))
        {
          resultat=1;
        }
        else
        {
          resultat=0;
        }
      }
      else
      {
        resultat=0;
      }
    }
    
    return resultat;
  }
};

struct punt
{
  int actiu;
  double x;
  double y;
};

#define TAMANY_MAXIM_SERIE 10

class serie
{
public:
  int desordenat[TAMANY_MAXIM_SERIE];
  int ordenat_de_major_a_menor[TAMANY_MAXIM_SERIE];
  int ordenat_de_menor_a_major[TAMANY_MAXIM_SERIE];
  int n_elements;
  
  void ordena()
  {
    int i, j;
    int temp;
    
    for(i=0; i<n_elements; i++)
    {
      temp=desordenat[i];
      for(j=i+1; j<n_elements; j++)
      {
        if(temp<desordenat[j])
        {
          desordenat[i]=desordenat[j];
          desordenat[j]=temp;
          temp=desordenat[i];
        }
      }
    }
    
    for(i=0; i<n_elements; i++)
    {
      ordenat_de_major_a_menor[i]=desordenat[i];
    }
    
    for(i=0; i<n_elements; i++)
    {
      ordenat_de_menor_a_major[n_elements-i-1]=desordenat[i];
    }
  }
};

class matriu
{
public:
  double m[2][2];
  
  double lambda1;
  double lambda2;
  punt v1;
  punt v2;
  
  void troba_valors_i_vectors_propis()
  {
    double temp;
    
    if(m[0][1]==m[1][0])
    {
      temp=sqrt((m[0][0]*m[0][0]) + 4.*m[0][1]*m[1][0] - 2.*m[0][0]*m[1][1] + (m[1][1]*m[1][1]));
      
      lambda1=(-temp + m[0][0] + m[1][1])/2.;
      lambda2=(temp + m[0][0] + m[1][1])/2.;
      
      v1.x=-(temp - m[0][0] + m[1][1])/(2.*m[1][0]);
      v1.y=1.;
      
      v2.x=-(-temp - m[0][0] + m[1][1])/(2.*m[1][0]);
      v2.y=1.;
    }
    else
    {
      std::cout << "S'ha produit un error al buscar els valors i vectors" << std::endl;
      std::cout << "propis de la matriu:" << std::endl << std::endl;
      std::cout << "(\t" << m[0][0] << "\t" << m[0][1] << "\t)" << std::endl;
      std::cout << "(\t" << m[1][0] << "\t" << m[1][1] << "\t)" << std::endl;
      exit(1);
    }
  }
  
  void ordena_valors_i_vectors_propis()
  {
    double temp;
    
    if(fabs(lambda1)<fabs(lambda2))
    {
      temp=lambda1;
      lambda1=lambda2;
      lambda2=temp;
      
      temp=v1.x;
      v1.x=v2.x;
      v2.x=temp;
      
      temp=v1.y;
      v1.y=v2.y;
      v2.y=temp;
    }
  }
  
  matriu operator + (matriu b)
  {
    matriu r;
    
    r.m[0][0] = m[0][0]+b.m[0][0];
    r.m[0][1] = m[0][1]+b.m[0][1];
    r.m[1][0] = m[1][0]+b.m[1][0];
    r.m[1][1] = m[1][1]+b.m[1][1];
    
    return r;
  }
};
