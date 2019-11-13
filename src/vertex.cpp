#include <iostream>
using namespace std;
#include "main.h"

void vertex::inicia_vertex(poblacio *pt, int i, double cx, double cy)
{
  //El vertex ha d'estar creat en ordre horari per complir l'ordre dintre de la celula
  //on es crea.
  
  p=pt;
  id=i;
  
  x=cx;
  y=cy;
  
  narestes=0;
  ncelules=0;
  
  f0_x=0.;
  f0_y=0.;
  
  
}

void vertex::calcula_base(int i)
{
  double tx,ty;
  
  tx=(vc1[i]->x)-(vc2[i]->x);
  ty=(vc1[i]->y)-(vc2[i]->y);
  
  tx=tx*tx;
  ty=ty*ty;
  
  base[i]=sqrt(tx+ty);
}

void vertex::calcula_totes_les_bases()
{
  int i;
  
  for(i=0;i<ncelules;i++)
  {
    calcula_base(i);
  }
}

void vertex::calcula_forca()
{
  int i;
  double temp;
  double derivada_aresta_x[MAXIM_ARESTES_V];
  double derivada_aresta_y[MAXIM_ARESTES_V];
  double tempx2, tempy2, tempx, tempy;
  double tf1x, tf1y;
  
  
  calcula_f0();
  
  forca_x=f0_x;
  forca_y=f0_y;
  
  //********************************************************************
  //Calcula el terme d'area
  //********************************************************************
  
  
  for(i=0;i<ncelules;i++)                         //Per cada celula, calcula la derivada de l'area.
        {                                               //I suma per totes les celules del vertex.
        temp=(c[i]->kappa_area_area0);
  forca_x=forca_x-(vc2[i]->y-vc1[i]->y)*temp;
  forca_y=forca_y+(vc2[i]->x-vc1[i]->x)*temp;
        }
        fx[0] = forca_x - f0_x;
  fy[0] = forca_y - f0_y;
  
  //********************************************************************
  //Calcula el terme lineal
  //********************************************************************
  
  fx[1] = 0.;
  fy[1] = 0.;
  
  for(i=0;i<narestes;i++)
  {
    derivada_aresta_x[i]=(x-v[i]->x)/(a[i]->l);
    derivada_aresta_y[i]=(y-v[i]->y)/(a[i]->l);
  }
  
  for(i=0;i<narestes;i++)
  {
    
    tempx=(v[i]->x-x);
    tempy=(v[i]->y-y);
    tempx2=tempx*tempx;
    tempy2=tempy*tempy;
    
    tf1x = forca_x;
    tf1y = forca_y;
    forca_x=forca_x+derivada_aresta_x[i]*(tempx2*a[i]->lambdah-tempy2*a[i]->templv2lh)/(tempx2+tempy2);
    forca_y=forca_y+derivada_aresta_y[i]*(tempy2*a[i]->lambdav-tempx2*a[i]->templh2lv)/(tempx2+tempy2);
    fx[1] = fx[1] + forca_x - tf1x;
    fy[1] = fy[1] + forca_y - tf1y;
    
    
    forca_x=forca_x+((x-v[i]->x)*a[i]->k_gamma_aresta);
    forca_y=forca_y+((y-v[i]->y)*a[i]->k_gamma_aresta);
  }
  fx[2] = forca_x - f0_x - fx[0] - fx[1];
  fy[2] = forca_y - f0_y - fy[0] - fy[1];
  
  //********************************************************************
  //Calcula el terme del perimetre
  //********************************************************************
  
  for(i=0;i<narestes;i++)
  {
    temp=((a[i]->c[0]->perimetre*a[i]->c[0]->k_gamma)+(a[i]->c[1]->perimetre*a[i]->c[1]->k_gamma));
    forca_x=forca_x+derivada_aresta_x[i]*temp;
    forca_y=forca_y+derivada_aresta_y[i]*temp;
  }
  fx[3] = forca_x - f0_x - fx[0] - fx[1] - fx[2];
  fy[3] = forca_y - f0_y - fy[0] - fy[1] - fy[2];
  
}

void vertex::desplaca_vertex()
{
  
  
  x=x-forca_x*DELTAT;
  y=y-forca_y*DELTAT;
  
}

void vertex::escriu_vertex()
{
  int i;
  
  std::cout << ncelules << " ";
  
  for(i=0;i<ncelules;i++){
    std::cout << c[i]->x << " " << c[i]->y << " ";
  }
  
  std::cout << x << " " << y << " ";
  
  std::cout << std::endl;
  
}

void vertex::troba_vertexs_veins()
{
  int i,j;
  int temp;
  int flag;
  
  if(id!=-1)
  {
    for(i=0; i<ncelules; i++)
    {
      flag=0;
      temp=c[i]->nvertexs;
      for(j=0; j<temp; j++)
      {
        if(c[i]->v[j]->id==id){
          flag=1;
          vc1[i]=c[i]->v[(j+temp-1)%temp];
          vc2[i]=c[i]->v[(j+1)%temp];
        }
      }
      if(flag==0){
        std::cout << endl << "*******************************" << endl;
        std::cout << endl << "Error: something failed when computing the neigboring vertexes." << endl;
        std::cout << endl << "*******************************" << endl;
        exit(1);
        
        vc1[i]=&p->vertex_buit;
        vc2[i]=&p->vertex_buit;
      }
    }
  }
}

void vertex::escriu_informacio_vertex(ofstream &arxiu)
{
  int i;
  
  arxiu << "Vertex: " << id << std::endl;
  arxiu << "-----------" << std::endl;
  arxiu << "      Numero de celules veines: " << ncelules << std::endl;
  for(i=0; i<ncelules; i++)
  {
    arxiu << "              Celula: " << c[i]->id << std::endl;
    arxiu << "              Amb vertexs veins en ordre horari: " << vc1[i]->id << " i " << vc2[i]->id << std::endl;
    arxiu << "              que formen una base de longitud: " << base[i] << std::endl;
  }
  arxiu << "      Numero d'arestes veines: " << narestes << std::endl;
  for(i=0; i<narestes; i++)
  {
    arxiu << "              Aresta: " << a[i]->id << " que va a parar al vertex: " << v[i]->id << std::endl;
    arxiu << "              que te com a celules veines: " << a[i]->c[0]->id << " i " << a[i]->c[1]->id << std::endl;
    arxiu << "              i que te una longitud de: " << a[i]->l << std::endl;
  }
}


void vertex::escriu_informacio_forces(ofstream &arxiu)
{
  int i;
  
  calcula_forces_per_separat();
  
  arxiu << x << " " << y << " " << -f0_x << " " << -f0_y << " ";
  for(i=0; i<4; i++)
  {
    arxiu << -fx[i] << " " << -fy[i] << " ";
  }
  
  arxiu << -forca_x << " " << -forca_y << std::endl;
}

void vertex::calcula_forces_per_separat()
{
  
  
  
}

void vertex::canvia_celules_veines(celula *cv, celula *cn)
{
  int i;
  int posicio;
  
  if(cv!=cn)
  {
    if(cn!=&p->celula_buida&&cv!=&p->celula_buida)
    {
      for(i=0; i<ncelules; i++)
      {
        if(c[i]==cv)
        {
          c[i]=cn;
        }
      }
    }
    else if(cv==&p->celula_buida)
    {
      c[ncelules]=cn;
      ncelules++;
    }
    else if(cn==&p->celula_buida)
    {
      posicio=-1;
      for(i=0; i<ncelules; i++)
      {
        if(c[i]==cv)
        {
          posicio=i;
        }
      }
      if(posicio==-1)
      {
        std::cout << endl << "*******************************" << endl;
        std::cout << endl << "Error: oh my....I tried to change a cell that actually does not exist at a vertex." << endl;
        std::cout << endl << "*******************************" << endl;
        exit(1);
      }
      else
      {
        posicio++;
        for(i=posicio; i<ncelules; i++)
        {
          c[i-1]=c[i];
          vc1[i-1]=vc1[i];
          vc2[i-1]=vc2[i];
        }
        ncelules--;
      }
    }
  }
  else
  {
    if(cv!=&p->celula_buida)
    {
      std::cout << endl << "*******************************" << endl;
      std::cout << endl << "Error: I tried to change two cells of the same vertex." << endl;
      std::cout << endl << "*******************************" << endl;
      exit(1);
    }
  }
}

void vertex::elimina_aresta(aresta *aeliminar)
{
  int i;
  int posicio;
  
  posicio=-1;
  for(i=0; i<narestes; i++)
  {
    if(a[i]==aeliminar)
    {
      posicio=i;
    }
  }
  if(posicio==-1)
  {
    std::cout << endl << "*******************************" << endl;
    std::cout << endl << "Error: ooooops....I was asked to kill an edge that does not exist." << endl;
    std::cout << endl << "*******************************" << endl;
    exit(1);
  }
  else
  {
    narestes--;
    a[posicio]=a[narestes];
    v[posicio]=v[narestes];
  }
}

void vertex::canvia_conexio_amb_vertex(vertex *vv, aresta *av, vertex *vn, aresta *an)
{
  int i;
  int posicio;
  
  if(vn->id!=-1)
  {
    posicio=-1;
    for(i=0; i<narestes; i++)
    {
      if(v[i]==vv)
      {
        posicio=i;
        v[i]=vn;
        if(a[i]!=av)
        {
          std::cout << endl << "*******************************" << endl;
          std::cout << endl << "Error: something went wrong when updating the information about the vertexes." << endl;
          std::cout << endl << "This is related to changing a connection between vertexes." << endl;
          std::cout << endl << "*******************************" << endl;
          exit(1);
        }
        else
        {
          a[i]=an;
        }
      }
    }
    if(posicio==-1)
    {
      std::cout << endl << "*******************************" << endl;
      std::cout << endl << "Error: something went wrong when updating the information about the vertexes." << endl;
      std::cout << endl << "This is related to changing a connection between vertexes." << endl;
      std::cout << endl << "The old vertex does not exist." << endl;
      std::cout << endl << "*******************************" << endl;
      exit(1);
    }
  }
  else
  {
    posicio=-1;
    for(i=0; i<narestes; i++)
    {
      if(v[i]==vv)
      {
        posicio=i;
      }
    }
    if(posicio!=-1)
    {
      posicio++;
      for(i=posicio; i<narestes; i++)
      {
        v[i-1]=v[i];
        a[i-1]=a[i];
      }
      narestes--;
    }
    else
    {
      std::cout << endl << "*******************************" << endl;
      std::cout << endl << "Error: something went wrong when updating the information about the vertexes." << endl;
      std::cout << endl << "This is related to changing a connection between vertexes." << endl;
      std::cout << endl << "The old vertex does not exist." << endl;
      std::cout << endl << "*******************************" << endl;
      exit(1);
    }
  }
}

vertex *vertex::troba_vertex_celula_celula(celula *c1, celula *c2)
{
  int i;
  int flag;
  vertex *vtornar;
  
  
  
  
  vtornar=&p->vertex_buit;
  flag=0;
  
  for(i=0; i<narestes; i++)
  {
    if((a[i]->c[0]==c1)&&(a[i]->c[1]==c2))
    {
      flag++;
      vtornar=((a[i]->v[0]==this)?a[i]->v[1]:a[i]->v[0]);
    }
    else if((a[i]->c[1]==c1)&&(a[i]->c[0]==c2))
    {
      flag++;
      vtornar=((a[i]->v[0]==this)?a[i]->v[1]:a[i]->v[0]);
    }
  }
  
  if((vtornar==&p->vertex_buit)||(flag!=1))
  {
    std::cout << endl << "*******************************" << endl;
    std::cout << endl << "Error: bummer...I found a problem at the function that look for the vertex connecting cells (cell " << c1->id << ", cell " << c2->id << ")" << endl;
    std::cout << endl << "At vertex " << id << endl;
    std::cout << endl << "it seems that there is no edge with neighboring cells: " << c1->id << " i " << c2->id << endl;
    for(i=0; i<narestes; i++)
    {
      std::cout << endl << "\tThe edge " << a[i]->id << " that connect to vertex " << v[i]->id << endl;
      std::cout << endl << "\thas a neighboring cells " << a[i]->c[0]->id << " and " << a[i]->c[1]->id << endl;
    }
    std::cout << endl << "*******************************" << endl;
    
    exit(1);
  }
  
  
  
  return vtornar;
}

void vertex::cambia_referencia_vertex(vertex *dv, vertex *dn)
{
  int i, j;
  
  if(dv!=dn)
  {
    for(i=0; i<narestes; i++)
    {
      for(j=0; j<2; j++)
      {
        if(a[i]->v[j]==dv)
        {
          a[i]->v[j]=dn;
        }
      }
    }
    
    for(i=0; i<narestes; i++)
    {
      for(j=0; j<v[i]->narestes; j++)
      {
        if(v[i]->v[j]==dv)
        {
          v[i]->v[j]=dn;
        }
      }
      for(j=0; j<v[i]->ncelules; j++)
      {
        if(v[i]->vc1[j]==dv)
        {
          v[i]->vc1[j]=dn;
        }
        if(v[i]->vc2[j]==dv)
        {
          v[i]->vc2[j]=dn;
        }
      }
      
    }
    
    for(i=0; i<ncelules; i++)
    {
      for(j=0; j<c[i]->nvertexs; j++)
      {
        if(c[i]->v[j]==dv)
        {
          c[i]->v[j]=dn;
        }
      }
    }
  }
}

void vertex::calcula_f0()
{
  int i;
  
  
  f0_x = 0.;
  f0_y = 0.;
  
  for(i=0; i<ncelules; i++)
  {
    
    f0_x = f0_x + c[i]->k_force_x;
    f0_y = f0_y + c[i]->k_force_y;
  }
}

void vertex::copia_des_de_vertex(vertex *voriginal)
{
  int i;
  
  narestes = voriginal->narestes;
  for(i=0; i<voriginal->narestes; i++)
  {
    a[i] = voriginal->a[i];
    v[i] = voriginal->v[i];
  }
  ncelules = voriginal->ncelules;
  for(i=0; i<voriginal->ncelules; i++)
  {
    c[i]=voriginal->c[i];
    vc1[i]=voriginal->vc1[i];
    vc2[i]=voriginal->vc2[i];
    base[i]=voriginal->base[i];
  }
}

aresta* vertex::troba_aresta_vertex_oposat(celula *cref, vertex *vref)
{
  int i;
  vertex *vresultat;
  aresta *aresultat;
  
  vresultat=NULL;
  
  for(i=0; i<ncelules; i++)
  {
    if(c[i]==cref)
    {
      if(vc1[i]==vref)
      {
        vresultat=vc2[i];
      }
      else
      {
        vresultat=vc1[i];
      }
    }
  }
  
  if(vresultat==NULL)
  {
    std::cout << endl << "*******************************" << endl;
    std::cout << endl << "Error: oooops...I was expecting an opposing vertex that actually does not exist." << endl;
    std::cout << endl << "*******************************" << endl;
    exit(1);
  }
  
  aresultat = NULL;
  for(i=0; i<narestes; i++)
  {
    if(v[i]==vresultat)
    {
      aresultat = a[i];
    }
  }
  
  if(aresultat==NULL)
  {
    std::cout << endl << "*******************************" << endl;
    std::cout << endl << "Error: oh my...there is a problem between the vertexes v1 and v2 of a vertex and its neighboring vertexes." << endl;
    std::cout << endl << "*******************************" << endl;
    exit(1);
  }
  
  return aresultat;
}

void vertex::calcula_energia()
{
  int i;
  double temp;
    
  energiaf0 = f0_x*x + f0_y*y;


  //********************************************************************
  //Calcula el terme d'area
  //********************************************************************
  temp = 0.;
  energia[0] = 0.;
  for(i=0;i<ncelules;i++)                         
  {                                               
      temp = (c[i]->area-c[i]->area0);
      temp *= temp;
      temp *= c[i]->k_kappa;
      energia[0] += temp;
  }
  
  //********************************************************************
  //Calcula el terme lineal
  //********************************************************************
  
  energia[1] = 0.;
  energia[2] = 0.;
  
  for(i=0;i<narestes;i++)
  {
    energia[1] += a[i]->lambdah*a[i]->l;
    energia[2] += a[i]->k_gamma_aresta*a[i]->l*a[i]->l;
  }
  
  //********************************************************************
  //Calcula el terme del perimetre
  //********************************************************************
  
  energia[3] = 0.;

  for(i=0;i<ncelules;i++)
  {
    energia[3] += c[i]->perimetre*c[i]->k_gamma;
  }
  
}

void vertex::escriu_informacio_energia(ofstream &arxiu)
{
  int i;
  double temptotal;
  
  temptotal = energiaf0;
  arxiu << x << " " << y << " " << energiaf0 << " ";
  for(i=0; i<4; i++)
  {
    arxiu << energia[i] << " ";
    temptotal += energia[i];
  }
  
  arxiu << temptotal << std::endl;
}