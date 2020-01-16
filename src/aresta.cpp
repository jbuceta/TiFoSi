#include <iostream>
using namespace std;
#include "main.h"

void aresta::calcula_longitud()
{
  double tx, ty;
  
  tx=(v[0]->x)-(v[1]->x);
  ty=(v[0]->y)-(v[1]->y);
  
  tx=tx*tx;
  ty=ty*ty;
  
  l=sqrt(tx+ty);
  
}

void aresta::calcula_longitud_dinamica()
{
  int i;
  int index, flag;
  double tx, ty;
  bool comp1, comp2;
  aresta *atemp;
  
  atemp=NULL;
  
  tx=(v[0]->x)-(v[1]->x);
  ty=(v[0]->y)-(v[1]->y);
  
  tx=tx*tx;
  ty=ty*ty;
  
  l=sqrt(tx+ty);
  
  canvi++;
  
  if((l<TOLERANCIA)&&(canvi>PASSOS_ABANS_CANVI))
  {
    canvi=0;
    
    troba_celules_puntes_aresta();
    c[0]->troba_celules_veines();
    c[1]->troba_celules_veines();
    cp[0]->troba_celules_veines();
    cp[1]->troba_celules_veines();
    
    if(1==0)
    {
      
    }
    else
    {
      if((c[0]==&p->celula_buida)&&(c[1]==&p->celula_buida))
      {
        std::cout << endl << "*******************************" << endl;
        std::cout << endl << "Erro: oooooops...there is an edge separating to empty cells!!" << endl;
        std::cout << endl << "Edge separating cells "<< c[0]->id << " " << c[1]->id << endl;
        std::cout << endl << "Vertexes "<< v[0]->id << " " << v[1]->id << endl;
        std::cout << endl << "Coordinates vertex 1: ("<< v[0]->x << ", " << v[0]->y << ")" << endl;
        std::cout << endl << "Coordinates vertex 2: ("<< v[1]->x << ", " << v[1]->y << ")" << endl;
        std::cout << endl << "Edge length: "<< l << "." << endl;
        std::cout << endl << "*******************************" << endl;
        exit(1);
      }
      else
      {
        comp1=((c[0]->ncellvertexes) > 3)||(c[0]==&p->celula_buida);
        comp2=((c[1]->ncellvertexes) > 3)||(c[1]==&p->celula_buida);
        if(comp1&&comp2)
        {
          //std::cout << endl << "Proces t1 a l'aresta amb id: " << id << " l: " << l << " canvi: " << canvi << endl;
          //std::cout << "que uneix els vertex: " << v[0]->id << " i " << v[1]->id << endl;
          //cout.flush();
          p->flog << p->stage << " " << p->idx_exterior << " " << p->idx_interior << ": " << "t1 process on edge " << id << " (" << v[0]->x << ", " << v[0]->y << ")-(" << v[1]->x << ", " << v[1]->y << ") that divides cells " << c[0]->c_track_id << " and " << c[1]->c_track_id << ". ";
          
          proces_t1();
          
          p->flog << "Coordinates after transition: (" << v[0]->x << ", " << v[0]->y << ")-(" << v[1]->x << ", " << v[1]->y << ") and divides cells " << c[0]->c_track_id << " and " << c[1]->c_track_id << "." << std::endl;
        }
        else
        {
          if(!comp1)
          {
            index=0;
          }
          else
          {
            index=1;
          }
          flag=0;
          for(i=0; i<3; i++)
          {
            if(c[index]->c[i]->narestes==3)
            {
              flag++;
              atemp=c[index]->a[i];
            }
          }
          if(flag==0)
          {
            //std::cout << endl << "Proces t2 a l'aresta amb id: " << id << " l: " << l << " canvi: " << canvi << endl;
            //fflush(stdout);
            p->flog << p->stage << " " << p->idx_exterior << " " << p->idx_interior << ": " << "t2 process on cell " << c[index]->c_track_id << " triggering edge " << id << " (" << v[0]->x << ", " << v[0]->y << ")-(" << v[1]->x << ", " << v[1]->y << ") that divides cells " << c[0]->c_track_id << " and " << c[1]->c_track_id  << "." << std::endl;
            proces_t2();
          }
          else if(flag==1)
          {
            /*std::cout << endl << "*******************************" << endl;
            std::cout << endl << "S'ha produit un proces t3." << endl;
            std::cout << endl << "S'ha parat la simulacio degut a un error en la implementació." << endl;
            std::cout << endl << "*******************************" << endl;
            exit(1);*/
            //std::cout << endl << "Proces t3 a l'aresta amb id: " << atemp->id << " l: " << atemp-> l << " canvi: " << atemp->canvi << endl;
            //fflush(stdout);
            p->flog << p->stage << " " << p->idx_exterior << " " << p->idx_interior << ": " << "t3 process on cell " << atemp->c[0]->c_track_id << " and cell " << atemp->c[0]->c_track_id << " triggering edge " << id << " (" << v[0]->x << ", " << v[0]->y << ")-(" << v[1]->x << ", " << v[1]->y << ") that divides cells " << c[0]->c_track_id << " and " << c[1]->c_track_id  << "." << std::endl;
            proces_t3(atemp);
          }
          else
          {
            std::cout << endl << "*******************************" << endl;
            std::cout << endl << "Error: an error has ocurred when identifying the topological process to trigger." << endl;
            std::cout << endl << "*******************************" << endl;
            exit(1);
          }
        }
      }
    }
    for(i=0; i<p->n_matriu_v; i++)
    {
      p->matriu_v[i].calcula_f0();
    }
  }
  
}

void aresta::actualitza_coordenades()
{
  int i;
  
  for(i=0;i<2;i++)
  {
    x[i]=v[i]->x;
    y[i]=v[i]->y;
  }
  
}

void aresta::crea_aresta(poblacio *pt, vertex *v1, vertex *v2, celula *c1, celula *c2)
{
  
  //ATENCIO!! Per cridar aquesta subrutina, es necessari que les celules c1 i c2
  //ja tinguin definits els vertexs v1 i v2. N'estas segur?
  //I les coordenades dels vertexs han d'estar definides tambe.
  
  int i,j;
  int flag1, flag2;
  
  p=pt;
  
  v[0]=v1;
  v[1]=v2;
  
  c[0]=c1;
  c[1]=c2;
  
  calcula_lambda();
  
  canvi=0;
  
  actualitza_coordenades();
  calcula_longitud();
  
  //Actualitzem en els objectes vertex que defineixen l'aresta. Introduim les celules a les cuals
  //pertany el vertex.
  
  for(j=0;j<2;j++){
    flag1=0;
    flag2=0;
    for(i=0;i<v[j]->ncelules;i++)
    {
      if(c1==v[j]->c[i]){
        flag1=1;
      }
      if(c2==v[j]->c[i]){
        flag2=1;
      }
    }
    if(flag1==0)
    {
      if(c1!=&p->celula_buida)
      {
        v[j]->c[v[j]->ncelules]=c1;
        v[j]->ncelules=v[j]->ncelules+1;
      }
    }
    if(flag2==0)
    {
      if(c2!=&p->celula_buida)
      {
        v[j]->c[v[j]->ncelules]=c2;
        v[j]->ncelules=v[j]->ncelules+1;
      }
    }
  }
  
  //Actualitzem en els objectes vertex que defineixen l'aresta. Afegim l'aresta actual, en el cas
  //que no hi fos. Tambe afegim el nou vertex al qual uneix l'aresta l'objecte vertex.
  //Si l'aresta ja estigues definida dintre de l'objecte vertex, actualitzem a el nou vertex.
  
  for(j=0;j<2;j++){
    flag1=0;
    for(i=0;i<v[j]->narestes;i++)
    {
      if(this==v[j]->a[i]){
        flag1=1;
        flag2=i;
      }
    }
    if(flag1==0)
    {
      v[j]->a[v[j]->narestes]=this;
      v[j]->v[v[j]->narestes]=v[(j+1)%2];
      v[j]->narestes=v[j]->narestes+1;
    }
    else if((flag1!=0)&&(v[j]->narestes>0))
    {
      v[j]->v[i]=v[(j+1)%2];
    }
  }
  
  //Aixo esta malament?
  //
  
  for(j=0;j<2;j++){
    if(c[j]!=&p->celula_buida)
    {
      flag1=0;
      for(i=0; i<c[j]->narestes; i++)
      {
        if(this==c[j]->a[i]){
          flag1=1;
        }
      }
      if(flag1==0)
      {
        c[j]->a[c[j]->narestes]=this;
        c[j]->narestes=c[j]->narestes+1;
      }
    }
  }
  
}

void aresta::escriu_aresta(ofstream &arxiu)
{
  int i;
  
  arxiu << "2" << " ";
  
  for(i=0; i<2; i++){
    arxiu << c[i]->ctype << " ";
  }
  
  arxiu << "2" << " ";
  
  for(i=0; i<2; i++){
    arxiu << v[i]->id << " " << v[i]->x << " " << v[i]->y << " ";
  }
  
  arxiu << std::endl;
  
}

void aresta::escriu_informacio_aresta(ofstream &arxiu)
{
  int i;
  
  if((id<p->n_matriu_a)&&(id>0))
  {
    arxiu << "Aresta: " << id << " l: " << l << std::endl;
    arxiu << "-----------" << std::endl;
    arxiu << "      2 celules veines: " << std::endl;
    for(i=0; i<2; i++)
    {
      arxiu << "              Celula: " << c[i]->id << std::endl;
    }
    arxiu << "      2 celules a les puntes: " << std::endl;
    for(i=0; i<2; i++)
    {
      arxiu << "              Celula: " << cp[i]->id << std::endl;
    }
    arxiu << "      2 vertexs adjacents: " << std::endl;
    arxiu << "              Vertex: " << v[0]->id << " que va a parar al vertex: " << v[1]->id << std::endl;
  }
  else
  {
    arxiu << "Aresta: " << id << " l: " << l << std::endl;
  }
}

void aresta::troba_celules_puntes_aresta()
{
  int i,j;
  int nctemp;
  celula *ctemp[MAXIM_CELULES_V];
  
  
  for(j=0; j<2; j++)
  {
    nctemp=0;
    for(i=0; i<v[j]->ncelules; i++)
    {
      if((v[j]->c[i]!=c[0])&&(v[j]->c[i]!=c[1]))
      {
        ctemp[nctemp]=v[j]->c[i];
        nctemp++;
      }
    }
    if(nctemp==0)
    {
      cp[j]=&p->celula_buida;
    }
    else if(nctemp==1)
    {
      cp[j]=ctemp[0];
    }
    else
    {
      std::cout << endl << "*******************************" << endl;
      std::cout << endl << "Error: something went wrong when computing the leading cells of the edges!!! (nctemp): " << nctemp << endl;
      std::cout << endl << "Edge: " << id << endl;
      std::cout << endl << "Neighbors: " << c[0]->id << " i " << c[1]->id << endl;
      std::cout << endl << "Vertexes: " << v[0]->id << " i " << v[1]->id << endl;
      std::cout << endl << "Specifically when looking the leading cell defined by vertex: " << v[j]->id << endl;
      std::cout << "having " << v[j]->ncelules << " neighbors cells:" << endl;
      for(i=0; i<v[j]->ncelules; i++)
      {
        std::cout << "\tcelula: " << v[j]->c[i]->id << endl;
      }
      std::cout << endl << "There are: " << p->n_matriu_a << " edges," << endl;
      std::cout << p->n_matriu_c << " cells and" << endl;
      std::cout << p->n_matriu_v << " vertexes." << endl;
      std::cout << endl << "*******************************" << endl;
      exit(1);
    }
    
  }
  
}

void aresta::canvia_conexio(vertex *vo, vertex *vv, vertex *vn)
{
  if(v[0]==vo)
  {
    if(v[1]==vv)
    {
      v[1]=vn;
    }
    else
    {
      std::cout << endl << "*******************************" << endl;
      std::cout << endl << "Error: when changing an edge connection!! (1: at the end of edge)" << endl;
      std::cout << endl << "*******************************" << endl;
      exit(1);
    }
  }
  else if(v[1]==vo)
  {
    if(v[0]==vv)
    {
      v[0]=vn;
    }
    else
    {
      std::cout << endl << "*******************************" << endl;
      std::cout << endl << "Error: when changing an edge connection!! (2: at the end of edge)" << endl;
      std::cout << endl << "*******************************" << endl;
      exit(1);
    }
  }
  else
  {
    std::cout << endl << "*******************************" << endl;
      std::cout << endl << "Error: when changing an edge connection!! (3: at the origin of edge)" << endl;
    std::cout << endl << "*******************************" << endl;
    exit(1);
  }
}

void aresta::canvia_celules_veines(celula *cv, celula *cn)
{
  if(c[1]==cv)
  {
    c[1]=cn;
  }
  else if(c[0]==cv)
  {
    c[0]=cn;
  }
  else
  {
    std::cout << endl << "*******************************" << endl;
    std::cout << endl << "Error: when swaping the cell neighbors of an edge!!" << endl;
    std::cout << endl << "The old cell does not belong to the edge." << endl;
    std::cout << endl << "*******************************" << endl;
    exit(1);
  }
}

punt aresta::comprova_si_intersecta_amb_recta(double alfa, double beta)
{
  punt p;
  double a, b;    //alfa i beta de la recta que passa per l'aresta.
  double tx[2];   //Variable que serveix per guardar els dos punts x del vertex, ordenats de menor a major.
  
  p.actiu=0;
  
  a=((v[1]->y)-(v[0]->y))/((v[1]->x)-(v[0]->x));
  b=(v[0]->y)-((v[0]->x)*a);
  
  if((v[1]->x)>(v[0]->x))
  {
    tx[0]=v[0]->x;
    tx[1]=v[1]->x;
  }
  else
  {
    tx[0]=v[1]->x;
    tx[1]=v[0]->x;
  }
  
  p.x=(beta-b)/(a-alfa);
  
  if((p.x<tx[1])&&(p.x>=tx[0]))
  {
    p.y=a*p.x+b;
    p.actiu=1;
  }
  
  
  
  return p;
}

void aresta::cambia_referencia_aresta(aresta *dv, aresta *dn)
{
  int i, j;
  
  //troba_celules_puntes_aresta();
  
  for(i=0; i<2; i++)
  {
    for(j=0; j<c[i]->narestes; j++)
    {
      if(c[i]->a[j]==dv)
      {
        c[i]->a[j]=dn;
      }
    }
  }
  
  for(i=0; i<2; i++)
  {
    for(j=0; j<cp[i]->narestes; j++)
    {
      if(cp[i]->a[j]==dv)
      {
        cp[i]->a[j]=dn;
      }
    }
  }
  
  for(i=0; i<2; i++)
  {
    for(j=0; j<v[i]->narestes; j++)
    {
      if(v[i]->a[j]==dv)
      {
        v[i]->a[j]=dn;
      }
    }
  }
}

void aresta::calcula_lambda()
{
  
  
  lambdav=p->k_lambda_v[c[0]->ctype][c[1]->ctype];
  lambdah=p->k_lambda_h[c[0]->ctype][c[1]->ctype];
  k_gamma_aresta=p->k_gamma_aresta[c[0]->ctype][c[1]->ctype];
  
  templh2lv=(lambdah-2.*lambdav);
  templv2lh=(lambdav-2.*lambdah);
  
  calcula_longitud();
}

void aresta::proces_t1()
{
  int i, j;
  vertex *vtemp[2];       //Vertexs a canviar de la celula "i" a part dels implicats en l'aresta
  aresta *atemp[2];
  double vx, vy;          //Vector que apunta al centre de l'aresta.
  double vax, vay;        //Vector perpendicular a l'aresta actual i de modul igual.
  double temp1, temp2;
  int index;
  double desfesx[2], desfesy[2];
  
  ofstream error;
  
  vtemp[0]=NULL;
  vtemp[1]=NULL;
  atemp[0]=NULL;
  atemp[1]=NULL;
  
  desfesx[0] = v[0]->x;
  desfesy[0] = v[0]->y;
  desfesx[1] = v[1]->x;
  desfesy[1] = v[1]->y;
  
  //Mirem quina de les dos possibles celules que separa l'aresta no esta buida i en calculem el centre (aixo s'ha de fer abans de moure els vertexs. Ara esta be.)
  
  if(c[0]->id!=-1)
  {
    index=0;
  }
  else
  {
    index=1;
  }
  
  c[index]->calcula_centre();
  
  //Movem els vertexs de tal forma que l'aresta actual passi a ser perpendicular a si mateixa passant pel seu centre.
  
  vx=v[0]->x+(v[1]->x-v[0]->x)/2.;
  vy=v[0]->y+(v[1]->y-v[0]->y)/2.;
  
  vax=-(v[1]->y-v[0]->y);
  vay=(v[1]->x-v[0]->x);
  
  v[0]->x=vx+vax/2.;
  v[0]->y=vy+vay/2.;
  
  v[1]->x=vx-vax/2.;
  v[1]->y=vy-vay/2.;
  
  //Comprovem que el vertex v[0] es el que ha quedat mes a prop de la celula c[0] i viceversa.
  
  temp1=((c[index]->x-v[index]->x)*(c[index]->x-v[index]->x)+(c[index]->y-v[index]->y)*(c[index]->y-v[index]->y));
  temp2=((c[index]->x-v[(index+1)%2]->x)*(c[index]->x-v[(index+1)%2]->x)+(c[index]->y-v[(index+1)%2]->y)*(c[index]->y-v[(index+1)%2]->y));
  
  if(temp1>temp2)
  {
    temp1=v[0]->x;
    temp2=v[0]->y;
    v[0]->x=v[1]->x;
    v[0]->y=v[1]->y;
    v[1]->x=temp1;
    v[1]->y=temp2;
  }
  
  //Fem diferents coses en funcio de si hi han celules buides i a on.
  
  if((cp[0]!=&p->celula_buida)&&(cp[1]!=&p->celula_buida))
  {
    //std::cout << "Cas 1" << std::endl;
    
    //Busquem els altres vertexs de les celules c[0] i c[1] que s'han de reconectar
    
    if((c[0]!=&p->celula_buida)||(cp[1]!=&p->celula_buida))
    {
      vtemp[0]=v[1]->troba_vertex_celula_celula(c[0], cp[1]);
    }
    else
    {
      vtemp[0]=&p->vertex_buit;
    }
    if((c[1]!=&p->celula_buida)||(cp[0]!=&p->celula_buida))
    {
      vtemp[1]=v[0]->troba_vertex_celula_celula(c[1], cp[0]);
    }
    else
    {
      vtemp[1]=&p->vertex_buit;
    }
    
    
    
    for(i=0; i<vtemp[0]->narestes; i++)
    {
      if(vtemp[0]->v[i]==v[1])
      {
        vtemp[0]->v[i]=v[0];
        vtemp[0]->a[i]->canvia_conexio(vtemp[0], v[1], v[0]);
        atemp[0]=vtemp[0]->a[i];
      }
    }
    
    for(i=0; i<v[0]->narestes; i++)
    {
      if(v[0]->v[i]==vtemp[1])
      {
        v[0]->v[i]=vtemp[0];
        v[0]->a[i]=atemp[0];
      }
    }
    
    for(i=0; i<vtemp[1]->narestes; i++)
    {
      if(vtemp[1]->v[i]==v[0])
      {
        vtemp[1]->v[i]=v[1];
        vtemp[1]->a[i]->canvia_conexio(vtemp[1], v[0], v[1]);
        atemp[1]=vtemp[1]->a[i];
      }
    }
    
    for(i=0; i<v[1]->narestes; i++)
    {
      if(v[1]->v[i]==vtemp[0])
      {
        v[1]->v[i]=vtemp[1];
        v[1]->a[i]=atemp[1];
      }
    }
    
    //Canviem les celules veines dels vertexs v[0] i v[1]
    
    //Constants. Recalcular parametres forces...
    v[0]->canvia_celules_veines(c[1], cp[1]);
    v[1]->canvia_celules_veines(c[0], cp[0]);
    
    cp[0]->introdueix_vertex(v[1], v[0], vtemp[1]);
    cp[1]->introdueix_vertex(v[0], v[1], vtemp[0]);
    c[0]->elimina_vertex(v[1]);
    c[1]->elimina_vertex(v[0]);
    
    
    
    c[0]=cp[0];
    c[1]=cp[1];
    
    troba_celules_puntes_aresta();
    cp[0]->calcula_propietats();
    cp[1]->calcula_propietats();
    c[0]->calcula_propietats();
    c[1]->calcula_propietats();
    
    v[0]->troba_vertexs_veins();
    v[1]->troba_vertexs_veins();
    vtemp[0]->troba_vertexs_veins();
    vtemp[1]->troba_vertexs_veins();
    
    for(i=0; i<v[0]->narestes; i++)
    {
      v[0]->a[i]->troba_celules_puntes_aresta();
      v[0]->a[i]->calcula_longitud();
    }
    for(i=0; i<v[1]->narestes; i++)
    {
      v[1]->a[i]->troba_celules_puntes_aresta();
      v[1]->a[i]->calcula_longitud();
    }
    
    calcula_lambda();
  }
  else if((cp[0]==&p->celula_buida)&&(cp[1]==&p->celula_buida))
  {
    //std::cout << "Cas 2" << std::endl;
    v[0]->x = desfesx[0];
    v[0]->y = desfesy[0];
    v[1]->x = desfesx[1];
    v[1]->y = desfesy[1];
    if((c[0]->id==-1)||(c[1]->id==-1))
    {
      /****************************************************************************************
       ****************************************************************************************/
      #ifdef ALLOWT1C2
      if(c[0]->id==-1)
      {
        index=1;
      }
      else
      {
        index=0;
      }
      
      if(v[1]->narestes==2)
      {
        vtemp[0]=((v[1]->v[0]==v[0])?v[1]->v[1]:v[1]->v[0]);
      }
      else
      {
        std::cout << endl << "*******************************" << endl;
        std::cout << endl << "Error: t1 process at empty cells:" << endl;
        std::cout << endl << "may be due to the increase in the number of possible neighbors vertexes." << endl;
        std::cout << endl << "*******************************" << endl;
        exit(1);
      }
      
      for(i=0; i<vtemp[0]->narestes; i++)
      {
        if(vtemp[0]->v[i]==v[1])
        {
          vtemp[0]->v[i]=v[0];
          vtemp[0]->a[i]->canvia_conexio(vtemp[0], v[1], v[0]);
          atemp[0]=vtemp[0]->a[i];
        }
      }
      
      for(i=0; i<v[0]->narestes; i++)
      {
        if(v[0]->v[i]==v[1]) 
        {
          v[0]->v[i]=vtemp[0];
          v[0]->a[i]=atemp[0];
        }
      }
      
      atemp[0]->calcula_longitud();
      
      c[index]->elimina_vertex(v[1]);
      c[index]->calcula_propietats();
      
      vtemp[0]->troba_vertexs_veins();
      v[0]->troba_vertexs_veins();
      
      p->vertexs_per_destruir[0]=v[1]->id;
      p->n_vertexs_destruir=1;
      
      p->arestes_per_destruir[0]=id;
      p->n_arestes_destruir=1;
      
      p->destrueix_elements();
      
      
      
      
      #else
      std::cout << endl << "*******************************" << endl;
      std::cout << endl << "Warning: t1 (case 2) was supposed to occur but it has been deactivated and will not take place." << endl;
      std::cout << endl << "*******************************" << endl;
      #endif
      /****************************************************************************************
       ****************************************************************************************/
    }
    else
    {
      
      /****************************************************************************************
       ****************************************************************************************/
      
      p->matriu_v[p->n_matriu_v].inicia_vertex(p,p->n_matriu_v,v[0]->x,v[0]->y);
      p->matriu_v[p->n_matriu_v+1].inicia_vertex(p,p->n_matriu_v+1,v[1]->x,v[1]->y);
      for(i=0; i<2; i++)
      {
        p->matriu_v[p->n_matriu_v + i].copia_des_de_vertex(v[i]);
      }
      
      p->matriu_a[p->n_matriu_a].p=p;
      p->matriu_a[p->n_matriu_a].id=p->n_matriu_a;
      p->matriu_a[p->n_matriu_a].copia_des_de_aresta(this);
      
      vtemp[0] = v[0]->troba_vertex_celula_celula(c[1], &(p->celula_buida));
      vtemp[1] = v[1]->troba_vertex_celula_celula(c[1], &(p->celula_buida));
      atemp[0] = v[0]->troba_aresta_vertex_oposat(c[1], v[1]);
      atemp[1] = v[1]->troba_aresta_vertex_oposat(c[1], v[0]);
      
      //A les celules s'han de substituir les referencies als vertexs i arestes. (FET)
      for(i=0; i<c[1]->ncellvertexes; i++)
      {
        for(j=0; j<2; j++)
        {
          if(c[1]->v[i] == v[j])
          {
            c[1]->v[i] = &(p->matriu_v[p->n_matriu_v + j]);
          }
        }
      }
      for(i=0; i<c[1]->narestes; i++)
      {
        if(c[1]->a[i] == this)
        {
          c[1]->a[i] = &(p->matriu_a[p->n_matriu_a]);
        }
      }
      
      
      
      //Als vertes s'han d'eliminar les arestes i celules que ja no hi son.
      
      v[0]->elimina_aresta(v[0]->troba_aresta_vertex_oposat(c[1], v[1]));
      v[0]->canvia_celules_veines(c[1], &(p->celula_buida));
      v[0]->troba_vertexs_veins();
      
      v[1]->elimina_aresta(v[1]->troba_aresta_vertex_oposat(c[1], v[0]));
      v[1]->canvia_celules_veines(c[1], &(p->celula_buida));
      v[1]->troba_vertexs_veins();
      
      p->matriu_v[p->n_matriu_v].elimina_aresta(p->matriu_v[p->n_matriu_v].troba_aresta_vertex_oposat(c[0], v[1]));
      p->matriu_v[p->n_matriu_v].canvia_conexio_amb_vertex(v[1], this, &(p->matriu_v[p->n_matriu_v + 1]), &(p->matriu_a[p->n_matriu_a]));
      p->matriu_v[p->n_matriu_v].canvia_celules_veines(c[0], &(p->celula_buida));
      p->matriu_v[p->n_matriu_v].troba_vertexs_veins();
      
      p->matriu_v[p->n_matriu_v + 1].elimina_aresta(p->matriu_v[p->n_matriu_v+1].troba_aresta_vertex_oposat(c[0], v[0]));
      p->matriu_v[p->n_matriu_v + 1].canvia_conexio_amb_vertex(v[0], this, &(p->matriu_v[p->n_matriu_v]), &(p->matriu_a[p->n_matriu_a]));
      p->matriu_v[p->n_matriu_v + 1].canvia_celules_veines(c[0], &(p->celula_buida));
      p->matriu_v[p->n_matriu_v + 1].troba_vertexs_veins();
      
      //Als vertexs nous s'ha de cambiar la referencia a l'aresta nova.
      
      //A les arestes s'han d'actualitzar les celules i els vertexs.
      //Un cop fet aixo, les constants de l'aresta nova i vella cambien, i per tant s'haurien de calcular.!!!!
      //No cal canviar celules punta aresta ja que son celules buides i ja s'han copiat.
      canvia_celules_veines(c[1], &(p->celula_buida));
      calcula_lambda();
      p->matriu_a[p->n_matriu_a].canvia_celules_veines(c[0], &(p->celula_buida));
      p->matriu_a[p->n_matriu_a].v[0] = &(p->matriu_v[p->n_matriu_v]);
      p->matriu_a[p->n_matriu_a].v[1] = &(p->matriu_v[p->n_matriu_v + 1]);
      p->matriu_a[p->n_matriu_a].calcula_lambda();
      
      atemp[0]->canvia_conexio(vtemp[0], v[0], &(p->matriu_v[p->n_matriu_v]));
      atemp[1]->canvia_conexio(vtemp[1], v[1], &(p->matriu_v[p->n_matriu_v+1]));
      vtemp[0]->canvia_conexio_amb_vertex(v[0], atemp[0], &(p->matriu_v[p->n_matriu_v]), atemp[0]);
      vtemp[1]->canvia_conexio_amb_vertex(v[1], atemp[1], &(p->matriu_v[p->n_matriu_v+1]), atemp[1]);
      
      c[0]->troba_celules_veines();
      c[1]->troba_celules_veines();
      p->matriu_a[p->n_matriu_a].c[0]->troba_celules_veines();
      p->matriu_a[p->n_matriu_a].c[1]->troba_celules_veines();
      
      p->n_matriu_v = p->n_matriu_v + 2;
      p->n_matriu_a = p->n_matriu_a + 1;
      /****************************************************************************************
       ****************************************************************************************/
    }
  }
  else
  {
    //std::cout << "Cas 3:" << " c[0] " << c[0]->id << " " << " c[1] " << c[1]->id << " "  << " cp[0] " << cp[0]->id << " "  << " cp[1] " << cp[1]->id << " "  << std::endl;
    //Busquem els altres vertexs de les celules c[0] i c[1] que s'han de reconectar
    
    if((c[0]!=&p->celula_buida)||(cp[1]!=&p->celula_buida))
    {
      if(v[1]->narestes>2)
      {
        vtemp[0]=v[1]->troba_vertex_celula_celula(c[0], cp[1]);
      }
      else if(v[1]->narestes==2)
      {
        vtemp[0]=((v[1]->v[0]==v[0])?v[1]->v[1]:v[1]->v[0]);
      }
      else
      {
        std::cout << "Error: there was a problem during triggering a t1 (case 3) process."  << std::endl;
        exit(1);
      }
    }
    else
    {
      vtemp[0]=&p->vertex_buit;
    }
    
    if((c[1]!=&p->celula_buida)||(cp[0]!=&p->celula_buida))
    {
      if(v[0]->narestes>2)
      {
        vtemp[1]=v[0]->troba_vertex_celula_celula(c[1], cp[0]);
      }
      else if(v[0]->narestes==2)
      {
        vtemp[1]=((v[0]->v[0]==v[1])?v[0]->v[1]:v[0]->v[0]);
      }
      else
      {
        std::cout << "Error: there was a problem during triggering a t1 (case 3) process."  << std::endl;
        exit(1);
      }
    }
    else
    {
      vtemp[1]=&p->vertex_buit;
    }
    
    
    
    if(vtemp[0]->id!=-1)
    {
      for(i=0; i<vtemp[0]->narestes; i++)
      {
        if(vtemp[0]->v[i]==v[1])
        {
          vtemp[0]->v[i]=v[0];
          vtemp[0]->a[i]->canvia_conexio(vtemp[0], v[1], v[0]);
          atemp[0]=vtemp[0]->a[i];
        }
      }
      //aixo esta malament!!! Pot ser que siguin els dos vertexs vtemp reals.
      if(vtemp[1]->id==-1)
      {
        v[0]->v[v[0]->narestes]=vtemp[0];
        v[0]->a[v[0]->narestes]=atemp[0];
        v[0]->narestes=v[0]->narestes+1;
      }
      else
      {
        for(i=0; i<v[0]->narestes; i++)
        {
          if(v[0]->v[i]==vtemp[1])
          {
            v[0]->v[i]=vtemp[0];
            v[0]->a[i]=atemp[0];
          }
        }
      }
    }
    else
    {
      v[0]->canvia_conexio_amb_vertex(vtemp[1], &p->aresta_buida, vtemp[0], &p->aresta_buida);
    }
    
    if(vtemp[1]->id!=-1)
    {
      for(i=0; i<vtemp[1]->narestes; i++)
      {
        if(vtemp[1]->v[i]==v[0])
        {
          vtemp[1]->v[i]=v[1];
          vtemp[1]->a[i]->canvia_conexio(vtemp[1], v[0], v[1]);
          atemp[1]=vtemp[1]->a[i];
        }
      }
      
      if(vtemp[0]->id==-1)
      {
        v[1]->v[v[1]->narestes]=vtemp[1];
        v[1]->a[v[1]->narestes]=atemp[1];
        v[1]->narestes=v[1]->narestes+1;
      }
      else
      {
        for(i=0; i<v[1]->narestes; i++)
        {
          if(v[1]->v[i]==vtemp[0])
          {
            v[1]->v[i]=vtemp[1];
            v[1]->a[i]=atemp[1];
          }
        }
      }
    }
    else
    {
      v[1]->canvia_conexio_amb_vertex(vtemp[0], &p->aresta_buida, vtemp[1], &p->aresta_buida);
    }
    
    
    //Canviem les celules veines dels vertexs v[0] i v[1]
    
    v[0]->canvia_celules_veines(c[1], cp[1]);
    v[1]->canvia_celules_veines(c[0], cp[0]);
    
    //Constants. Recalcular parametres forces...
    
    cp[0]->introdueix_vertex(v[1], v[0], vtemp[1]);
    cp[1]->introdueix_vertex(v[0], v[1], vtemp[0]);
    c[0]->elimina_vertex(v[1]);
    c[1]->elimina_vertex(v[0]);
    
    
    
    c[0]=cp[0];
    c[1]=cp[1];
    
    troba_celules_puntes_aresta();
    cp[0]->calcula_propietats();
    cp[1]->calcula_propietats();
    c[0]->calcula_propietats();
    c[1]->calcula_propietats();
    
    v[0]->troba_vertexs_veins();
    v[1]->troba_vertexs_veins();
    vtemp[0]->troba_vertexs_veins();
    vtemp[1]->troba_vertexs_veins();
    
    for(i=0; i<v[0]->narestes; i++)
    {
      v[0]->a[i]->troba_celules_puntes_aresta();
      v[0]->a[i]->calcula_longitud();
    }
    for(i=0; i<v[1]->narestes; i++)
    {
      v[1]->a[i]->troba_celules_puntes_aresta();
      v[1]->a[i]->calcula_longitud();
    }
    
    calcula_lambda();
  }
  
}

void aresta::proces_t2()
{
  int i;
  vertex *vtemp[2];       //Vertexs a canviar de la celula "i" a part dels implicats en l'aresta
  aresta *atemp[2];
  vertex *vt;             //Vertex que conservarem despres del proces t2. Aquest es l'oposat a l'aresta actual de la celula triangular.
  int aresta_id_eliminar[2];
  bool cond1, cond2;
  
  vtemp[0]=NULL;
  vtemp[1]=NULL;
  atemp[0]=NULL;
  atemp[1]=NULL;
  vt=NULL;
  aresta_id_eliminar[0]=-1;
  aresta_id_eliminar[1]=-1;
  
  cond1=(((c[0]->ncellvertexes) < 4)&&(c[0]->id!=-1));
  cond2=(((c[1]->ncellvertexes) < 4)&&(c[1]->id!=-1));
  
  if(cond1&&cond2)
  {
    std::cout << endl << "*******************************" << endl;
    std::cout << endl << "Error: well...this is...very embarrassing..." << endl;
    std::cout << endl << "a 4 corners vertex??????..." << endl;
    std::cout << endl << "ummm....not really sure..." << endl;
    std::cout << endl << "*******************************" << endl;
    exit(1);
    
  }
  else if(cond1||cond2)
  {
    int index;
    
    troba_celules_puntes_aresta();
    
    //Mirem quina de les dos celules es la que te tres vertexs.
    if(cond1)
    {
      index=0;
    }
    else
    {
      index=1;
    }
    
    //Busca el vertex que conservarem i l'asignem a vt
    for(i=0; i<c[index]->ncellvertexes; i++)
    {
      if((c[index]->v[i]!=v[0])&&(c[index]->v[i]!=v[1]))
      {
        vt=c[index]->v[i];
      }
    }
    //Posiciona vt al centre de la celula que desapareix
    c[index]->calcula_centre();
    vt->x=c[index]->x;
    vt->y=c[index]->y;
    //Busca els vertexs que reconectarem amb vt
    vtemp[0]=v[0]->troba_vertex_celula_celula(c[(index+1)%2], cp[0]);
    vtemp[1]=v[1]->troba_vertex_celula_celula(c[(index+1)%2], cp[1]);
    
    //Reconecta vt amb vtemp[] i cambia les propietats de les arestes que reconectem.
    for(i=0; i<vtemp[0]->narestes; i++)
    {
      if(vtemp[0]->v[i]==v[0])
      {
        vtemp[0]->v[i]=vt;
        vtemp[0]->a[i]->canvia_conexio(vtemp[0], v[0], vt);
        atemp[0]=vtemp[0]->a[i];
      }
    }
    for(i=0; i<vt->narestes; i++)
    {
      if(vt->v[i]==v[0])
      {
        vt->v[i]=vtemp[0];
        aresta_id_eliminar[0]=vt->a[i]->id;
        vt->a[i]=atemp[0];
      }
    }
    
    for(i=0; i<vtemp[1]->narestes; i++)
    {
      if(vtemp[1]->v[i]==v[1])
      {
        vtemp[1]->v[i]=vt;
        vtemp[1]->a[i]->canvia_conexio(vtemp[1], v[1], vt);
        atemp[1]=vtemp[1]->a[i];
      }
    }
    for(i=0; i<vt->narestes; i++)
    {
      if(vt->v[i]==v[1])
      {
        vt->v[i]=vtemp[1];
        aresta_id_eliminar[1]=vt->a[i]->id;
        vt->a[i]=atemp[1];
      }
    }
    
    //Elimina els vertexs sobrants v[] en les celules que queden vives
    cp[0]->elimina_vertex(v[0]);
    cp[1]->elimina_vertex(v[1]);
    for(i=0; i<c[(index+1)%2]->ncellvertexes; i++)
    {
      if(c[(index+1)%2]->v[i]==v[1])
      {
        c[(index+1)%2]->v[i]=vt;
      }
    }
    c[(index+1)%2]->elimina_vertex(v[0]);
    
    /*Per calcular propietats segurament cal definir be les celules punta de l'aresta!!! S'ha de comprovar*/
    
    
    
    //Cambiem les celules veines del vertex vt.
    
    vt->canvia_celules_veines(c[index], c[(index+1)%2]);
    
    /*Per calcular propietats segurament cal definir be les celules punta de l'aresta!!! S'ha de comprovar*/
    
    troba_celules_puntes_aresta();
    cp[0]->calcula_propietats();
    cp[1]->calcula_propietats();
    c[(index+1)%2]->calcula_propietats();
    
    vt->troba_vertexs_veins();
    vtemp[0]->troba_vertexs_veins();
    vtemp[1]->troba_vertexs_veins();
    
    //Trobem les celules punta de les arestes.
    
    for(i=0; i<vt->narestes; i++)
    {
      vt->a[i]->troba_celules_puntes_aresta();
      vt->a[i]->calcula_longitud();
    }
    
    //Eliminem els vertexs, arestes i celules que han desaparegut amb la transicio t2
    //de les matrius.
    
    p->celules_per_destruir[0]=c[index]->id;
    p->n_celules_destruir=1;
    
    p->vertexs_per_destruir[0]=v[0]->id;
    p->vertexs_per_destruir[1]=v[1]->id;
    p->n_vertexs_destruir=2;
    
    p->arestes_per_destruir[0]=aresta_id_eliminar[0];
    p->arestes_per_destruir[1]=aresta_id_eliminar[1];
    p->arestes_per_destruir[2]=id;
    p->n_arestes_destruir=3;
    
    p->destrueix_elements();
    
    
    
  }
  else
  {
    std::cout << endl << "*******************************" << endl;
    std::cout << endl << "Error: paranormal activity....this is really weird!!!" << endl;
    std::cout << endl << "Edge dividing cells "<< c[0]->id << " " << c[1]->id << endl;
    std::cout << endl << "with the number of vertexes: "<< c[0]->ncellvertexes << " " << c[1]->ncellvertexes << endl;
    std::cout << endl << "*******************************" << endl;
    exit(1);
  }
  
}

void aresta::proces_t3(aresta *acentral)
{
  int i,j;
  int flag;
  bool cond1, cond2, cond3, cond4;
  
  vertex *vmantenir[2];
  aresta *amantenir[2];
  
  vertex *veliminar1[2];
  vertex *veliminar2[2];
  aresta *aeliminar[2][2];
  celula *celiminar[2];
  
  //Iniciem algunes variables:
  
  //Les celules que s'eliminaran:
  
  celiminar[0]=acentral->c[0];
  celiminar[1]=acentral->c[1];
  
  //Els vertexs que s'eliminaran i que estan conectats per l'aresta central, acentral:
  
  veliminar1[0]=acentral->v[0];
  veliminar1[1]=acentral->v[1];
  
  //Els vertexs que s'eliminaran:
  
  flag=0;
  for(j=0; j<2; j++)
  {
    for(i=0; i<3; i++)
    {
      if((celiminar[j]->v[i]!=veliminar1[0])&&(celiminar[j]->v[i]!=veliminar1[1]))
      {
        veliminar2[j]=celiminar[j]->v[i];
        flag++;
      }
    }
  }
  if(flag!=2)
  {
    std::cout << endl << "*******************************" << endl;
    std::cout << endl << "Error: when implementing a t3 process!!" << endl;
    std::cout << endl << "*******************************" << endl;
    exit(1);
  }
  
  cond1=(veliminar1[0]->narestes!=3);
  cond2=(veliminar1[1]->narestes!=3);
  cond3=(veliminar2[0]->narestes!=3);
  cond4=(veliminar2[1]->narestes!=3);
  
  if(cond1||cond2||cond3||cond4)
  {
    std::cout << endl << "*******************************" << endl;
    std::cout << endl << "Error: when implementing a t3 process!!" << endl;
    std::cout << endl << "Each of the vertexes to kill does not have three edges." << endl;
    std::cout << endl << "Ummmm, this is really weird!!!" << endl;
    std::cout << endl << "*******************************" << endl;
    exit(1);
  }
  
  //Busquem les arestes que s'eliminaran, i els vertexs i arestes que es "mantindran".
  
  flag=0;
  for(j=0; j<2; j++)
  {
    for(i=0; i<3; i++)
    {
      if(veliminar2[j]->v[i]==veliminar1[0])
      {
        aeliminar[j][0]=veliminar2[j]->a[i];
        flag++;
      }
      else if(veliminar2[j]->v[i]==veliminar1[1])
      {
        aeliminar[j][1]=veliminar2[j]->a[i];
        flag++;
      }
      else if((veliminar2[j]->v[i]!=veliminar1[0])&&(veliminar2[j]->v[i]!=veliminar1[1]))
      {
        vmantenir[j]=veliminar2[j]->v[i];
        amantenir[j]=veliminar2[j]->a[i];
        flag++;
      }
    }
  }
  if(flag!=6)
  {
    std::cout << endl << "*******************************" << endl;
    std::cout << endl << "Error: when implementing a t3 process!!" << endl;
    std::cout << endl << "edges to kill missing!!" << endl;
    std::cout << endl << "*******************************" << endl;
    exit(1);
  }
  
  //Reconectem vmantenir[0] i vmantenir[1]
  
  flag=0;
  for(i=0; i<vmantenir[0]->narestes; i++)
  {
    if(vmantenir[0]->v[i]==veliminar2[0])
    {
      vmantenir[0]->v[i]=vmantenir[1];
      flag++;
    }
  }
  for(i=0; i<vmantenir[1]->narestes; i++)
  {
    if(vmantenir[1]->v[i]==veliminar2[1])
    {
      vmantenir[1]->v[i]=vmantenir[0];
      vmantenir[1]->a[i]=amantenir[0];
      flag++;
    }
  }
  if(flag!=2)
  {
    std::cout << endl << "*******************************" << endl;
    std::cout << endl << "Error: when implementing a t3 process!!" << endl;
    std::cout << endl << "trying to reconnect the remaining vertexes!!" << endl;
    std::cout << endl << "*******************************" << endl;
    exit(1);
  }
  
  if(amantenir[0]->v[0]==vmantenir[0])
  {
    amantenir[0]->v[1]=vmantenir[1];
  }
  else
  {
    amantenir[0]->v[0]=vmantenir[1];
  }
  
  //Eliminem de les celules veines, els vertexs que s'eliminaran. I els hi trobem les noves arestes.
  
  acentral->troba_celules_puntes_aresta();
  
  acentral->cp[0]->elimina_vertex(veliminar1[0]);
  acentral->cp[1]->elimina_vertex(veliminar1[1]);
  
  acentral->cp[0]->elimina_vertex(veliminar2[0]);
  acentral->cp[0]->elimina_vertex(veliminar2[1]);
  acentral->cp[1]->elimina_vertex(veliminar2[0]);
  acentral->cp[1]->elimina_vertex(veliminar2[1]);
  
  //Recalculem longituds i propietats de les celules i arestes afectades.
  
  amantenir[0]->troba_celules_puntes_aresta();
  amantenir[0]->calcula_longitud();
  
  for(i=0; i<2; i++)
  {
    if(acentral->cp[i]->id!=-1)
    {
      acentral->cp[i]->calcula_propietats();
      
    }
  }
  
  //Es calculen els vertex veins dels vertexs que es mantenen;
  
  vmantenir[0]->troba_vertexs_veins();
  vmantenir[1]->troba_vertexs_veins();
  
  //Eliminem de les matrius les celules, arestes i vertexs que han desaparegut.
  
  p->celules_per_destruir[0]=celiminar[0]->id;
  p->celules_per_destruir[1]=celiminar[1]->id;
  p->n_celules_destruir=2;
  
  p->vertexs_per_destruir[0]=veliminar1[0]->id;
  p->vertexs_per_destruir[1]=veliminar1[1]->id;
  p->vertexs_per_destruir[2]=veliminar2[0]->id;
  p->vertexs_per_destruir[3]=veliminar2[1]->id;
  p->n_vertexs_destruir=4;
  
  p->arestes_per_destruir[0]=aeliminar[0][0]->id;
  p->arestes_per_destruir[1]=aeliminar[0][1]->id;
  p->arestes_per_destruir[2]=aeliminar[1][0]->id;
  p->arestes_per_destruir[3]=aeliminar[1][1]->id;
  p->arestes_per_destruir[4]=acentral->id;
  p->arestes_per_destruir[5]=amantenir[1]->id;
  p->n_arestes_destruir=6;
  
  p->destrueix_elements();
  
  
  
}

matriu aresta::troba_moment_inercia(double x, double y)
{
  matriu resultat;
  double x1, y1, x2, y2;
  double temp;
  
  //Trobem la posicio dels vertexs de l'aresta respecte l'origen de coordenades definit per (x,y)
  x1=(v[0]->x)-x;
  y1=(v[0]->y)-y;
  x2=(v[1]->x)-x;
  y2=(v[1]->y)-y;
  
  //Calculem el moment d'inercia de l'aresta respecte l'origen de coordenades definit per (x,y)
  temp=sqrt(1. + (((y1 - y2)*(y1 - y2))/((x1 - x2)*(x1 - x2))));
  
  resultat.m[0][0]=-(temp*(x1 - x2)*((y1*y1) + y1*y2 + (y2*y2)))/3.;
  resultat.m[1][1]=(temp*(-(x1*x1*x1) + (x2*x2*x2)))/3.;
  resultat.m[0][1]=(temp*(x1 - x2)*(x1*(2*y1 + y2) + x2*(y1 + 2*y2)))/6.;
  
  //Comprova que els signes siguin els correctes
  if(resultat.m[0][0]<0.)
  {
    resultat.m[0][0]=-resultat.m[0][0];
    resultat.m[1][1]=-resultat.m[1][1];
    resultat.m[0][1]=-resultat.m[0][1];
  }
  resultat.m[1][0]=resultat.m[0][1];
  
  if(resultat.m[1][1]<0.)
  {
    std::cout << "Error: when computing the inertia tensor!!" << std::endl;
    exit(1);
  }
  
  return resultat;
}

void aresta::copia_des_de_aresta(aresta *aoriginal)
{
  int i;
  
  for(i=0; i<2; i++)
  {
    x[i] = aoriginal->x[i];
    y[i] = aoriginal->y[i];
  }
  
  lambdah = aoriginal->lambdah;
  lambdav = aoriginal->lambdav;
  templh2lv = aoriginal->templh2lv;
  templv2lh = aoriginal->templv2lh;
  k_gamma_aresta = aoriginal->k_gamma_aresta;
  canvi = aoriginal->canvi;
  l = aoriginal->l;
  
  for(i=0; i<MAXIM_CELULES_A; i++)
  {
    c[i] = aoriginal->c[i];
    cp[i] = aoriginal->cp[i];
  }
  
  for(i=0; i<MAXIM_VERTEXS_A; i++)
  {
    v[i] = aoriginal->v[i];
  }
}

