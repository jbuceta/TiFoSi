#include <iostream>
using namespace std;
#include "main.h"


void celula::inicia_celula(poblacio *pt, int i, int tipus_celula, double a0, int nc, int nv, int na)
{
  double temp;
  
  p=pt;
  id=i;
  ctype=tipus_celula;
  
  
  k_gamma=p->k_gamma[ctype];
  k_kappa=p->k_kappa[ctype];
  k_force_x=p->k_force_x[ctype];
  k_force_y=p->k_force_y[ctype];
  
  
  
  
  temp=(p->dispersio_cicle[ctype]*2.*MIG_CICLE)+((1.-p->dispersio_cicle[ctype])*expdev(2.*MIG_CICLE, &(p->llavor)));
  for(int idx=0; idx < p->nfases[ctype]; idx++)
  {
    duracio_cicle[idx] = temp*p->proporcio_cicle[ctype][idx];
  }
  posicio_cicle = 0.;
  fase = 0;
  pas_cicle=p->pas_cicle[ctype];
  
  if(p->proporcio_cicle[ctype][0]<=0.)
  {
    allow_change_phase = 0;
    area_growth = 0;
  }
  else
  {
    allow_change_phase = 1;
    
    area_growth=(p->area0_final[ctype][0]-p->area0[ctype][0])/duracio_cicle[0];
  }
  
  area0=a0;
  
  ncellvertexes=nv;                    //Numero de vertexs de la celula
  narestes=na;                    //Numero d'arestes de la celula
  neighboringcells=nc;                    //Numero de celules veines, que es igual al numero d'arestes.
  
  proteines.c=this;
  proteines.inicia_constants(p->stage);
}

void celula::recalcula_fase()
{
  pas_cicle=p->pas_cicle[ctype];
  
  if(p->proporcio_cicle[ctype][fase]<=0.)
  {
    allow_change_phase = 0;
    area_growth = 0;
  }
  else
  {
    allow_change_phase = 1;
    
    area_growth=(p->area0_final[ctype][fase]-p->area0[ctype][fase])/duracio_cicle[fase];
  }
  
}

void celula::actualitza_constants()
{
  k_gamma=p->k_gamma[ctype];
  k_kappa=p->k_kappa[ctype];
  k_force_x=p->k_force_x[ctype];
  k_force_y=p->k_force_y[ctype];
}

void celula::introdueix_vertex(vertex *vi, vertex *v1, vertex *v2)
{
  int i;
  int flag;
  int posicio;
  
  if(id!=-1)
  {
    if((v1!=&p->vertex_buit)&&(v2!=&p->vertex_buit))
    {
      for(i=0;i<ncellvertexes;i++)
      {
        if(v[i]==vi)
        {
          std::cout << endl << "*******************************" << endl;
          std::cout << endl << "Error: the vertex " << vi->id << " to insert already exist at the cell " << id << " that has " << ncellvertexes << " vertexes." << endl;
          std::cout << endl << "*******************************" << endl;
          exit(1);                        
        }
      }
      
      flag=0;
      i=0;
      do
      {
        if((v[i%ncellvertexes]==v1)||(v[i%ncellvertexes]==v2))
        {
          if((v[(i+1)%ncellvertexes]==v1)||(v[(i+1)%ncellvertexes]==v2))
          {
            posicio=(i+1)%ncellvertexes;
            flag=1;
          }
          else if((v[(i+ncellvertexes-1)%ncellvertexes]==v1)||(v[(i+ncellvertexes-1)%ncellvertexes]==v2))
          {
            posicio=i%ncellvertexes;
            flag=1;
          }
          else
          {
            std::cout << endl << "*******************************" << endl;
            std::cout << endl << "Error: v1 and v2 are not joined at the cell." << endl;
            std::cout << endl << "*******************************" << endl;
            exit(1);
          }
        }
        i++;
      }while((flag<1)&&(i<ncellvertexes));
      
      if(flag==0)
      {
        std::cout << endl << "*******************************" << endl;
        std::cout << endl << "Error: the vertexes where I have to insert the new vertex in between do not exist." << endl;
        std::cout << endl << "*******************************" << endl;
        exit(1);
      }
      else if(flag==1)
      {
        for(i=ncellvertexes; i>posicio; i--)
        {
          v[i]=v[i-1];
        }
        v[posicio]=vi;
        ncellvertexes++;
      }
      else
      {
        std::cout << endl << "*******************************" << endl;
        std::cout << endl << "Error: paranormal activity (sometimes I see dead people...)." << endl;
        std::cout << endl << "*******************************" << endl;
        exit(1);
      }
    }
    else
    {
      flag=0;
      for(i=0;i<ncellvertexes;i++)
      {
        if(v[i]==vi)
        {
          flag=1;
        }
      }
      if(flag==0)
      {
        v[ncellvertexes]=vi;
        ncellvertexes++;
      }
    }
  }
}

void celula::calcula_propietats()
{
  if(id!=-1)
  {
    troba_arestes();
    troba_celules_veines();
    calcula_area();
    calcula_perimetre();
    calcula_kappa_area_area0();
  }
}

void celula::calcula_kappa_area_area0()
{
  kappa_area_area0=k_kappa*(area-area0);
}

void celula::calcula_perimetre()
{
  int i;
  
  perimeter=a[0]->l;
  
  for(i=1;i<narestes;i++)
  {
    perimeter=perimeter+a[i]->l;
  }
  
}

void celula::escriu_celula(ofstream &arxiu)
{
  int i;
  
  
  arxiu << c_track_id << " " << ctype << " " << ncellvertexes << " " << area << " ";
  proteines.escriu_especies(arxiu);
  
  #if GUARDA_VEINS_CELULA==True
  calcula_centre();
  arxiu << x << " " << y << " ";
  for(i=0;i<ncellvertexes;i++){
    arxiu << c[i]->id << " ";
  }
  #endif
  
  for(i=0;i<ncellvertexes;i++){
    arxiu << v[i]->x << " " << v[i]->y << " ";
  }
  
  arxiu << std::endl;
  
}

void celula::calcula_centre()
{
  int i;
  double temp, temp2;
  
  //Sembla que ja esta be. Comprovar-ho!!! Mirar la wikipedia: centroid irregular polygon.
  //Podria ser que hi hagues un error degut a que l'area ha perdut el signe.
  
  if(id!=-1)
  {
    temp=6.*area_centre/2.;
    
    x=0.;
    y=0.;
    for(i=0; i<ncellvertexes-1; i++)
    {
      temp2=(((v[i]->x)*(v[i+1]->y))-((v[i+1]->x)*(v[i]->y)));
      x=x+((v[i]->x)+(v[i+1]->x))*temp2;
      y=y+((v[i]->y)+(v[i+1]->y))*temp2;
    }
    temp2=((v[ncellvertexes-1]->x)*(v[0]->y))-((v[0]->x)*(v[ncellvertexes-1]->y));
    x=x+((v[ncellvertexes-1]->x)+(v[0]->x))*temp2;
    y=y+((v[ncellvertexes-1]->y)+(v[0]->y))*temp2;
    x=x/temp;
    y=y/temp;
  }
  
}

void celula::calcula_area()
{
  int i;
  
  if(id!=-1)
  {
    area_centre=0.;
    for(i=0; i<ncellvertexes-1; i++)
    {
      area_centre=area_centre+((v[i]->x)*(v[i+1]->y))-((v[i+1]->x)*(v[i]->y));
    }
    area_centre=area_centre+((v[ncellvertexes-1]->x)*(v[0]->y))-((v[0]->x)*(v[ncellvertexes-1]->y));
    area=fabs(area_centre/2.);
  }
}

void celula::calcula_area_dinamica()
{
  int i;
  
  area_centre=0.;
  for(i=0; i<ncellvertexes-1; i++)
  {
    area_centre=area_centre+((v[i]->x)*(v[i+1]->y))-((v[i+1]->x)*(v[i]->y));
  }
  area_centre=area_centre+((v[ncellvertexes-1]->x)*(v[0]->y))-((v[0]->x)*(v[ncellvertexes-1]->y));
  area=fabs(area_centre/2.);
  
  area0=area0+(area_growth*pas_cicle);
  area0=(area0<0.)?0.:area0;
  
  if(
    ((area>((p->area0_final[ctype][fase])*(p->reldiv[ctype][fase]))) || (area_growth<0.))
    &&
    (posicio_cicle>duracio_cicle[fase])
    &&
    (allow_change_phase==1)
  )
  {
    if(fase == (p->nfases[ctype] - 1))
    {
      divideix_celula();
    }
    else
    {
      fase++;
      posicio_cicle = 0.;
      area0=p->area0[ctype][fase];
      if(p->proporcio_cicle[ctype][fase]<=0.)
      {
        allow_change_phase = 0;
        area_growth=0;
      }
      else
      {
        allow_change_phase = 1;
        
        area_growth=(p->area0_final[ctype][fase]-p->area0[ctype][fase])/duracio_cicle[fase];
      }
    }
  }
}

void celula::calcula_rellotge()
{
  posicio_cicle=posicio_cicle+pas_cicle;
  
  
}

void celula::troba_celules_veines()
{
  int i;
  
  if(id!=-1)
  {
    neighboringcells=narestes;
    
    for(i=0; i<narestes; i++)
    {
      c[i]=(this!=a[i]->c[0]?a[i]->c[0]:a[i]->c[1]);
    }
  }
  
}

void celula::escriu_informacio_celula(ofstream &arxiu)
{
  int i;
  
  arxiu << "Celula: " << id << " amb area: " << area << " i area0: " << area0 << std::endl;
  arxiu << "que te un perimeter de: " << perimeter << " i centre (x,y): (" << x << ", " << y << ")" << std::endl;
  arxiu << "i kappa_area_area0: " << kappa_area_area0 << std::endl;
  arxiu << "-----------" << std::endl;
  arxiu << "      Numero de vertexs: " << ncellvertexes << std::endl;
  for(i=0; i<ncellvertexes; i++)
  {
    arxiu << "              Vertex: " << v[i]->id << " (" << v[i]->x << ", " << v[i]->y << ")" << std::endl;
  }
  arxiu << "      Numero d'arestes veines: " << narestes << std::endl;
  for(i=0; i<narestes; i++)
  {
    arxiu << "              Aresta: " << a[i]->id << " que fa frontera amb: " << (a[i]->c[0]->id==id?a[i]->c[1]->id:a[i]->c[0]->id) << std::endl;
    arxiu << "              formada pels vertexs: " << a[i]->v[0]->id << " i " << a[i]->v[1]->id << std::endl;
    arxiu << "              i que te una longitud de: " << a[i]->l << std::endl;
  }
  arxiu << "      Numero de celules: " << neighboringcells << std::endl;
  for(i=0; i<neighboringcells; i++)
  {
    arxiu << "              Celula: " << c[i]->id << " amb " << c[i]->neighboringcells << " celules veines " << std::endl;
  }
}

void celula::elimina_vertex(vertex *ve)
{
  int i;
  int posicio;
  
  if(id!=-1)
  {
    posicio=-1;
    for(i=0;i<ncellvertexes;i++)
    {
      if(v[i]==ve)
      {
        posicio=i;
      }
    }
    if(posicio!=-1)
    {
      posicio++;
      for(i=posicio;i<ncellvertexes;i++)
      {
        v[i-1]=v[i];
      }
      ncellvertexes--;
      narestes=ncellvertexes;
    }
    else
    {
      std::cout << endl << "*******************************" << endl;
      std::cout << endl << "Error: I tried to kill a vertex that does not exist." << endl;
      std::cout << endl << "*******************************" << endl;
      exit(1);
    }
  }
}

void celula::troba_arestes()
{
  int i,j;
  int flag;
  
  if(id!=-1)
  {
    for(i=0; i<ncellvertexes; i++)
    {
      flag=0;
      for(j=0; j<v[i]->narestes; j++)
      {
        if(v[i]->a[j]->v[0]==v[(i+1)%ncellvertexes])
        {
          flag++;
          a[i]=v[i]->a[j];
        }
        else if(v[i]->a[j]->v[1]==v[(i+1)%ncellvertexes])
        {
          flag++;
          a[i]=v[i]->a[j];
        }
      }
      if(flag!=1)
      {
        
        std::cout << endl << "*******************************" << endl;
        std::cout << endl << "Error: when calling the function that looks for the edges that belong to cell " << id << ", I found that the vertexes do not have well defined edges: " << flag << endl;
        std::cout << "In particular, looking for the edge that connect the vertex " << v[i]->id << " to vertex " << v[(i+1)%ncellvertexes]->id << "."<< std::endl;
        std::cout << "At vertex " << v[i]->id << " there are " << v[i]->narestes << " edges:" << std::endl;
        for(j=0; j<v[i]->narestes; j++)
        {
          std::cout << "\tVertex: " << v[i]->a[j]->id << " that connects to vertex: " << v[i]->v[j]->id << std::endl;
          std::cout << "\t\tthat connects to the vertexes: " << v[i]->a[j]->v[0]->id << " and " << v[i]->a[j]->v[1]->id << std::endl;
          std::cout << "\t\tthat separates the cells: " << v[i]->a[j]->c[0]->id << " and " << v[i]->a[j]->c[1]->id << std::endl;
        }
        std::cout << "At vertex " << v[(i+1)%ncellvertexes]->id << " there are " << v[i]->narestes << " edges:" << std::endl;
        for(j=0; j<v[(i+1)%ncellvertexes]->narestes; j++)
        {
          std::cout << "\tVertex: " << v[(i+1)%ncellvertexes]->a[j]->id << " that connect to vertex: " << v[(i+1)%ncellvertexes]->v[j]->id << std::endl;
          std::cout << "\t\tthat connects to the vertexes: " << v[(i+1)%ncellvertexes]->a[j]->v[0]->id << " and " << v[(i+1)%ncellvertexes]->a[j]->v[1]->id << std::endl;
          std::cout << "\t\tthat separated the cells: " << v[(i+1)%ncellvertexes]->a[j]->c[0]->id << " and " << v[(i+1)%ncellvertexes]->a[j]->c[1]->id << std::endl;
        }
        std::cout << endl << "*******************************" << endl;
        exit(1);
      }
    }
    narestes=ncellvertexes;
    troba_celules_veines();
  }
}

vertex *celula::torna_vertex_anterior(vertex *vc)
{
  int i;
  vertex *resultat;
  
  resultat=&p->vertex_buit;
  for(i=0; i<ncellvertexes; i++)
  {
    if(v[i]==vc)
    {
      resultat=v[(i+ncellvertexes-1)%ncellvertexes];
    }
  }
  
  if(resultat==&p->vertex_buit)
  {
    std::cout << endl << "*******************************" << endl;
    std::cout << endl << "Error: I was asked for a vertex previous to a vertex that does not exist." << endl;
    std::cout << endl << "*******************************" << endl;
    exit(1);
  }
  
  return resultat;
}

vertex *celula::torna_vertex_posterior(vertex *vc)
{
  int i;
  vertex *resultat;
  
  resultat=&p->vertex_buit;
  for(i=0; i<ncellvertexes; i++)
  {
    if(v[i]==vc)
    {
      resultat=v[(i+1)%ncellvertexes];
    }
  }
  
  if(resultat==&p->vertex_buit)
  {
    std::cout << endl << "*******************************" << endl;
    std::cout << endl << "Error: I was asked for a vertex previous to a vertex that does not exist." << endl;
    std::cout << endl << "*******************************" << endl;
    exit(1);
  }
  
  return resultat;
}

void celula::divideix_celula()
{
  int i;
  
  p->matriu_c[p->n_matriu_c].inicia_celula(p,p->n_matriu_c,ctype,p->area0[ctype][0]);
  inicia_celula(p,id,ctype,p->area0[ctype][0],neighboringcells,ncellvertexes,narestes);
  
  
  p->matriu_c[p->n_matriu_c].c_track_id = c_track_id + "-" + toString( p->c_track_idx);
  p->c_track_idx++;
  c_track_id = c_track_id + "-" + toString( p->c_track_idx);
  p->c_track_idx++;
  
  
  
  proteines.divideix_celula(&(p->matriu_c[p->n_matriu_c]));
  
  
  
  p->matriu_v[p->n_matriu_v].inicia_vertex(p,p->n_matriu_v,0.,0.);
  p->matriu_v[p->n_matriu_v+1].inicia_vertex(p,p->n_matriu_v+1,0.,0.);
  
  p->matriu_a[p->n_matriu_a].p=p;
  p->matriu_a[p->n_matriu_a].id=p->n_matriu_a;
  p->matriu_a[p->n_matriu_a].canvi=0;
  p->matriu_a[p->n_matriu_a+1].p=p;
  p->matriu_a[p->n_matriu_a+1].id=p->n_matriu_a+1;
  p->matriu_a[p->n_matriu_a+1].canvi=0;
  p->matriu_a[p->n_matriu_a+2].p=p;
  p->matriu_a[p->n_matriu_a+2].id=p->n_matriu_a+2;
  p->matriu_a[p->n_matriu_a+2].canvi=0;
  
  //std::cout << "celula: " << id << " amb area0: " << area0 << " i area: " << area << std::endl;
  p->flog << p->stage << " " << p->idx_exterior << " " << p->idx_interior << ": cell division. Daughter cells: " << c_track_id << " and " << p->matriu_c[p->n_matriu_c].c_track_id << std::endl;
  
  aux_divideix_celula(&p->matriu_c[p->n_matriu_c],&p->matriu_v[p->n_matriu_v],&p->matriu_v[p->n_matriu_v+1],&p->matriu_a[p->n_matriu_a],&p->matriu_a[p->n_matriu_a+1],&p->matriu_a[p->n_matriu_a+2]);
  
  p->matriu_a[p->n_matriu_a].calcula_lambda();
  p->matriu_a[p->n_matriu_a+1].calcula_lambda();
  p->matriu_a[p->n_matriu_a+2].calcula_lambda();
  
  p->n_matriu_c=p->n_matriu_c+1;
  p->n_matriu_v=p->n_matriu_v+2;
  p->n_matriu_a=p->n_matriu_a+3;
  //Comprovar que no es passa del maxim de celules, arestes o vertexs disponibles.
  for(i=0; i<p->n_matriu_v; i++)
  {
    p->matriu_v[i].calcula_f0();
  }
}

void celula::aux_divideix_celula(celula *cn, vertex *vn1, vertex *vn2, aresta *anc, aresta *an1, aresta *an2)
{
  int i;
  int index;
  int nvv, nvn;           //Variables on guardem el numero de vertexs en la celula vella i nova respectivament.
  double alfa, beta;
  double x0,x1,y0,y1;
  aresta *ai[2];
  int posicio_v[2];
  vertex *vtemp[2];
  vertex *veins_vtemp[2][2];      //Es guardaran els veins de vtemp[] en ordre horari. El primer index es el que relaciona amb vtemp.
  celula *ctemp[2];
  punt punt_temp;
  
  ofstream error;
  
  ai[0]=NULL;
  ai[1]=NULL;
  vtemp[0]=NULL;
  vtemp[1]=NULL;
  veins_vtemp[0][0]=NULL;
  veins_vtemp[0][1]=NULL;
  veins_vtemp[1][0]=NULL;
  veins_vtemp[1][1]=NULL;
  ctemp[0]=NULL;
  ctemp[1]=NULL;
  
  //Comprovem que la celula tingui ben definides les arestes i les celules veines.
  
  calcula_propietats();
  
  //Calculem el centre geometric de la per calcular com sera la divisio.
  
  calcula_centre();
  
  //Definim la direccio del segment que separara la celula en dos.
  
  x0=x;
  y0=y;
  
  punt_temp=calcula_direccio_divisio();
  if(punt_temp.actiu==1)
  {
    x1=punt_temp.x;
    y1=punt_temp.y;
  }
  else
  {
    std::cout << endl << "*******************************" << endl;
    std::cout << endl << "Error: when dividing a cell I could not find a cleavage orientation." << endl;
    std::cout << endl << "*******************************" << endl;
    exit(1);
  }
  
  alfa=(y1)/(x1);
  beta=y0-(x0*alfa);
  
  //Reassignem els punters dels vertexs nous que es crearan amb la divisio a la matriu vtemp[].
  
  vtemp[0]=vn1;
  vtemp[1]=vn2;
  
  //Aquest bucle ens trobara els dos punts on intersecta la recta divisora amb el perimetre de la celula.
  //Tambe ens definira algunes variables.
  //Ens trobara el següent:
  //      vtemp[].x, vtemp[].y            sera la posicio dels nous vertexs.
  //      veins_vtemp[][]                 seran els vertexs veins dels nous vertexs en aquesta celula.
  //      ai[]                            sera un punter a les arestes on es troben les interseccions. Estaran en ordre horari.
  //      ctemp[]                         seran punters a les celules veines, que separa l'aresta ai[] de la celula actual.
  //      posicio_v[]                     guarda l'index (dins de la celula) del primer v[] que pertany.a l'aresta ai[] en ordre horari.
  
  index=0;
  for(i=0; i<narestes; i++)
  {
    punt_temp=a[i]->comprova_si_intersecta_amb_recta(alfa, beta);
    if(punt_temp.actiu==1)
    {
      if(index==2)
      {
        std::cout << endl << "*******************************" << endl;
        std::cout << endl << "Error: when diving a cell the cleavage plane intersects more than two edges.." << endl;
        std::cout << endl << "*******************************" << endl;
        exit(1);
      }
      else
      {
        vtemp[index]->x=punt_temp.x;
        vtemp[index]->y=punt_temp.y;
        veins_vtemp[index][0]=v[i];
        veins_vtemp[index][1]=v[(i+1)%narestes];
        ai[index]=a[i];
        ctemp[index]=((a[i]->c[0]!=this)?a[i]->c[0]:a[i]->c[1]);
        posicio_v[index]=i;
        index++;
      }
    }
  }
  
  if(index<2)
  {
    std::cout << endl << "*******************************" << endl;
    std::cout << endl << "Error: I could not find the intersection points when dividing a cell: " << id << endl;
    std::cout << "with its center at: " << "(" << x << ", " << y << ")" << endl;
    std::cout << endl << "I found " << index << " intersections" << endl;
    for(i=0; i<index; i++)
    {
      std::cout << "\tAt the edge: " << ai[i]->id << " at the point (" << vtemp[i]->x << ", " << vtemp[i]->y << ")" << endl;
    }
    std::cout << "there are " << narestes << " edges." << endl;
    std::cout << "The rest of the edges are:" << endl;
    for(i=0; i<narestes; i++)
    {
      std::cout << "\tEdges: " << a[i]->id << " with a length: " << a[i]->l <<  endl;
      std::cout << "\ti Vertexes: " << a[i]->v[0]->id << " (" << a[i]->v[0]->x << ", " << a[i]->v[0]->y << ") i " << a[i]->v[1]->id << " (" << a[i]->v[1]->x << ", " << a[i]->v[1]->y << ")" << endl;
    }
    std::cout << endl << "*******************************" << endl;
    exit(1);
  }
  
  
  
  
  //Introdueix els vertexs nous vtemp[] en les celules veines que es trobaran afectades per la divisio: ctemp[].
  
  if(ctemp[0]!=&p->celula_buida)
  {
    ctemp[0]->introdueix_vertex(vtemp[0],ai[0]->v[0],ai[0]->v[1]);          //Calcula arestes, mes avall.
  }
  if(ctemp[1]!=&p->celula_buida)
  {
    ctemp[1]->introdueix_vertex(vtemp[1],ai[1]->v[0],ai[1]->v[1]);          //Calcula arestes, mes avall.
  }
  
  //Els següents bucles reparteixen els vertexs antics entre la celula que es divideix i la celula nova.
  //      Aquest primer assigna els vertexs i arestes que van entre v[posicio_v[1]+1] i v[posicio>v[0]] a la celula nova.
  //      Tambe canvia, ens els vertexs nous, la informacio de les celules veines. Canvia l'actual per les celules noves.
  
  index=posicio_v[1];
  nvn=0;
  cn->ncellvertexes=0;
  do
  {
    index=(index+1)%ncellvertexes;
    cn->v[nvn]=v[index];
    cn->v[nvn]->canvia_celules_veines(this,cn);
    cn->a[nvn]=a[index];
    cn->a[nvn]->canvia_celules_veines(this, cn);
    cn->ncellvertexes=cn->ncellvertexes+1;
    nvn++;
  }while(index!=posicio_v[0]);
  
  cn->a[nvn-1]->canvia_conexio(cn->v[nvn-1], veins_vtemp[0][1], vtemp[0]);
  
  cn->v[cn->ncellvertexes]=vtemp[0];
  cn->a[cn->ncellvertexes]=anc;
  cn->ncellvertexes=cn->ncellvertexes+1;
  cn->v[cn->ncellvertexes]=vtemp[1];
  cn->a[cn->ncellvertexes]=an2;
  cn->ncellvertexes=cn->ncellvertexes+1;
  
  cn->narestes=cn->ncellvertexes;
  
  //      El segon bucle realigna els indexs dels vertexs que s'han qeudat en aquesta celula, de tal forma que el primer tingui l'index 0
  //      Tambe redefineix la variable nvertexs, donant-li el valor del numero de vertexs que te de moment la celula antiga.
  
  nvv=ncellvertexes-nvn;
  nvn=cn->ncellvertexes;
  
  for(i=0; i<nvv; i++)
  {
    v[i]=v[i+posicio_v[0]+1];
    a[i]=a[i+posicio_v[0]+1];
    if((i+posicio_v[0]+1)>(ncellvertexes-1))
    {
      std::cout << endl << "*******************************" << endl;
      std::cout << endl << "Error: I found a problem when distributing the vertexes when dividing a cell." << endl;
      std::cout << endl << "*******************************" << endl;
      exit(1);
    }
  }
  ncellvertexes=nvv;
  
  a[ncellvertexes-1]->canvia_conexio(v[ncellvertexes-1], veins_vtemp[1][1], vtemp[1]);
  
  v[ncellvertexes]=vtemp[1];
  a[ncellvertexes]=anc;
  ncellvertexes++;
  v[ncellvertexes]=vtemp[0];
  a[ncellvertexes]=an1;
  ncellvertexes++;
  
  nvv=ncellvertexes;
  narestes=ncellvertexes;
  
  
  //Arreglem les celules veines i els vertexs i arestes conectades dels vertexs afectats:
  
  //      Dels dos nous vertexs que s'han creat:
  
  vtemp[0]->v[0]=vtemp[1];
  vtemp[0]->a[0]=anc;
  vtemp[0]->v[1]=veins_vtemp[0][0];
  vtemp[0]->a[1]=ai[0];
  vtemp[0]->v[2]=veins_vtemp[0][1];
  vtemp[0]->a[2]=an1;
  vtemp[0]->narestes=3;
  
  vtemp[0]->c[0]=cn;
  vtemp[0]->c[1]=this;
  vtemp[0]->ncelules=2;
  if(ctemp[0]!=&p->celula_buida)
  {
    vtemp[0]->c[2]=ctemp[0];
    vtemp[0]->ncelules=vtemp[0]->ncelules+1;
  }
  
  vtemp[1]->v[0]=vtemp[0];
  vtemp[1]->a[0]=anc;
  vtemp[1]->v[1]=veins_vtemp[1][0];
  vtemp[1]->a[1]=ai[1];
  vtemp[1]->v[2]=veins_vtemp[1][1];
  vtemp[1]->a[2]=an2;
  vtemp[1]->narestes=3;
  
  vtemp[1]->c[0]=cn;
  vtemp[1]->c[1]=this;
  vtemp[1]->ncelules=2;
  if(ctemp[1]!=&p->celula_buida)
  {
    vtemp[1]->c[2]=ctemp[1];
    vtemp[1]->ncelules=vtemp[1]->ncelules+1;
  }
  
  
  //      Dels vertexs que ja existien.
  
  veins_vtemp[0][0]->canvia_conexio_amb_vertex(veins_vtemp[0][1], ai[0], vtemp[0], ai[0]);
  veins_vtemp[0][1]->canvia_conexio_amb_vertex(veins_vtemp[0][0], ai[0], vtemp[0], an1);
  veins_vtemp[0][0]->canvia_celules_veines(this,cn);
  
  veins_vtemp[1][0]->canvia_conexio_amb_vertex(veins_vtemp[1][1], ai[1], vtemp[1], ai[1]);
  veins_vtemp[1][1]->canvia_conexio_amb_vertex(veins_vtemp[1][0], ai[1], vtemp[1], an2);
  veins_vtemp[1][1]->canvia_celules_veines(this,cn);
  
  //Posem la informacio dels vertexs i les celules a les arestes noves.
  
  anc->v[0]=vtemp[1];
  anc->v[1]=vtemp[0];
  anc->c[0]=this;
  anc->c[1]=cn;
  
  an1->v[0]=vtemp[0];
  an1->v[1]=veins_vtemp[0][1];
  an1->c[0]=this;
  an1->c[1]=ctemp[0];
  
  an2->v[0]=vtemp[1];
  an2->v[1]=cn->v[0];
  an2->c[0]=cn;
  an2->c[1]=ctemp[1];
  
  //Trobem les noves arestes i celules veines de les celules afectades per la divisio. Tambe calculem les noves propietats d'arestes i celules.
  
  if(ctemp[0]!=&p->celula_buida)
  {
    ctemp[0]->calcula_propietats();
    for(i=0; i<ctemp[0]->neighboringcells; i++)
    {
      ctemp[0]->c[i]->calcula_propietats();
    }
    
  }
  
  if(ctemp[1]!=&p->celula_buida)
  {
    ctemp[1]->calcula_propietats();
    for(i=0; i<ctemp[1]->neighboringcells; i++)
    {
      ctemp[1]->c[i]->calcula_propietats();
    }
    
  }
  
  calcula_propietats();
  for(i=0; i<neighboringcells; i++)
  {
    c[i]->calcula_propietats();
  }
  
  
  cn->calcula_propietats();
  for(i=0; i<cn->neighboringcells; i++)
  {
    cn->c[i]->calcula_propietats();
  }
  
  
  vtemp[0]->troba_vertexs_veins();
  vtemp[1]->troba_vertexs_veins();
  veins_vtemp[0][0]->troba_vertexs_veins();
  veins_vtemp[0][1]->troba_vertexs_veins();
  veins_vtemp[1][0]->troba_vertexs_veins();
  veins_vtemp[1][1]->troba_vertexs_veins();
  
  ai[0]->calcula_longitud();
  ai[0]->troba_celules_puntes_aresta();
  ai[1]->calcula_longitud();
  ai[1]->troba_celules_puntes_aresta();
  anc->calcula_longitud();
  anc->troba_celules_puntes_aresta();
  an1->calcula_longitud();
  an1->troba_celules_puntes_aresta();
  an2->calcula_longitud();
  an2->troba_celules_puntes_aresta();
  
  
  
  
}

punt celula::calcula_direccio_divisio()
{
  punt punt_temp;
  double angle_divisio;
  double angle_divisio_hertwig;
  double angle_divisio_shift;
  double angle_divisio_shift_dispersion;
  double temp;
  double fractpart;
  double intpart;
  int i;
  matriu tensor_inercia_celula;
  const double PI = 3.141592654;
  const double PI_half = 1.570796327;
  const double two_PI = 6.283185308;
  const double three_PI_half = 4.71238898;
  
  calcula_centre();
  
  tensor_inercia_celula=a[0]->troba_moment_inercia(x,y);
  for(i=1; i<narestes; i++)
  {
    tensor_inercia_celula=tensor_inercia_celula+a[i]->troba_moment_inercia(x,y);
  }
  
  tensor_inercia_celula.troba_valors_i_vectors_propis();
  tensor_inercia_celula.ordena_valors_i_vectors_propis();
  
  //Calculem l'angle de divisio tenint en compte que la funcio asin() torna un resultat entre (pi/2,-pi/2)
  temp=sqrt(tensor_inercia_celula.v1.x*tensor_inercia_celula.v1.x+tensor_inercia_celula.v1.y*tensor_inercia_celula.v1.y);
  temp=tensor_inercia_celula.v1.y/temp;
  if(tensor_inercia_celula.v1.x<0)
  {
    temp=-temp;
  }
  angle_divisio=asin(temp);
  angle_divisio_hertwig=angle_divisio;
  //angle_divisio_hertwig corresponds to perpendicular to the longest cell axis
  
  angle_divisio = angle_divisio + p->divisionshift[ctype];
  angle_divisio_shift=angle_divisio;
  if(angle_divisio_shift>PI_half){angle_divisio_shift=angle_divisio_shift - PI;};
  if(angle_divisio_shift<-PI_half){angle_divisio_shift=angle_divisio_shift + PI;};
  //angle_divisio_shift corresponds to perpendicular to the longest cell axis plus the shift in the configuration file: 0 if hertwig and pi/2 if perpendicular (parallel to the longest cell axis)
  
  do
  {
    temp=box_muller(angle_divisio, p->divisiondispersion[ctype], &(p->llavor));
  }while((temp<(angle_divisio-p->divisiondispersionlimit[ctype]))||(temp>(angle_divisio+p->divisiondispersionlimit[ctype])));
  angle_divisio=temp;
  angle_divisio_shift_dispersion=angle_divisio;
  fractpart=modf(angle_divisio_shift_dispersion/two_PI,&intpart);
  angle_divisio_shift_dispersion=angle_divisio_shift_dispersion-two_PI*intpart;//This restrict the angle to a 0-2 pi range
  if(three_PI_half>angle_divisio_shift_dispersion && angle_divisio_shift_dispersion>PI_half){angle_divisio_shift_dispersion= angle_divisio_shift_dispersion - PI;};
  if(two_PI>angle_divisio_shift_dispersion && angle_divisio_shift_dispersion>three_PI_half){angle_divisio_shift_dispersion= angle_divisio_shift_dispersion - two_PI;};
  if(-three_PI_half<angle_divisio_shift_dispersion && angle_divisio_shift_dispersion<-PI_half){angle_divisio_shift_dispersion= angle_divisio_shift_dispersion + PI;};
  if(-two_PI<angle_divisio_shift_dispersion && angle_divisio_shift_dispersion<-three_PI_half){angle_divisio_shift_dispersion= angle_divisio_shift_dispersion + two_PI;};
  //angle_divisio_shift_dispersion corresponds to perpendicular to the longest cell axis plus the shift in the configuration file plus the random component: this is the actual division angle
  
  punt_temp.actiu=1;
  punt_temp.x=cos(angle_divisio);
  punt_temp.y=sin(angle_divisio);
  
  
  //  p->direccio_divisio << p->stage << " " << p->idx_exterior << " " << p->idx_interior << " " << id << " " << c_track_id << " " << p->matriu_c[p->n_matriu_c].c_track_id << " " << tipus << " " << angle_divisio << " " << x << " " << y << " " << area << " " << p->n_matriu_c << " " << narestes;
  p->direccio_divisio << p->stage << " " << p->idx_exterior << " " << p->idx_interior << " " << id << " " << c_track_id << " " << p->matriu_c[p->n_matriu_c].c_track_id << " " << ctype << " " << angle_divisio_shift_dispersion << " " << angle_divisio_hertwig << " " << angle_divisio_shift << " " << x << " " << y << " " << area << " " << p->n_matriu_c << " " << narestes;

  for(i=0; i<narestes; i++)
  {
    p->direccio_divisio << " " << v[i]->x << " " << v[i]->y;
  }
  p->direccio_divisio << std::endl;
  
  return punt_temp;
  
}

void celula::cambia_referencia_celula(celula *dv, celula *dn)
{
  int i, j;
  
  if((id!=-1)&&(dv!=dn))
  {
    for(i=0; i<ncellvertexes; i++)
    {
      for(j=0; j<v[i]->ncelules; j++)
      {
        if(v[i]->c[j]==dv)
        {
          v[i]->c[j]=dn;
        }
      }
    }
    for(i=0; i<narestes; i++)
    {
      for(j=0; j<2; j++)
      {
        if(a[i]->c[j]==dv)
        {
          a[i]->c[j]=dn;
        }
      }
      for(j=0; j<2; j++)
      {
        if(a[i]->cp[j]==dv)
        {
          a[i]->cp[j]=dn;
        }
      }
    }
    for(i=0; i<neighboringcells; i++)
    {
      for(j=0; j<c[i]->neighboringcells; j++)
      {
        if(c[i]->c[j]==dv)
        {
          c[i]->c[j]=dn;
        }
      }
    }
  }
  
  proteines.c=dn;
}



void celula::actualitza_speed_1()
{
  if(ctype == 2)
  {
    pas_cicle = 5.0;
  }
  else if(ctype == 1)
  {
    pas_cicle = 5.0;
  }
  else
  {
    std::cout << endl << "*******************************" << endl;
    std::cout << endl << "Error. There is a bug in the cycle speed computation." << endl;
    std::cout << endl << "*******************************" << endl;
    exit(1);
  }
}

