#include <iostream>
using namespace std;

#include "main.h"
#include "auxiliar.h"

<poblacio_cpp_1>

void poblacio::guarda_dades()
{
  int i;
  
  forces << n_matriu_v  << " " << stage << " " << std::endl;
  energia << n_matriu_v  << " " << stage << " " << std::endl;
  for(i=0; i<n_matriu_v; i++)
  {
    
    matriu_v[i].escriu_informacio_forces(forces);
    matriu_v[i].escriu_informacio_energia(energia);
  }
  
  celules << n_matriu_c  << " " << stage << " " << std::endl;
  for(i=0; i<n_matriu_c; i++)
  {
    
    matriu_c[i].escriu_celula(celules);
  }
  
  ftime << stage << " " << idx_exterior << " Partial: " << double( clock() - partial_start_time ) / (double)CLOCKS_PER_SEC << "s Total: " << double( clock() - start_time ) / (double)CLOCKS_PER_SEC << "s" << std::endl;
  ftime.flush();
  partial_start_time = clock();
}

void poblacio::guarda_dades_final(std::string stageIdx)
{
  int i;
  ofstream forces_final;
  ofstream celules_final;
  
  std::string temp;
  temp = "dforces_final_s" + stageIdx + ".dat";
  forces_final.open(temp.c_str());
  temp = "dcells_final_s" + stageIdx + ".dat";
  celules_final.open(temp.c_str());
  
  forces_final << n_matriu_v  << " " << stage << " " << std::endl;
  for(i=0; i<n_matriu_v; i++)
  {
    matriu_v[i].escriu_informacio_forces(forces_final);
  }
  
  celules_final << n_matriu_c  << " " << stage << " " << std::endl;
  for(i=0; i<n_matriu_c; i++)
  {
    
    matriu_c[i].escriu_celula(celules_final);
  }
    
  forces_final.close();
  celules_final.close();
}

void poblacio::crea_poblacio()
{
  n_matriu_c=0;
  n_matriu_v=0;
  n_matriu_a=0;
  
  n_celules_destruir=0;
  n_vertexs_destruir=0;
  n_arestes_destruir=0;
  
  llavor=500;
  
  
  stage=1;
  timepos=0.;
  defineix_constants_stage(1);
  c_track_idx = 0;
  
  start_time = clock();
  partial_start_time = clock();
  
  //Definim les constants del potencial
  /*k_kappa=0.5;          //Equival a k_kappa/2 de l'article;*/
  
  
  //Definim com estaran distribuides les celules de la població
  
  <poblacio_cpp_regions>
  
  //Creem la població
  
  matriu_c=new celula[N_CELULES];
  matriu_v=new vertex[N_VERTEXS];
  matriu_a=new aresta[N_ARESTES];
  
  if((matriu_c==0)||(matriu_v==0)||(matriu_a==0))
  {
    std::cout << endl << "*******************************" << endl;
    //    std::cout << endl << "No hi ha prou memòria per crear la població" << endl;
    std::cout << endl << "Error: there is not enough memory to create the tissue." << endl;
    std::cout << endl << "*******************************" << endl;
    exit(1);
  }
  
  forces.open("dforces.dat");
  energia.open("denergy.dat");
  celules.open("dcells.dat");
  
  constantspotencial.open("dconstants.dat");
  
  direccio_divisio.open("divisions.dat");
  flog.open("tifosi.log");
  ftime.open("time.dat");
  
  crea_celula_buida_vertex_buit_i_aresta_buida();
  
  <poblacio_cpp_celules_inicials>
  
  //fes_vertexs_random();
  fes_fase_random();
  
  if(n_matriu_c != 0)
  {
    ofstream ordre_especies;
    ordre_especies.open("protein_order.dat");
    ordre_especies << matriu_c[0].proteines.especies_en_ordre() << std::endl;
    ordre_especies.close();
  }
  else
  {
    std::cout << endl << "*******************************" << endl;
    std::cout << endl << "Error: The number of cells in the tissue is 0." << endl;
    std::cout << endl << "*******************************" << endl;
    exit(1);
  }
}

void poblacio::destrueix_poblacio()
{
  
  
  forces.close();
  energia.close();
  celules.close();
  
  constantspotencial.close();
  
  direccio_divisio.close();
  flog.close();
  ftime.close();
  
  delete [] matriu_c;
  delete [] matriu_v;
  delete [] matriu_a;
  
  //  std::cout << endl << "***********************S'ha acabat la simulació amb éxit***********************" << std::endl;
  std::cout << endl << "***********************Smile! the simulation is over!***********************" << std::endl;
}

void poblacio::crea_mapa_tipus_celules(int nx, int ny)
{
  int i, j, k;
  int resultat;
  
  for(i=0; i<nx; i++)
  {
    for(j=0; j<ny; j++)
    {
      mapa_tipus_celules[i][j]=-1;
    }
  }
  
  for(i=0; i<nx; i++)
  {
    for(j=0; j<ny; j++)
    {
      for(k=0; k<N_MAXIM_REGIONS; k++)
      {
        resultat=posicio_tipus[k].prova_si_esta_dins(i,j);
        if(resultat==1)
        {
          if(mapa_tipus_celules[i][j]==-1)
          {
            mapa_tipus_celules[i][j]=posicio_tipus[k].tipus;
          }
          else
          {
            std::cout << "Warning: you are probably aware of this (initializing the tissue) but there are regions with cells of type " << mapa_tipus_celules[i][j] << " that overlap with regions with cells of type " << posicio_tipus[k].tipus << std::endl;
            mapa_tipus_celules[i][j]=posicio_tipus[k].tipus;
          }
        }
      }
      if(mapa_tipus_celules[i][j]==-1)
      {
        std::cout << "Error: there are regions of the tissue with cells that do not belong to any type." << std::endl;
        exit(1);
      }
    }
  }
}

void poblacio::crea_celules_inicials(int nx, int ny)
{
  
  int i,j,k;
  int temp_tipus;
  double a,h;
  double x[6],y[6];
  double dx, dy;
  vertex *vertex_temp;
  
  //a=HEXAGONAL_CELL_EDGE;
  <poblacio_cpp_definicioCantoHexagon>
  h=a*sqrt(3.)/2.;
  
  for(i=0; i<nx; i++)
  {
    for(j=0; j<ny; j++)
    {
      dx=(double)(3./2.)*a*((double)i);
      dy=(double)(2.*h*((double)j))+((double)h*((double)(i%2)));
      
      y[0]=h+dy;
      y[1]=0.+dy;
      y[2]=-h+dy;
      y[3]=-h+dy;
      y[4]=0.+dy;
      y[5]=h+dy;
      
      x[0]=a/2.+dx;
      x[1]=a+dx;
      x[2]=a/2.+dx;
      x[3]=-a/2.+dx;
      x[4]=-a+dx;
      x[5]=-a/2.+dx;
      
      
      temp_tipus=mapa_tipus_celules[i][j];
      //<poblacio_cpp_iniciCelula>
      matriu_c[n_matriu_c].inicia_celula(this, n_matriu_c,temp_tipus,area0[temp_tipus][0]);
      
      
      matriu_c[n_matriu_c].proteines.inicia_especies(&(matriu_c[n_matriu_c]));
      matriu_c[n_matriu_c].c_track_id = toString( c_track_idx);
      c_track_idx++;
      
      vertex_temp=&vertex_buit;
      for(k=0; k<7; k++)
      {
        vertex_temp=crea_vertex(&matriu_c[n_matriu_c], vertex_temp, x[k%6], y[k%6]);
      }
      
      matriu_c[n_matriu_c].calcula_perimetre();
      matriu_c[n_matriu_c].calcula_area();
      matriu_c[n_matriu_c].calcula_kappa_area_area0();
      matriu_c[n_matriu_c].calcula_centre();
      matriu_c[n_matriu_c].troba_celules_veines();
      
      for(k=0; k<n_matriu_v; k++)
      {
        matriu_v[k].troba_vertexs_veins();
      }
      n_matriu_c++;   
    }
  }
  
  for(i=0; i<n_matriu_a; i++){
    matriu_a[i].calcula_longitud();
    matriu_a[i].troba_celules_puntes_aresta();
    matriu_a[i].calcula_lambda();
  }
  for(i=0; i<n_matriu_c; i++){
    matriu_c[i].calcula_area();
    matriu_c[i].calcula_perimetre();
    matriu_c[i].calcula_kappa_area_area0();
    matriu_c[i].calcula_propietats();
  }
  for(k=0; k<n_matriu_v; k++)
  {
    matriu_v[k].calcula_forca();
  }
  
}

void poblacio::crea_celula_buida_vertex_buit_i_aresta_buida()
{
  int i;
  
  celula_buida.p=this;
  celula_buida.id=-1;
  celula_buida.ctype=0;
  celula_buida.c_track_id="empty";
  
  for(i=0; i<MAXIM_FASES; i++)
  {
    celula_buida.duracio_cicle[i]=0.;
  }
  
  celula_buida.posicio_cicle=0.;
  celula_buida.pas_cicle=0.;
  celula_buida.area_growth=0.;
  
  celula_buida.area=0.;
  celula_buida.area0=0.;
  celula_buida.perimeter=0.;
  
  celula_buida.calcula_kappa_area_area0();
  
  celula_buida.x=0.;
  celula_buida.y=0.;
  
  celula_buida.ncellvertexes=0;
  
  celula_buida.narestes=0;
  
  celula_buida.neighboringcells=0;
  
  
  vertex_buit.p=this;
  vertex_buit.id=-1;                      //Identificador del vertex
  
  vertex_buit.x=0.;                       //Posicio x i y del vertex
  vertex_buit.y=0.;
  
  vertex_buit.ncelules=0;                 //Numero de celules adjacents al vertex
  vertex_buit.narestes=0;                 //Numero d'arestes adjacents al vertex
  
  for(i=0; i<MAXIM_CELULES_V; i++)        
  {
    vertex_buit.base[i]=0.;                 //Longitud de la base del triangle definit pels dos vertexs
    vertex_buit.vc1[i]=&vertex_buit;        //Punters als vertexs veins a aquest, en ordre horari
    vertex_buit.vc2[i]=&vertex_buit;        //I que pertanyen a la célula "i"
  }
  
  vertex_buit.forca_x=0.;
  vertex_buit.forca_y=0.;
  
  aresta_buida.p=this;    
  aresta_buida.id=-1;
  
  aresta_buida.l=0.;
  
  aresta_buida.v[0]=&vertex_buit;
  aresta_buida.v[1]=&vertex_buit;
  
  aresta_buida.c[0]=&celula_buida;
  aresta_buida.c[1]=&celula_buida;
  
  aresta_buida.cp[0]=&celula_buida;
  aresta_buida.cp[1]=&celula_buida;
  
}

vertex *poblacio::crea_vertex(celula *c, vertex *vanterior, double cx, double cy)
{
  
  int i;
  int flag;
  double temp;
  vertex *vtrobat;
  celula *cveina;
  int aid;
  
  //Iniciem algunes variables per evitar Warnings i possibles bugs.
  
  vtrobat=NULL;
  cveina=NULL;
  aid=-1;
  
  //Comença la funcio.
  
  flag=0;
  for(i=0; i<n_matriu_v; i++)
  {
    temp=sqrt((matriu_v[i].x-cx)*(matriu_v[i].x-cx)+(matriu_v[i].y-cy)*(matriu_v[i].y-cy));
    if(temp<LATICE_BUILD_E)
    {
      flag=1;
      vtrobat=&matriu_v[i];
    }
  }
  if(flag==0)
  {
    vtrobat=&matriu_v[n_matriu_v];
    vtrobat->inicia_vertex(this,n_matriu_v,cx,cy);
    n_matriu_v++;
    
    c->introdueix_vertex(vtrobat, &vertex_buit, &vertex_buit);
  }
  else
  {
    c->introdueix_vertex(vtrobat, &vertex_buit, &vertex_buit);
  }
  
  if(vanterior!=&vertex_buit)
  {
    flag=0;
    for(i=0; i<n_matriu_a; i++)
    {
      if((matriu_a[i].v[0]==vanterior)&&(matriu_a[i].v[1]==vtrobat))
      {
        flag=1;
        aid=i;
      }
      if((matriu_a[i].v[1]==vanterior)&&(matriu_a[i].v[0]==vtrobat))
      {
        flag=1;
        aid=i;
      }
    }
    if(flag==0)
    {
      aid=n_matriu_a;
      matriu_a[aid].id=aid;
      matriu_a[aid].crea_aresta(this, vanterior, vtrobat, c, &celula_buida); //Estic aqui!!!
      n_matriu_a++;
    }
    else
    {
      if(matriu_a[aid].c[0]==&celula_buida)
      {
        cveina=matriu_a[aid].c[1];
      }
      else if(matriu_a[aid].c[1]==&celula_buida)
      {
        cveina=matriu_a[aid].c[0];
      }
      else
      {
        std::cout << endl << "*******************************" << endl;
        std::cout << endl << "Error: something went wrong when computing the neighbors at the edges." << endl;
        std::cout << endl << "*******************************" << endl;
        exit(1);
      }
      matriu_a[aid].crea_aresta(this, vanterior, vtrobat, c, cveina);
    }
  }
  
  return vtrobat;
  
}

void poblacio::fes_vertexs_random()
{
  int i;
  long llavor;
  double temp1, temp2, temp3;
  
  llavor=500;
  
  temp1=HEXAGONAL_CELL_EDGE/3.;
  temp2=temp1;
  
  for(i=0; i<n_matriu_v; i++)
  {
    temp3=box_muller(HEXAGONAL_CELL_EDGE/4., temp2, &llavor);
    if(temp3>temp1)
    {
      temp3=temp1;
    }
    if(temp3<-temp1)
    {
      temp3=-temp1;
    }
    matriu_v[i].x=matriu_v[i].x+temp3;
    
    temp3=box_muller(HEXAGONAL_CELL_EDGE/4., temp2, &llavor);
    if(temp3>temp1)
    {
      temp3=temp1;
    }
    if(temp3<-temp1)
    {
      temp3=-temp1;
    }
    matriu_v[i].y=matriu_v[i].y+temp3;
  }
  
  for(i=0; i<n_matriu_v; i++){
    matriu_v[i].calcula_totes_les_bases();
  }
  for(i=0; i<n_matriu_a; i++){
    matriu_a[i].calcula_longitud();
  }
  for(i=0; i<n_matriu_c; i++){
    matriu_c[i].calcula_area();
    matriu_c[i].calcula_perimetre();
    matriu_c[i].calcula_kappa_area_area0();
  }
  for(i=0; i<n_matriu_v; i++)
  {
    matriu_v[i].calcula_forca();
  }
}

void poblacio::fes_fase_random()
{
  int i;
  double temp;
  long llavor=4768;
  
  for(i=0; i<n_matriu_c; i++){
    temp=ran3(&llavor)*MIG_CICLE*2.;
    temp=temp-MIG_CICLE;
    matriu_c[i].posicio_cicle=temp*proporcio_cicle[matriu_c[i].ctype][0];
  }
  
  
}

void poblacio::defineix_lambda(double valorh, double valorv, int tipus1, int tipus2)
{
  k_lambda_h[tipus1][tipus2]=valorh;
  k_lambda_h[tipus2][tipus1]=valorh;
  k_lambda_v[tipus1][tipus2]=valorv;
  k_lambda_v[tipus2][tipus1]=valorv;
}

void poblacio::destrueix_elements()
{
  serie elimina;
  int i,j;
  
  //Elimina les celules
  
  for(i=0; i<n_celules_destruir; i++)
  {
    elimina.desordenat[i]=celules_per_destruir[i];
  }
  elimina.n_elements=n_celules_destruir;
  elimina.ordena();
  
  for(i=0; i<n_celules_destruir; i++)
  {
    j=elimina.ordenat_de_major_a_menor[i];
    matriu_c[j]=matriu_c[n_matriu_c-(i+1)];
    matriu_c[j].id=j;
    matriu_c[j].cambia_referencia_celula(&matriu_c[n_matriu_c-(i+1)],&matriu_c[j]);
  }
  
  n_matriu_c=n_matriu_c-n_celules_destruir;
  n_celules_destruir=0;
  
  //Elimina els vertexs
  
  for(i=0; i<n_vertexs_destruir; i++)
  {
    elimina.desordenat[i]=vertexs_per_destruir[i];
  }
  elimina.n_elements=n_vertexs_destruir;
  elimina.ordena();
  
  for(i=0; i<n_vertexs_destruir; i++)
  {
    j=elimina.ordenat_de_major_a_menor[i];
    matriu_v[j]=matriu_v[n_matriu_v-(i+1)];
    matriu_v[j].id=j;
    matriu_v[j].cambia_referencia_vertex(&matriu_v[n_matriu_v-(i+1)],&matriu_v[j]);
  }
  
  n_matriu_v=n_matriu_v-n_vertexs_destruir;
  n_vertexs_destruir=0;
  
  //Elimina les arestes
  
  for(i=0; i<n_arestes_destruir; i++)
  {
    elimina.desordenat[i]=arestes_per_destruir[i];
  }
  elimina.n_elements=n_arestes_destruir;
  elimina.ordena();
  
  for(i=0; i<n_arestes_destruir; i++)
  {
    j=elimina.ordenat_de_major_a_menor[i];
    matriu_a[j]=matriu_a[n_matriu_a-(i+1)];
    matriu_a[j].id=j;
    matriu_a[j].cambia_referencia_aresta(&matriu_a[n_matriu_a-(i+1)],&matriu_a[j]);
  }
  
  n_matriu_a=n_matriu_a-n_arestes_destruir;
  n_arestes_destruir=0;
}

<poblacio_cpp_catualitza_constants>

<poblacio_cpp_escriu_constants>

double poblacio::function_hill(double x, double k, int n)
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

double poblacio::function_hill_inverse(double x, double k, int n)
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

double poblacio::function_hill_f(double x, double k, double n)
{
  double resultat;
  
  resultat=pow((x/k), n);
  
  resultat=1./(1.+resultat);
  
  return resultat;
}

double poblacio::function_hill_f_inverse(double x, double k, double n)
{
  double resultat;
  
  resultat=pow((x/k), n);
  
  resultat=resultat/(1.+resultat);
  
  return resultat;
}

double poblacio::white_noise()
{
  return (box_muller(0., 1., &llavor)/SQRTDELTAT);
}

double poblacio::f_step(double x, double l)
{
  return (x>l?1.:0.);
}

double poblacio::f_step_inverse(double x, double l)
{
  return (x>l?0.:1.);
}

void poblacio::crea_celules_inicials_des_de_arxiu()
{
  
  int i,j,k;
  int numero_celules, ctipus, cnvertexs;
  double x[NVERTEXS_ARXIU_ESTAT_INICIAL],y[NVERTEXS_ARXIU_ESTAT_INICIAL];
  vertex *vertex_temp;
  double nullvar;
  string strnullvar;
  
  ifstream finput;
  
  <poblacio_cpp_fitxer_inicial>
  
  finput >> numero_celules;
  finput >> nullvar;
  
  for(i=0; i<numero_celules; i++)
  {
    finput >> strnullvar; //Read and ignore c_track_id;
    finput >> ctipus;
    finput >> cnvertexs;
    finput >> nullvar; //Read and ignore area;
    
    if(cnvertexs > NVERTEXS_ARXIU_ESTAT_INICIAL)
    {
      std::cout << "Error: the cell " << i << " from the input file has " << cnvertexs << " vertexes that is above the maximum allowed number: " << NVERTEXS_ARXIU_ESTAT_INICIAL << "." << std::endl;
      std::cout << "You can change this limit value by modifying the value of the constant NVERTEXS_ARXIU_ESTAT_INICIAL (poblacio.h)." << std::endl;
      exit(1);
    }
    
    matriu_c[n_matriu_c].inicia_celula(this, n_matriu_c,ctipus,area0[ctipus][0]);
    
    matriu_c[n_matriu_c].proteines.llegeix_especies(finput);
    matriu_c[n_matriu_c].c_track_id = toString( c_track_idx);
    c_track_idx++;
    
    finput >> nullvar >> nullvar; //Read and ignore center.
    for(j=0; j<cnvertexs; j++) //Read and ignore neighbor ids.
                {
                  finput >> nullvar;
                }
                
                for(j=0; j<cnvertexs; j++)
                {
                  finput >> x[j] >> y[j];
                }
                
                vertex_temp=&vertex_buit;
    for(k=0; k<(cnvertexs + 1); k++)
    {
      vertex_temp=crea_vertex(&matriu_c[n_matriu_c], vertex_temp, x[k%cnvertexs], y[k%cnvertexs]);
    }
    
    matriu_c[n_matriu_c].calcula_perimetre();
    matriu_c[n_matriu_c].calcula_area();
    matriu_c[n_matriu_c].calcula_kappa_area_area0();
    matriu_c[n_matriu_c].calcula_centre();
    matriu_c[n_matriu_c].troba_celules_veines();
    
    for(k=0; k<n_matriu_v; k++)
    {
      matriu_v[k].troba_vertexs_veins();
    }
    n_matriu_c++;   
  }
  
  for(i=0; i<n_matriu_a; i++){
    matriu_a[i].calcula_longitud();
    matriu_a[i].troba_celules_puntes_aresta();
    matriu_a[i].calcula_lambda();
  }
  for(i=0; i<n_matriu_c; i++){
    matriu_c[i].calcula_area();
    matriu_c[i].calcula_perimetre();
    matriu_c[i].calcula_kappa_area_area0();
    matriu_c[i].calcula_propietats();
  }
  for(k=0; k<n_matriu_v; k++)
  {
    matriu_v[k].calcula_forca();
  }
  
}
