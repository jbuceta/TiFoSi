#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include <ctime>

#define MAXIM_FASES 20

#include "especies.h"
#include "celula.h"
#include "vertex.h"
#include "aresta.h"
#include "poblacio.h"

#include "aleatoris.h"

/*******************************************************************************************************
 * 
 *       Parametres que ens defineixen la geometria basica del teixit
 *
 ********************************************************************************************************/

#define A0_DEFAULT 1.
#define HEXAGONAL_CELL_EDGE sqrt(2.*A0_DEFAULT/(3.*sqrt(3.)))

/*******************************************************************************************************
 * 
 *       Parametres que son pàrt del potencial
 *
 ********************************************************************************************************/

#define GAMMA_POTENCIAL 0.04



/*******************************************************************************************************
 * 
 *       Parametres que utilitzem per definir les escales del sistema
 *
 ********************************************************************************************************/

#define RELACIO_ENTRE_TEMPS_RELAXACIO_I_CREIXEMENT 100.
#define PROPORCIO_PER_DIVIDIR 0.85
#define TEMPS_RELAXACIO (HEXAGONAL_CELL_EDGE/GAMMA_POTENCIAL)
#define VELOCITAT_INCREMENT_AREA ((A0_DEFAULT)/(TEMPS_RELAXACIO*(RELACIO_ENTRE_TEMPS_RELAXACIO_I_CREIXEMENT/2.)))

/*******************************************************************************************************
 * 
 *       Parametres que ens defineixen com sera el cicle celular
 *
 ********************************************************************************************************/

#define INCREMENT_AREA (DELTAT*VELOCITAT_INCREMENT_AREA)
#define PAS_CICLE 1.





#define MIG_CICLE ((TEMPS_RELAXACIO*(RELACIO_ENTRE_TEMPS_RELAXACIO_I_CREIXEMENT/2.))/DELTAT)

/********************************************************************************************************
 * 
 *       Parametres que indiquen si es realitzen certes parts de la
 *       simulacio o quin tipus de dades es guarden.
 *
 ********************************************************************************************************/

#define GUARDA_VEINS_CELULA True

/********************************************************************************************************
 * 
 *       Parametres que indiquen com es distribuiran els diferents tipus de celules
 *       en el teixit.
 *
 ********************************************************************************************************/
