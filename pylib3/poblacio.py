#!/bin/env/python
# -*- coding: utf-8 -*-

import math
from auxiliar import stageKeysInOrder

#define MAXIM_ELEMENTS_DESTRUIR 10
#define N_TIPUS_CELULA 4
#define N_MAXIM_REGIONS 4

def escriuError(llista):
  print("*****************************************")
  print("*       Error parsing config.xml        *")
  print("*****************************************\n")
  for e in llista:
    print(e)
  print("\n*****************************************")
  sys.exit(1)

def testaConstants(constants, defConstants):
  result = 0
  for c in defConstants:
    if defConstants[c] == 1:
      try:
        test = [float(i) for i in constants[c]]
      except:
        result = result + 1
    else:
      try:
        test = [ [float(j) for j in i] for i in constants[c] ]
      except:
        result = result + 1
  return result
        
def escriuCodiPoblacio(stages, defConstants, defdglobal):
  tipus = defdglobal['types']
  poblacioCPP = dict()
  poblacioH = dict()
  """*********************************************************************************
    Funció dinàmica principal
  *********************************************************************************"""
  temp = "void poblacio::bucle_dinamica_principal()\n{\n"
  temp = temp + "\tint i;\n\n"
  for s in stageKeysInOrder(stages):
    temp = temp + "\tdefineix_constants_stage("+str(s)+");\n"
    temp = temp + "\tfor(i=0; i<PASSOS_STAGE"+str(s)+"; i++)\n\t{\n"
    temp = temp + "\t\tidx_exterior = i;\n"
    temp = temp + "\t\tguarda_dades();\n"
    if testaConstants(stages[s]['potential'], defConstants) != 0:
      temp = temp + "\t\tescriu_constants_stage_%s();\n" % s
    temp = temp + "\t\tbucle_dinamica_parcial_stage_"+str(s)+"();\n"
    temp = temp + "\t}\n"
    temp = temp + "\tguarda_dades_final(\"" + str(s) + "\");\n"
  temp = temp + "}\n\n"

  """*********************************************************************************
    Funcions dinàmica parcial
  *********************************************************************************"""

  for s in stageKeysInOrder(stages):
    temp = temp + "void poblacio::bucle_dinamica_parcial_stage_"+str(s)+"()\n"
    temp = temp + "{\n"
    temp = temp + "\tint i,j;\n"
    temp = temp + "\tfor(j=0; j<PASSOS_BUCLE_PARCIAL%s; j++)\n" % (str(s))
    temp = temp + "\t{\n"
    temp = temp + "\t\tidx_interior = j;\n"
    temp = temp + "\t\ttimepos += DELTAT;\n"
    
    if testaConstants(stages[s]['potential'], defConstants) != 0:
      temp = temp + "\t\tactualitza_constants_stage_%s();\n" % s
    
    if ('growth' in stages[s]['funcions']) or ('mechanics' in stages[s]['funcions']):
      temp = temp + "\t\tfor(i=0; i<n_matriu_v; i++){\n"
      temp = temp + "\t\t\tmatriu_v[i].calcula_forca();\n"
      temp = temp + "\t\t}\n"
      temp = temp + "\t\tfor(i=0; i<n_matriu_v; i++){\n"
      temp = temp + "\t\t\tmatriu_v[i].desplaca_vertex();\n"
      temp = temp + "\t\t}\n"
      temp = temp + "\t\tfor(i=0; i<n_matriu_a; i++){\n"
      temp = temp + "\t\t\tmatriu_a[i].calcula_longitud_dinamica();\n"
      temp = temp + "\t\t}\n"
      temp = temp + "\t\tfor(i=0; i<n_matriu_c; i++){\n"
      if stages[s]['ConstantSpeed'] == False:
        temp = temp + "\t\t\tmatriu_c[i].actualitza_speed_%s();\n" % s
      if 'growth' in stages[s]['funcions']:
        temp = temp + "\t\t\tmatriu_c[i].calcula_rellotge();\n"
      temp = temp + "\t\t\tmatriu_c[i].calcula_area_dinamica();\n"
      temp = temp + "\t\t\tmatriu_c[i].calcula_perimetre();\n"
      temp = temp + "\t\t\tmatriu_c[i].calcula_kappa_area_area0();\n" #Ha de ser una nova funció on la constant depengui de la concentració. També s'ha d'afegir funcions noves...
      if 'proteins' in stages[s]['funcions'] or testaConstants(stages[s]['potential'], defConstants) != 0:
        temp = temp + "\t\t\tmatriu_c[i].calcula_centre();\n"
      temp = temp + "\t\t}\n"
    else:
      if 'proteins' in stages[s]['funcions'] or testaConstants(stages[s]['potential'], defConstants) != 0:
        temp = temp + "\t\tfor(i=0; i<n_matriu_c; i++){\n"
        temp = temp + "\t\t\tmatriu_c[i].calcula_centre();\n"
        temp = temp + "\t\t}\n"
                  
    if 'proteins' in stages[s]['funcions']:
      #temp = temp + "\t\tradi_mitja=0.;\n"
      temp = temp + "\t\tfor(i=0; i<n_matriu_c; i++){\n"
      temp = temp + "\t\t\tmatriu_c[i].proteines.calcula_dinamica_especies();\n"
      #temp = temp + "\t\t\tmatriu_c[i].calcula_radi();\n"
      #temp = temp + "\t\t\tradi_mitja=radi_mitja+matriu_c[i].radi;\n"
      temp = temp + "\t\t}\n"
      #temp = temp + "\t\tradi_mitja=(double)(radi_mitja/((double)n_matriu_c));\n"
      #temp = temp + "\t\tkradi_mitja2=(double)(radi_mitja*radi_mitja)*(3./2.);\n"
      temp = temp + "\t\tfor(i=0; i<n_matriu_c; i++){\n"
      temp = temp + "\t\t\tmatriu_c[i].proteines.actualitza_especies();\n"
      temp = temp + "\t\t}\n"
      
    temp = temp + "\t}\n"
    
    if testaConstants(stages[s]['potential'], defConstants) != 0:
      temp = temp + "\tactualitza_constants_stage_%s();\n" % s
    temp = temp + "\tfor(i=0; i<n_matriu_v; i++){\n"
    temp = temp + "\t\tmatriu_v[i].calcula_energia();\n"
    temp = temp + "\t}\n"
    
    temp = temp + "}\n\n"

  """*********************************************************************************
    Funció canvi stage
  *********************************************************************************"""
  
  relacioConstants = {'KAPPA': ['k_kappa'], 'GAMMA': ['k_gamma'], 'LAMBDA': ['k_lambda_h','k_lambda_v'], 'lambda': ['k_gamma_aresta'],'force_x': ['k_force_x'],'force_y': ['k_force_y']}
  #relacioConstants = {'GAMMA': 1, 'LAMBDA': 2, 'lambda': 2,'force': 1}
  temp = temp + "void poblacio::defineix_constants_stage(int s)\n{\n"
  temp = temp + "\tint i;\n\n"
  for key in tipus:
    temp = temp + "\t//Tipus %s te l'index: %s\n" % (key, str(tipus[key]))
  temp = temp + "\nstage=s;\n"
  for s in stageKeysInOrder(stages):
    temp = temp + "\tif(s=="+str(s)+"){\n"
    """ Definim constants potencial """
    for key in defConstants:
      if defConstants[key] == 1:
        for i in range(len(tipus)):
          for nom_k in relacioConstants[key]:
            tempcodi = str(stages[s]['potential'][key][i])
            if tempcodi.find("%") != -1:
              tempcodi = "0."
              #tempcodi = tempcodi.replace('%c1', 'matriu_c[i].proteines')
            if key == "KAPPA":
              temp = temp + "\t\t%s[%s]=0.5*(%s);\n" % (nom_k,str(i),tempcodi)
            else:
              temp = temp + "\t\t%s[%s]=%s;\n" % (nom_k,str(i),tempcodi)
      if defConstants[key] == 2:
        for i in range(len(tipus)):
          for j in range(len(tipus)):
            for nom_k in relacioConstants[key]:
              tempcodi = str(stages[s]['potential'][key][i][j])
              if tempcodi.find("%") != -1:
                tempcodi = "0."
                #tempcodi = tempcodi.replace('%c1', 'matriu_a[i].c[0]->proteines')
                #tempcodi = tempcodi.replace('%c2', 'matriu_a[i].c[1]->proteines')
              temp = temp + "\t\t%s[%s][%s]=%s;\n" % (nom_k,str(i),str(j),tempcodi)
              
    if testaConstants(stages[s]['potential'], defConstants) != 0:
        temp = temp + "\t\tactualitza_constants_stage_%s();\n" % s               #Aquesta sauria de poder treure ja que esta mes avall.
    """ Definim cicle celular """
    for key in tipus:
      if stages[s]['ConstantSpeed'] == False:
        try:
          temp = temp + "\t\tpas_cicle[%s]=PAS_CICLE*%s;\n" % (str(tipus[key]), str(float(stages[s]['cycle'][key]['speed'])))
        except:
          temp = temp + "\t\tpas_cicle[%s]=PAS_CICLE;\n" % (str(tipus[key]), )
      else:
        temp = temp + "\t\tpas_cicle[%s]=PAS_CICLE*%s;\n" % (str(tipus[key]), str(float(stages[s]['cycle'][key]['speed'])))
      #temp = temp + "\t\tcreixement_area[%s]=%s*INCREMENT_AREA*%s;\n" % (str(tipus[key]), str(stages[s]['cycle'][key]['area0']), stages[s]['cycle'][key]['vgrow'])
      temp = temp + "\t\tdispersio_cicle[%s]=%s;\n" % (str(tipus[key]), stages[s]['cycle'][key]['dispersion'])
      #temp = temp + "\t\tarea0[%s]=%s;\n" % (str(tipus[key]), str(stages[s]['cycle'][key]['area0']))
      temp = temp + "\t\tdivisiondispersion[%s]=%s;\n" % (str(tipus[key]), str(stages[s]['cycle'][key]['divisiondispersion']))
      temp = temp + "\t\tdivisiondispersionlimit[%s]=%s;\n" % (str(tipus[key]), str(stages[s]['cycle'][key]['divisiondispersionlimit']))
      #temp = temp + "\t\tdivisionshift[%s]=%s;\n" % (str(tipus[key]), str(stages[s]['cycle'][key]['divisionshift']))
      nphases = len(stages[s]['cycle'][key]['phase'])
      temp = temp + '\t\tnfases[%s]=%s;\n' % (str(tipus[key]), str(nphases))
      for idx in range(1, nphases + 1):
        temp = temp + '\t\tarea0[%s][%s]=%s;\n' % (str(tipus[key]), str(idx-1), str(stages[s]['cycle'][key]['phase'][idx]['a0i']))
        temp = temp + '\t\tarea0_final[%s][%s]=%s;\n' % (str(tipus[key]), str(idx-1), str(stages[s]['cycle'][key]['phase'][idx]['a0f']))
        temp = temp + '\t\treldiv[%s][%s]=%s;\n' % (str(tipus[key]), str(idx-1), str(stages[s]['cycle'][key]['phase'][idx]['reldiv']))
        temp = temp + '\t\tproporcio_cicle[%s][%s]=%s;\n' % (str(tipus[key]), str(idx-1), str(stages[s]['cycle'][key]['phase'][idx]['prop']))
      
    temp = temp + "\t}\n"
    
  if testaConstants(stages[s]['potential'], defConstants) != 0:
    for s in stageKeysInOrder(stages):
      temp = temp + "\tif(s=="+str(s)+"){\n"
      temp = temp + "\t\tactualitza_constants_stage_%s();\n" % s
      temp = temp + "\t}\n"
  #---
  temp = temp + "\tfor(i=0; i<n_matriu_c; i++){\n"
  temp = temp + "\t\tmatriu_c[i].proteines.inicia_constants(s);\n"
  temp = temp + "\t\tmatriu_c[i].recalcula_fase();\n"
  temp = temp + "\t\tmatriu_c[i].actualitza_constants();\n"
  temp = temp + "\t}\n"
  
  temp = temp + "\tfor(i=0; i<n_matriu_v; i++){\n"
  temp = temp + "\t\tmatriu_v[i].calcula_f0();\n"
  temp = temp + "\t}\n"
  
  temp = temp + "\tfor(i=0; i<n_matriu_a; i++){\n"
  temp = temp + "\t\tmatriu_a[i].calcula_lambda();\n"
  temp = temp + "\t}\n"
  
  temp = temp + "\tfor(i=0; i<n_matriu_v; i++){\n"
  temp = temp + "\t\tmatriu_v[i].calcula_forca();\n"
  temp = temp + "\t\tmatriu_v[i].calcula_energia();\n"
  temp = temp + "\t}\n"

  temp = temp + "}\n\n"
  
  poblacioCPP['<poblacio_cpp_1>'] = temp

  """*********************************************************************************
    Funció actualitza constants
  *********************************************************************************"""
  
  relacioConstantsDepenenPotencial = {'KAPPA': ['k_kappa'], 'GAMMA': ['k_gamma'], 'LAMBDA': ['lambdah','lambdav'], 'lambda': ['k_gamma_aresta'],'force_x': ['k_force_x'],'force_y': ['k_force_y']}
  #relacioConstants = {'GAMMA': 1, 'LAMBDA': 2, 'lambda': 2,'force': 1}

  temp = ''

  for s in stageKeysInOrder(stages):
    llistatemp1 = list()
    llistatemp2 = list()
  
    temp = temp + "void poblacio::actualitza_constants_stage_%s()\n{\n" % s
    temp = temp + "\tint i;\n\n"
    #temp = temp + "\tconstantspotencial << \"STAGE %s STEP \" << idx_interior << std::endl;\n\n" % s

    """ Actualitzem constants potencial """
    for key in defConstants:
      if defConstants[key] == 1:
        for i in range(len(tipus)):
          for nom_k in relacioConstantsDepenenPotencial[key]:
            try:
              esfloat = float(stages[s]['potential'][key][i])
            except:
              temp1 = "\t\tif(matriu_c[i].ctype==%s)\n\t\t{\n" % str(i)
              tempcodi = str(stages[s]['potential'][key][i])
              tempcodi = tempcodi.replace('%c1property', 'matriu_c[i]')
              tempcodi = tempcodi.replace('%c1', 'matriu_c[i].proteines')
              tempcodi = tempcodi.replace('%t', 'timepos')
              if key == "KAPPA":
                temp1 = temp1 + "\t\t\tmatriu_c[i].%s=0.5*(%s);\n" % (nom_k, tempcodi)
              else:
                temp1 = temp1 + "\t\t\tmatriu_c[i].%s=%s;\n" % (nom_k, tempcodi)
              #temp1 = temp1 + "\t\t\tconstantspotencial << \"%s \" << matriu_c[i].c_track_id << \" \" << matriu_c[i].%s << std::endl;\n" % (key, nom_k)
              temp1 = temp1 + "\t\t}\n"
              llistatemp1.append(temp1)
      if defConstants[key] == 2:
        for i in range(len(tipus)):
          for j in range(i, len(tipus)):
            try:
              esfloat = float(stages[s]['potential'][key][i][j])
            except:
              temp1 = "\t\tif(((matriu_a[i].c[0]->ctype==%s)&&(matriu_a[i].c[1]->ctype==%s))||((matriu_a[i].c[0]->ctype==%s)&&(matriu_a[i].c[1]->ctype==%s)))\n\t\t{\n" % (str(i), str(j), str(j), str(i))
              tempcodi = str(stages[s]['potential'][key][i][j])
              tempcodi = tempcodi.replace('%c1property.', 'matriu_a[i].c[0]->')
              tempcodi = tempcodi.replace('%c2property.', 'matriu_a[i].c[1]->')
              tempcodi = tempcodi.replace('%c1', 'matriu_a[i].c[0]->proteines')
              tempcodi = tempcodi.replace('%c2', 'matriu_a[i].c[1]->proteines')
              tempcodi = tempcodi.replace('%cp1property.', 'matriu_a[i].cp[0]->')
              tempcodi = tempcodi.replace('%cp2property.', 'matriu_a[i].cp[1]->')
              tempcodi = tempcodi.replace('%cedge', 'matriu_a[i]')
              tempcodi = tempcodi.replace('%t', 'timepos')
              temp1 = temp1 + "\t\t\tmatriu_a[i].%s=%s;\n" % (relacioConstantsDepenenPotencial[key][0], tempcodi)
              #temp1 = temp1 + "\t\t\tconstantspotencial << \"%s \" << matriu_a[i].c[0]->c_track_id << \" \" << matriu_a[i].c[1]->c_track_id << \" \" << matriu_a[i].%s << std::endl;\n" % (key, relacioConstantsDepenenPotencial[key][0])
              if key == 'LAMBDA':
                temp1 = temp1 + "\t\t\tmatriu_a[i].%s=matriu_a[i].%s;\n" % (relacioConstantsDepenenPotencial[key][1],relacioConstantsDepenenPotencial[key][0])
                temp1 = temp1 + "\t\t\tmatriu_a[i].templh2lv=(matriu_a[i].lambdah-2.*matriu_a[i].lambdav);\n"
                temp1 = temp1 + "\t\t\tmatriu_a[i].templv2lh=(matriu_a[i].lambdav-2.*matriu_a[i].lambdah);\n"
              temp1 = temp1 + "\t\t}\n"
              llistatemp2.append(temp1)
    temp = temp + "\tfor(i=0; i<n_matriu_c; i++){\n"
    temp = temp + "\n".join(llistatemp1)
    temp = temp + "\t}\n"
    temp = temp + "\tfor(i=0; i<n_matriu_a; i++){\n"
    temp = temp + "\n".join(llistatemp2)
    temp = temp + "\t}\n"
    temp = temp + "}\n\n"
  
  poblacioCPP['<poblacio_cpp_catualitza_constants>'] = temp
  
  """*********************************************************************************
    Funció escriu constants
  *********************************************************************************"""
  
  relacioConstantsDepenenPotencial = {'KAPPA': ['k_kappa'], 'GAMMA': ['k_gamma'], 'LAMBDA': ['lambdah','lambdav'], 'lambda': ['k_gamma_aresta'],'force_x': ['k_force_x'],'force_y': ['k_force_y']}
  #relacioConstants = {'GAMMA': 1, 'LAMBDA': 2, 'lambda': 2,'force': 1}

  temp = ''

  for s in stageKeysInOrder(stages):
    llistatemp1 = list()
    llistatemp2 = list()
    
    temp = temp + "void poblacio::escriu_constants_stage_%s()\n{\n" % s
    temp = temp + "\tint i;\n\n"
    temp = temp + "\tconstantspotencial << \"%s\" << std::endl;\n\n" % s

    """ Actualitzem constants potencial """
    for key in defConstants:
      if defConstants[key] == 1:
        for i in range(len(tipus)):
          for nom_k in relacioConstantsDepenenPotencial[key]:
            try:
              esfloat = float(stages[s]['potential'][key][i])
            except:
              temp1 = "\t\tif(matriu_c[i].ctype==%s)\n\t\t{\n" % str(i)
              temp1 = temp1 + "\t\t\tconstantspotencial << \"%s \" << matriu_c[i].c_track_id << \" \" << matriu_c[i].%s << std::endl;\n" % (key, nom_k)
              temp1 = temp1 + "\t\t}\n"
              llistatemp1.append(temp1)
      if defConstants[key] == 2:
        for i in range(len(tipus)):
          for j in range(i, len(tipus)):
            try:
              esfloat = float(stages[s]['potential'][key][i][j])
            except:
              temp1 = "\t\tif(((matriu_a[i].c[0]->ctype==%s)&&(matriu_a[i].c[1]->ctype==%s))||((matriu_a[i].c[0]->ctype==%s)&&(matriu_a[i].c[1]->ctype==%s)))\n\t\t{\n" % (str(i), str(j), str(j), str(i))
              temp1 = temp1 + "\t\t\tconstantspotencial << \"%s \" << matriu_a[i].c[0]->c_track_id << \" \" << matriu_a[i].c[1]->c_track_id << \" \" << matriu_a[i].%s << std::endl;\n" % (key, relacioConstantsDepenenPotencial[key][0])
              temp1 = temp1 + "\t\t}\n"
              llistatemp2.append(temp1)
    temp = temp + "\tfor(i=0; i<n_matriu_c; i++){\n"
    temp = temp + "\n".join(llistatemp1)
    temp = temp + "\t}\n"
    temp = temp + "\tfor(i=0; i<n_matriu_a; i++){\n"
    temp = temp + "\n".join(llistatemp2)
    temp = temp + "\t}\n"
    temp = temp + "}\n\n"
  
  poblacioCPP['<poblacio_cpp_escriu_constants>'] = temp

  """*********************************************************************************
    Definim <poblacio_cpp_regions> de poblacio.cpp.template
  *********************************************************************************"""
  
  if defdglobal['itissue']['file'].strip() == '':
    i = 0
    temp = ""
    
    temp = temp + "\tposicio_tipus[%s].tipus=%s;\n" % (str(i), str(defdglobal['types'][defdglobal['itissue']['cbackground']]))
    temp = temp + "\tposicio_tipus[%s].x=%s;\n" % (str(i),str(0))
    temp = temp + "\tposicio_tipus[%s].y=%s;\n" % (str(i),str(0))
    temp = temp + "\tposicio_tipus[%s].amplada=%s;\n" % (str(i),str(defdglobal['itissue']['ncelulesx']))
    temp = temp + "\tposicio_tipus[%s].alcada=%s;\n\n" % (str(i),str(defdglobal['itissue']['ncelulesy']))
    i = i + 1
    
    for regio in defdglobal['itissue']['regions']['square']:
      temp = temp + "\tposicio_tipus[%s].tipus=%s;\n" % (str(i), str(defdglobal['types'][regio['type']]))
      temp = temp + "\tposicio_tipus[%s].x=%s;\n" % (str(i),str(regio['ix']))
      temp = temp + "\tposicio_tipus[%s].y=%s;\n" % (str(i),str(regio['iy']))
      temp = temp + "\tposicio_tipus[%s].amplada=%s;\n" % (str(i),str(regio['sizex']))
      temp = temp + "\tposicio_tipus[%s].alcada=%s;\n\n" % (str(i),str(regio['sizey']))
      i = i + 1
      
    temp = temp + "\n\tcrea_mapa_tipus_celules(%s,%s);\n" % (str(defdglobal['itissue']['ncelulesx']), str(defdglobal['itissue']['ncelulesy']))
          
    poblacioCPP['<poblacio_cpp_regions>'] = temp
  else:
    poblacioCPP['<poblacio_cpp_regions>'] = ''
  
  if defdglobal['itissue']['file'].strip() == '':
    poblacioCPP['<poblacio_cpp_celules_inicials>'] = "\tcrea_celules_inicials(%s, %s);" % (str(defdglobal['itissue']['ncelulesx']), str(defdglobal['itissue']['ncelulesy']))
    poblacioCPP['<poblacio_cpp_fitxer_inicial>'] = "\tfinput.open(\"\");"
  else:
    poblacioCPP['<poblacio_cpp_celules_inicials>'] = "\tcrea_celules_inicials_des_de_arxiu();"
    poblacioCPP['<poblacio_cpp_fitxer_inicial>'] = "\tfinput.open(\"%s\");" % defdglobal['itissue']['file'].strip()

  """*********************************************************************************
    Definim <poblacio_cpp_iniciCelula> de poblacio.cpp.template
  *********************************************************************************"""
  if defdglobal['itissue']['file'].strip() == '':
    #poblacioCPP['<poblacio_cpp_definicioCantoHexagon>'] = "\ta=sqrt(2.*%s/(3.*sqrt(3.)));\n" % (str(stages[s]['cycle'][defdglobal['itissue']['cbackground']]['area0']))
    poblacioCPP['<poblacio_cpp_definicioCantoHexagon>'] = "\ta=sqrt(2.*%s/(3.*sqrt(3.)));\n" % (str(stages[s]['cycle'][defdglobal['itissue']['cbackground']]['phase'][1]['a0i']))
  else:
    poblacioCPP['<poblacio_cpp_definicioCantoHexagon>'] = "\ta=0.;\n"
  #poblacioCPP['<poblacio_cpp_iniciCelula>']  = "\t\t\tmatriu_c[n_matriu_c].inicia_celula(this, n_matriu_c,temp_tipus,%s);\n" % (str(stages[s]['cycle'][defdglobal['itissue']['cbackground']]['area0']))

  """*********************************************************************************
    Definim definicions poblacio.h.template
  *********************************************************************************"""
  
  temp = ""
  for s in stageKeysInOrder(stages):
    temp = temp + "#define PASSOS_STAGE%s %s\n" % (str(s),stages[s]['duracio'])
  for s in stageKeysInOrder(stages):
    temp = temp + "#define PASSOS_BUCLE_PARCIAL%s %s\n" % (str(s),stages[s]['intermediate'])
  poblacioH['<poblacio_h_duracioStages>'] = temp
  
  if defdglobal['AllowT1C2'] == True:
    poblacioH['<poblacio_h_AllowT1C2>'] = '#define ALLOWT1C2 1\n'
  else:
    poblacioH['<poblacio_h_AllowT1C2>'] = ''
  
  poblacioH['<poblacio_h_deltat>'] = '#define DELTAT ' + str(float(defdglobal['deltat'])) + '\n'
  poblacioH['<poblacio_h_sqrtdeltat>'] = '#define SQRTDELTAT ' + str(math.sqrt(float(defdglobal['deltat']))) + '\n'
  
  if defdglobal['itissue']['file'].strip() == '':
    poblacioH['<poblacio_h_maxNXInicial>'] = '#define MAXIM_NX_INICIAL ' + str(defdglobal['itissue']['ncelulesx']) + '\n'
    poblacioH['<poblacio_h_maxNYInicial>'] = '#define MAXIM_NY_INICIAL ' + str(defdglobal['itissue']['ncelulesy']) + '\n'
  else:
    poblacioH['<poblacio_h_maxNXInicial>'] = '#define MAXIM_NX_INICIAL 0\n'
    poblacioH['<poblacio_h_maxNYInicial>'] = '#define MAXIM_NY_INICIAL 0\n'
  
  poblacioH['<poblacio_h_nTipusCelula>'] = '#define N_TIPUS_CELULA ' + str(len(defdglobal['types'])) + '\n'

  temp = 0
  if defdglobal['itissue']['file'].strip() == '':
    for regio in defdglobal['itissue']['regions']:
      temp = temp + len(defdglobal['itissue']['regions'][regio])
  poblacioH['<poblacio_h_nMaximRegions>'] = '#define N_MAXIM_REGIONS ' + str(temp + 1) + '\n'
  
  temp = '\tdouble k_gamma_aresta[N_TIPUS_CELULA][N_TIPUS_CELULA];\n'
  temp = temp + '\tdouble k_force_x[N_TIPUS_CELULA];\n'
  temp = temp + '\tdouble k_force_y[N_TIPUS_CELULA];\n'
  temp = temp + '\tdouble k_kappa[N_TIPUS_CELULA];\n'
  temp = temp + '\tdouble k_gamma[N_TIPUS_CELULA];\n'
  temp = temp + '\tdouble k_lambda_h[N_TIPUS_CELULA][N_TIPUS_CELULA];\n'
  temp = temp + '\tdouble k_lambda_v[N_TIPUS_CELULA][N_TIPUS_CELULA];\n\n'
  temp = temp + '\tdouble pas_cicle[N_TIPUS_CELULA];\n'
  #temp = temp + '\tdouble creixement_area[N_TIPUS_CELULA];\n'
  temp = temp + '\tint nfases[N_TIPUS_CELULA];\n'
  temp = temp + '\tdouble dispersio_cicle[N_TIPUS_CELULA];\n'
  temp = temp + '\tdouble area0[N_TIPUS_CELULA][MAXIM_FASES];\n'
  temp = temp + '\tdouble area0_final[N_TIPUS_CELULA][MAXIM_FASES];\n'
  temp = temp + '\tdouble reldiv[N_TIPUS_CELULA][MAXIM_FASES];\n'
  temp = temp + '\tdouble proporcio_cicle[N_TIPUS_CELULA][MAXIM_FASES];\n'
  temp = temp + '\tdouble divisiondispersion[N_TIPUS_CELULA];\n'
  temp = temp + '\tdouble divisiondispersionlimit[N_TIPUS_CELULA];\n'
  #temp = temp + '\tdouble divisionshift[N_TIPUS_CELULA];\n'
  poblacioH['<poblacio_h_constants>'] = temp
  
  temp = '\tvoid bucle_dinamica_principal();\n'
  for s in stageKeysInOrder(stages):
    temp = temp + '\tvoid bucle_dinamica_parcial_stage_' + str(s) + '();\n'
  for s in stageKeysInOrder(stages):
    temp = temp + "\tvoid actualitza_constants_stage_%s();\n" % s
  for s in stageKeysInOrder(stages):
    temp = temp + "\tvoid escriu_constants_stage_%s();\n" % s
  temp = temp + '\tvoid defineix_constants_stage(int s);\n'
  poblacioH['<poblacio_h_funcions>'] = temp
  
  
  """*********************************************************************************
    Escribim el source
  *********************************************************************************"""

  #Afegir la definicio de la A0, si despres amb el fitxer d'entrada es canvia, ja es canviara.......

  ft = open("../templates/poblacio.template.cpp", "r")
  fs = open("../src/poblacio.cpp", "w")
  for line in ft:
    if line.strip() in poblacioCPP:
      fs.write(poblacioCPP[line.strip()].replace('\t', '  '))
    else:
      fs.write(line)
  ft.close()
  fs.close()
  
  ft = open("../templates/poblacio.template.h", "r")
  fs = open("../src/poblacio.h", "w")
  for line in ft:
    if line.strip() in poblacioH:
      fs.write(poblacioH[line.strip()].replace('\t', '  '))
    else:
      fs.write(line)
  ft.close()
  fs.close()
