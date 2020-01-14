#!/bin/env/python
# -*- coding: utf-8 -*-

from auxiliar import stageKeysInOrder

def escriuCodiEspecies(data, dglobal, stages):
  indexs = dglobal['types']
  especiesCPP = dict()
  especiesH = dict()
  """*********************************************************************************
    Funció dinàmica
  *********************************************************************************"""
  temp = "void especies::calcula_dinamica_especies()\n{\n"
  for e in data:
    if data[e]['negval'] == True:
      temp = temp + "\t" + e + "_temp = " + e + " + ("+data[e]['edinamica']+")*DELTAT;\n"
    else:
      temp = temp + "\t" + e + "_temp = " + e + " + ("+data[e]['edinamica']+")*DELTAT;\n"
      temp = temp + "\t" + e + "_temp = ((" + e + "_temp<0.)?0.:" + e + "_temp);\n"
  temp = temp + "}\n\n"
  
  """*********************************************************************************
    Funció actualitza concentracions
  *********************************************************************************"""
  temp = temp + "void especies::actualitza_especies()\n{\n"
  for e in data:
    temp = temp + "\t" + e + " = " + e + "_temp;\n"
  temp = temp + "}\n\n"
  
  """*********************************************************************************
    Funció defineix concentracions inicials
  *********************************************************************************"""
  temp = temp + "void especies::inicia_especies(celula *ctemp)\n{\n"
  temp = temp + "\tc = ctemp;\n"
  
  for ct in indexs:
    if str(ct) != 'empty':
      temp = temp + "\tif(c->ctype == " + str(indexs[ct]) + ")\n\t{\n"
      for e in data:
        if data[e]['inicial'][ct]['stochastic'] == 'yes':
          temp = temp + "\t\t" + e + " = box_muller(" + str(data[e]['inicial'][ct]['valor']) + ","+str(data[e]['inicial'][ct]['dispersion'])+", &(c->p->llavor));\n"
          temp = temp + "\t\tif(%s < 0.)\n\t\t{\n" % e
          temp = temp + "\t\t\t%s = 0.;\n" % e
          temp = temp + "\t\t}\n"
        else:
          temp = temp + "\t\t" + e + " = " + str(data[e]['inicial'][ct]['valor']) + ";\n"
      temp = temp + "\t}\n"
  temp = temp + "}\n\n"

  """*********************************************************************************
    Funció divideix especies entre celules
  *********************************************************************************"""

  if dglobal['ProteinBinomialDistribution'] == True:
    temp = temp + "void especies::divideix_celula(celula *cfilla)\n{\n"

    for e in data:
      temp = temp + "\tcfilla->proteines." + e + " = " + e + ";\n"
      temp = temp + "\t" + e + " = box_muller(" + e + "/2., sqrt(" + e + "*0.5*0.5), &(c->p->llavor));\n"
      temp = temp + "\tif(" + e + " < 0.)\n\t{\n"
      temp = temp + "\t\t" + e + " = 0.;\n"
      temp = temp + "\t}\n"
      temp = temp + "\tif(" + e + " > cfilla->proteines." + e + ")\n\t{\n"
      temp = temp + "\t\t" + e + " = cfilla->proteines." + e + ";\n"
      temp = temp + "\t\tcfilla->proteines." + e + " = 0.;\n"
      temp = temp + "\t}\n"
      temp = temp + "\telse\n\t{\n"
      temp = temp + "\t\tcfilla->proteines." + e + " = cfilla->proteines." + e + " - " + e + ";\n"
      temp = temp + "\t}\n"
      temp = temp + "\t" + e + "_temp = " + e + ";\n"
      temp = temp + "\tcfilla->proteines." + e + "_temp = cfilla->proteines." + e + ";\n\n"
    temp = temp + "}\n\n"
  else:
    temp = temp + "void especies::divideix_celula(celula *cfilla)\n{\n"

    for e in data:
      temp = temp + "\t" + e + " = " + e + "/2.;\n"
      temp = temp + "\t" + e + "_temp = " + e + "_temp/2.;\n"
      temp = temp + "\tcfilla->proteines." + e + " = " + e + ";\n"
      temp = temp + "\tcfilla->proteines." + e + "_temp = " + e + "_temp;\n"
    temp = temp + "}\n\n"

  """*********************************************************************************
    Funció escriu espècies
  *********************************************************************************"""

  especiesEnOrdre = ""

  temp = temp + "void especies::escriu_especies(ofstream &arxiu)\n{\n"

  temp = temp + "\tarxiu << " + str(len(data)) + " << \" \";\n"
  for e in data:
      temp = temp + "\tarxiu << " + str(e) + " << \" \";\n"
      especiesEnOrdre += str(e) + " "
  temp = temp + "}\n\n"

  """*********************************************************************************
    Funció retorna nom espècies en ordre
  *********************************************************************************"""

  temp = temp + "std::string especies::especies_en_ordre()\n{\n"
  temp = temp + "\treturn \"%s\";\n" % especiesEnOrdre.strip()
  temp = temp + "}\n\n"

  """*********************************************************************************
    Funció llegeix espècies
  *********************************************************************************"""

  temp = temp + "void especies::llegeix_especies(ifstream &arxiu)\n{\n"
  temp = temp + "\tint temp;\n"
  temp = temp + "\tarxiu >> temp;"
  for e in data:
      temp = temp + "\tarxiu >> c->proteines." + str(e) + ";"
  temp = temp + "\n}\n\n"
  
  """*********************************************************************************
    Funció defineix les constants
  *********************************************************************************"""
  temp = temp + "void especies::inicia_constants(int stage)\n{\n"
  for s in stageKeysInOrder(stages):
    temp = temp + "\tif(stage == " + str(s) + ")\n\t{\n"
    for k in stages[s]['proteinConstants']:
      temp = temp + "\t\t" + k + "=" + str(stages[str(s)]['proteinConstants'][k]) + ";\n"
    temp = temp + "\t}\n"
  temp = temp + "}\n\n"
  
  
  especiesCPP['<especies_cpp_1>'] = temp

  """*********************************************************************************
    Funció difusió
  *********************************************************************************"""
  temp = ""
  for e in data:
    temp = temp + "double especies::difusion_%s()\n{\n" % e
      
    temp = temp + "\tdouble  resultat;\n"
    temp = temp + "\tint i;\n"
  
    temp = temp + "\tresultat = 0.;\n"
    temp = temp + "\tfor(i=0; i<c->ncelules; i++)\n"
    temp = temp + "\t{\n"
    temp = temp + "\t\tif(c->c[i]->id!=-1)\n"
    temp = temp + "\t\t{\n"
    temp = temp + "\t\t\tresultat=resultat+((c->c[i]->proteines.%s*c->area/c->c[i]->area)-%s);\n" % (e, e)
    temp = temp + "\t\t}\n"
    temp = temp + "\t}\n"
  
    temp = temp + "\tresultat = resultat/(c->p->kradi_mitja2);\n"

    temp = temp + "\treturn resultat;\n"
    temp = temp + "}\n\n"
    
  especiesCPP['<especies_cpp_difusio_old>'] = temp
  
  """*********************************************************************************
    Funció difusió 2
  *********************************************************************************"""
  temp = ""
  for e in data:
    temp = temp + "double especies::difusion_%s()\n{\n" % e
      
    temp = temp + "\tdouble  resultat;\n"
    temp = temp + "\tint i;\n"
  
    temp = temp + "\tresultat = 0.;\n"
    temp = temp + "\tfor(i=0; i<c->ncelules; i++)\n"
    temp = temp + "\t{\n"
    temp = temp + "\t\tif(c->c[i]->id!=-1)\n"
    temp = temp + "\t\t{\n"
    temp = temp + "\t\t\tresultat=resultat+(c->a[i]->l/calcula_r(c, c->c[i]))*((c->c[i]->proteines.%s/c->c[i]->area)-(%s/c->area));\n" % (e, e)
    temp = temp + "\t\t}\n"
    temp = temp + "\t}\n"
  
    temp = temp + "\treturn resultat;\n"
    temp = temp + "}\n\n"
    
  especiesCPP['<especies_cpp_difusio>'] = temp

  """*********************************************************************************
    Funció senyal
  *********************************************************************************"""

  temp = ""
  for e in data:
    temp = temp + "double especies::signal_%s()\n{\n" % e
    temp = temp + "\tdouble temp, resultat;\n"
    temp = temp + "\tint i;\n"

    temp = temp + "\ttemp=0.;\n"
    temp = temp + "\tfor(i=0; i<c->ncelules; i++)\n"
    temp = temp + "\t{\n"
    temp = temp + "\t\tif(c->c[i]->id!=-1)\n"
    temp = temp + "\t\t{\n"
    temp = temp + "\t\t\ttemp=temp+((c->a[i]->l/c->c[i]->perimetre)*(c->c[i]->proteines.%s));\n" % e
    temp = temp + "\t\t}\n"
    temp = temp + "\t}\n"

    temp = temp + "\tresultat=temp;\n"
  
    temp = temp + "\treturn resultat;\n"
    temp = temp + "}\n\n"
    
  especiesCPP['<especies_cpp_senyal>'] = temp
  
  """*********************************************************************************
    Funció senyal raw
  *********************************************************************************"""

  temp = ""
  for e in data:
    temp = temp + "double especies::signal_raw_%s()\n{\n" % e
    temp = temp + "\tdouble temp, resultat;\n"
    temp = temp + "\tint i;\n"

    temp = temp + "\ttemp=0.;\n"
    temp = temp + "\tfor(i=0; i<c->ncelules; i++)\n"
    temp = temp + "\t{\n"
    temp = temp + "\t\tif(c->c[i]->id!=-1)\n"
    temp = temp + "\t\t{\n"
    temp = temp + "\t\t\ttemp=temp+(c->c[i]->proteines.%s);\n" % e
    temp = temp + "\t\t}\n"
    temp = temp + "\t}\n"

    temp = temp + "\tresultat=temp;\n"
  
    temp = temp + "\treturn resultat;\n"
    temp = temp + "}\n\n"
    
  especiesCPP['<especies_cpp_senyal_raw>'] = temp


  """*********************************************************************************
    Escriu especies.h
  *********************************************************************************"""
  temp = ""
  
  for e in data:
    temp = temp + "\tdouble " + e + ";\n"
    temp = temp + "\tdouble " + e + "_temp;\n"

  temp = temp + "\n"
  for s in stageKeysInOrder(stages):
    for k in stages[s]['proteinConstants']:
      temp = temp + "\tdouble " + str(k) + ";\n"""
    break


  especiesH['<especies_h_constants>'] = temp
  
  temp = ""
  for e in data:
    temp = temp + "\tdouble difusion_%s();\n" % e
  temp = temp + "\n"
  for e in data:
    temp = temp + "\tdouble signal_%s();\n" % e
  temp = temp + "\n"
  for e in data:
    temp = temp + "\tdouble signal_raw_%s();\n" % e
  
  especiesH['<especies_h_funcions>'] = temp

  """*********************************************************************************
    Escribim el source
  *********************************************************************************"""

  ft = open("../templates/especies.template.cpp", "r")
  fs = open("../src/especies.cpp", "w")
  for line in ft:
    if line.strip() in especiesCPP:
      fs.write(especiesCPP[line.strip()].replace('\t', '  '))
    else:
      fs.write(line)
  ft.close()
  fs.close()

  ft = open("../templates/especies.template.h", "r")
  fs = open("../src/especies.h", "w")
  for line in ft:
    if line.strip() in especiesH:
      fs.write(especiesH[line.strip()].replace('\t', '  '))
    else:
      fs.write(line)
  ft.close()
  fs.close()
