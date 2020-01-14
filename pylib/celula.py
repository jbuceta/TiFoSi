#!/bin/env/python
# -*- coding: utf-8 -*-

from auxiliar import stageKeysInOrder

def escriuCodiCelula(dglobal, stages):
  indexs = dglobal['types']
  celulaCPP = dict()
  celulaH = dict()


  temp = ''
  """*********************************************************************************
    Funció actualitza speed del cicle
  *********************************************************************************"""
  for s in stageKeysInOrder(stages):
    temp = temp + "void celula::actualitza_speed_%s()\n{\n" % str(s)
    tflag = 0
    for ct in indexs:
      if str(ct) != 'empty':
        if tflag == 0:
          temp = temp + "\tif(ctype == " + str(indexs[ct]) + ")\n\t{\n"
          tflag = 1
        else:
          temp = temp + "\telse if(ctype == " + str(indexs[ct]) + ")\n\t{\n"
        stringtemp = str(stages[str(s)]['cycle'][ct]['speed'])
        stringtemp = stringtemp.replace('%c1', 'proteines')
        temp = temp + "\t\tpas_cicle = %s;\n" % stringtemp
        temp = temp + "\t}\n"

    temp = temp + "\telse\n\t{\n"
    temp = temp + "\t\tstd::cout << endl << \"*******************************\" << endl;\n"
    temp = temp + "\t\tstd::cout << endl << \"Error. There is a bug in the cycle speed computation.\" << endl;\n"
    temp = temp + "\t\tstd::cout << endl << \"*******************************\" << endl;\n"
    temp = temp + "\t\texit(1);\n"
    temp = temp + "\t}\n"
    
    temp = temp + "}\n\n"
    
  celulaCPP['<celula_cpp_actualitza_speed>'] = temp
  
  """*********************************************************************************
    Escriu especies.h
  *********************************************************************************"""
  temp = ''
  for s in stageKeysInOrder(stages):
    temp = temp + "\t\t\t\tvoid actualitza_speed_%s();\n" % str(s)
  
  celulaH['<celula_h_funcions>'] = temp

    
  """*********************************************************************************
    Escribim el source
  *********************************************************************************"""

  ft = open("../templates/celula.template.cpp", "r")
  fs = open("../src/celula.cpp", "w")
  for line in ft:
    if line.strip() in celulaCPP:
      fs.write(celulaCPP[line.strip()].replace('\t', '  '))
    else:
      fs.write(line)
  ft.close()
  fs.close()

  ft = open("../templates/celula.template.h", "r")
  fs = open("../src/celula.h", "w")
  for line in ft:
    if line.strip() in celulaH:
      fs.write(celulaH[line.strip()].replace('\t', '  '))
    else:
      fs.write(line)
  ft.close()
  fs.close()
