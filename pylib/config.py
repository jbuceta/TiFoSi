#!/bin/env/python
# -*- coding: utf-8 -*-

#Error linia 690 aresta.cpp = -> == -> CORREGIT!!

import sys
import xml.etree.ElementTree as etree

from especies import escriuCodiEspecies
from poblacio import escriuCodiPoblacio
from celula import escriuCodiCelula
from estructura import escriuEstructura

def escriuError(llista):
  print "*****************************************"
  print "*       Error parsing config.xml        *"
  print "*****************************************\n"
  for e in llista:
    print e
  print "\n*****************************************"
  sys.exit(1)

def getTag(parent, tag, n):
  result = parent.findall(tag)
  if n != -1:
    if len(result) != n:
      escriuError(["Inside <" + parent.tag + ">, the <" + tag + "> tag is not properly defined or there are a different number of <" + tag + "> tags than expected."])
  if n == 1:
    result = result[0]
  return result

def comprova_propietats_tag(ntag, dictatributs, ftipus):
  for key in dictatributs:
    if key not in ntag.attrib:
      escriuError(["The tag <" + ntag.tag + "> does not have properly defined the attribute '" + key + "'."])
    else:
      if isinstance(dictatributs[key], list): #Comprovem si es una lllista i en aquest cas mirem si el valor de l'atribut en forma part.
        if ntag.attrib[key] not in dictatributs[key]:
          escriuError(["The value of attribute '" + key + "' of tag <" + ntag.tag + "> does not have a valid value. Valid values are:", dictatributs[key]])
      elif dictatributs[key] == None: #Si es None no es fa res
        pass
      else: #Si no es suposa que l'atribut ha de ser del tipus dictatributs[key]
        try:
          prova = dictatributs[key](ntag.attrib[key])
        except:
          escriuError(["The value of attribute '" + key + "' of tag <" + ntag.tag + "> Is not of a valid type. Valid types are:", dictatributs[key]])
  
  if isinstance(ftipus, list): #Comprovem si es una lllista i en aquest cas mirem si el valor de l'atribut en forma part.
    if ntag.text not in ftipus:
      escriuError(["The tag <" + ntag.tag + "> does not have a valid value. Valid values are:", ftipus])
  elif ftipus == None:
    pass
  else:
    try:
      prova = ftipus(ntag.text)
    except:
      escriuError(["The value of tag <" + ntag.tag + "> is of invalid type. Valid types are:",ftipus])
      

def defineixConstant(tconstant, tstages, tistage, nt, ttypes):
  #tconstant->constant, tstages-> stages, tistage->index stage, nt->ntypes constant, ttypes->dglobal['types']
  
  if tconstant.tag not in tstages[tistage]['potential']:
    if nt == 1:
      stages[tistage]['potential'][tconstant.tag] = [0 for i in range(len(ttypes))]
    elif nt == 2:
      stages[tistage]['potential'][tconstant.tag] = [[0 for j in range(len(ttypes))] for i in range(len(ttypes))]
  if nt == 1:
    if 't' in tconstant.attrib:
      t1 = tconstant.attrib['t']
      if t1 not in ttypes:
        escriuError(["The cell type defined on attribute 't' of constant '" + tconstant.tag + "' is not properly defined"])
    else:
      escriuError(["The attribute 't' of constant '" + tconstant.tag + "' is not properly defined"])
    try:
      stages[tistage]['potential'][tconstant.tag][ttypes[t1]] = float(tconstant.text)
    except:
      stages[tistage]['potential']['funcio'] = "si"
      stages[tistage]['potential'][tconstant.tag][ttypes[t1]] = tconstant.text
  elif nt == 2:
    if ('t1' in tconstant.attrib) and ('t2' in tconstant.attrib):
      t1 = tconstant.attrib['t1']
      t2 = tconstant.attrib['t2']
      if (t1 not in ttypes) or (t2 not in ttypes):
        escriuError(["The cell type defined on attributes 't1' or 't2' of constant '" + tconstant.tag + "' is not properly defined"])
    else:
      escriuError(["The attributes 't1' or 't2' of constant '" + tconstant.tag + "' is not properly defined"])
    try:
      stages[tistage]['potential'][tconstant.tag][ttypes[t1]][ttypes[t2]] = float(tconstant.text)
    except:
      stages[tistage]['potential']['funcio'] = "si"
      stages[tistage]['potential'][tconstant.tag][ttypes[t1]][ttypes[t2]] = tconstant.text
    try:
      stages[tistage]['potential'][tconstant.tag][ttypes[t2]][ttypes[t1]] = float(tconstant.text)
    except:
      stages[tistage]['potential']['funcio'] = "si"
      stages[tistage]['potential'][tconstant.tag][ttypes[t2]][ttypes[t1]] = tconstant.text

"""*********************************************************************************
  Main
*********************************************************************************"""

print "Configuring source files..."

tree = etree.parse('../config.xml')
poblacio = tree.getroot()

eGlobal = getTag(poblacio, 'global', 1)
eStages = getTag(poblacio, 'stages', 1)
ePotentials = getTag(poblacio, 'potentials', 1)
eCycles = getTag(poblacio, 'cycles', 1)
eProteins = getTag(poblacio, 'proteins', 1)

"""*********************************************************************************
  Global
*********************************************************************************"""

dglobal = dict()
dglobal['deltat'] = getTag(eGlobal, 'deltat', 1)
comprova_propietats_tag(dglobal['deltat'], dict(), float)
dglobal['deltat'] = dglobal['deltat'].text

dglobal['types'] = dict()
temp = getTag(eGlobal, 'types', 1)
itemp = 1
for t in temp:
  if t.tag not in dglobal['types']:
    if t.tag != "empty":
      dglobal['types'][t.tag] = itemp
      itemp = itemp + 1
    else:
      escriuError(["The 'empty' cell type name is reserved and can not be defined."])
  else:
    escriuError(["There are at least two cell type names repeated."])
dglobal['types']['empty'] = 0

tempitissue = getTag(eGlobal, 'itissue', 1)
temp = getTag(tempitissue, 'file', 1)
dglobal['itissue'] = dict()
try:
  dglobal['itissue']['file'] = str(temp.attrib['f'])
except:
  escriuError(["The tag <file> does not have the attribute 'f' properly defined."])
if len(dglobal['itissue']['file']) == 0:
  temp = getTag(tempitissue, 'ncellsx', 1)
  comprova_propietats_tag(temp, dict(), int)
  dglobal['itissue']['ncelulesx'] = int(temp.text)
  temp = getTag(tempitissue, 'ncellsy', 1)
  comprova_propietats_tag(temp, dict(), int)
  dglobal['itissue']['ncelulesy'] = int(temp.text)
  
  temp = getTag(tempitissue, 'backgroundcells', 1)
  comprova_propietats_tag(temp, dict(), dglobal['types'].keys())
  dglobal['itissue']['cbackground'] = temp.text
  #comprova_propietats_tag(tempx, {'t': [key for key in dglobal['types']]}, int)
  
  regions = getTag(tempitissue, 'square', -1)
  dglobal['itissue']['regions'] = dict()
  dglobal['itissue']['regions']['square'] = list()
  #S'ha de comprovar que el tipus empty no es possible...
  for regio in regions:
    temp = dglobal['types'].keys()
    temp.remove('empty')
    comprova_propietats_tag(regio, {'t': temp, 'ix': int, 'iy': int, 'sizex': int, 'sizey': int}, None)
    dglobal['itissue']['regions']['square'].append({'type': regio.attrib['t'], 'ix': regio.attrib['ix'], 'iy': regio.attrib['iy'], 'sizex': regio.attrib['sizex'], 'sizey': regio.attrib['sizey']})
  
  #Cal agafar les condicions inicials.
"""
      <square t="bulk" ix="bulk" iy="2" sizex="10" sizey="2" />
      <square t="organizer" ix="bulk" iy="2" sizex="10" sizey="2" />
"""
temp = getTag(eGlobal, 'AllowT1C2', -1)
if len(temp) == 0:
  dglobal['AllowT1C2'] = False
else:
  dglobal['AllowT1C2'] = True

temp = getTag(eGlobal, 'StoreEnergyInformation', -1)
if len(temp) == 0:
  dglobal['StoreEnergyInformation'] = False
else:
  dglobal['StoreEnergyInformation'] = True
  
temp = getTag(eGlobal, 'ProteinBinomialDistribution', -1)
if len(temp) == 0:
  dglobal['ProteinBinomialDistribution'] = False
else:
  dglobal['ProteinBinomialDistribution'] = True
  
del temp
del itemp

"""*********************************************************************************
  Stages
*********************************************************************************"""

#funcionsStage = {'proteins':"""fes1
#fes2
#fes3""",
#'tissue':"""tfes1
#tfes2""",
#'grow':"""gfes1
#"""
#}

stages = dict()
for s in eStages:
  if 'order' in s.attrib:
    stages[s.attrib['order']] = dict()
    stages[s.attrib['order']]['funcions'] = [i.tag.strip() for i in s]
    try:
      stages[s.attrib['order']]['duracio'] = int(s.attrib['duration'].strip())
    except:
      escriuError(["Error: The 'duracio' attribute in tag <stage> with 'order'=" + s.attrib['order'] + " is not properly defined. It must be an integer value."])
    try:
      stages[s.attrib['order']]['intermediate'] = int(s.attrib['intermediate'].strip())
    except:
      escriuError(["Error: The 'intermediate' attribute in tag <stage> with 'order'=" + s.attrib['order'] + " is not properly defined. It must be an integer value."])
    stages[s.attrib['order']]['potential'] = dict()
    stages[s.attrib['order']]['cycle'] = dict()
    stages[s.attrib['order']]['proteinConstants'] = dict()
  else:
    escriuError(["Error: The 'order' attribute must be set in all <stage> tags."])

for i in range(len(stages)):
  if str(i+1) not in stages:
    escriuError(["There are " + str(len(stages)) + " stages. However, the stage " + str(i+1) + " is not defined."])


"""*********************************************************************************
  Potential
*********************************************************************************"""
dictConstants = {'KAPPA': 1, 'GAMMA': 1, 'LAMBDA': 2, 'lambda': 2,'force_x': 1, 'force_y': 1}

#S'inicien totes les possibles constants a 0
for ist in stages:
  for ico in dictConstants:
    if dictConstants[ico] == 1:
      stages[ist]['potential'][ico] = [0 for i in range(len(dglobal['types']))]
    elif dictConstants[ico] == 2:
      stages[ist]['potential'][ico] = [[0 for j in range(len(dglobal['types']))] for i in range(len(dglobal['types']))]

temp = getTag(ePotentials, 'potential', -1)
for p in temp:
  if 'stage' in p.attrib:
    stemp = p.attrib['stage'].strip().lower()
  else:
    escriuError(["A <potential> tag have not properly set the 'stage' attribute. The defined attributes are:", p.attrib.keys()])
  if stemp == "all":
    rangetemp = range(1,len(stages)+1)
  else:
    if stemp in stages.keys():
      rangetemp = range(int(stemp),int(stemp)+1)
    else:
      escriuError(["The value of 'stage' attribute of a <potential> tag is not valid. The valid attributes are:", stages.keys().append('all')])

  for s in rangetemp:
    for c in p:
      if c.tag not in dictConstants:
        escriuError(["'" + c.tag + "' is not a valid constant. The valid constants are:", dictConstants.keys()])
      defineixConstant(c, stages, str(s), dictConstants[c.tag], dglobal['types'])

"""*********************************************************************************
  Cell cycle
*********************************************************************************"""

#Iniciem les variables al seu valor per defecte.
for s in range(1,len(stages)+1):
  #dglobal['ConstantSpeed'] = True
  stages[str(s)]['ConstantSpeed'] = True
  for t in dglobal['types']:
    stages[str(s)]['cycle'][t] = dict()
    #stages[str(s)]['cycle'][t]['area0'] = 1.
    #stages[str(s)]['cycle'][t]['vgrow'] = 1.
    stages[str(s)]['cycle'][t]['speed'] = 0.
    stages[str(s)]['cycle'][t]['dispersion'] = 0.8
    stages[str(s)]['cycle'][t]['divisiondispersion'] = 0.2
    stages[str(s)]['cycle'][t]['divisiondispersionlimit'] = 0.5
    stages[str(s)]['cycle'][t]['divisionshift'] = 0.

for c in eCycles:
  if 'stage' in c.attrib:
    stemp = c.attrib['stage'].strip().lower()
  else:
    escriuError(["A <cycle> tag have not properly set the 'stage' attribute. The defined attributes are:", c.attrib.keys()])
  if stemp == "all":
    rangetemp = range(1,len(stages)+1)
  else:
    if stemp in stages.keys():
      rangetemp = range(int(stemp),int(stemp)+1)
    else:
      svalides = stages.keys()
      svalides.append('all')
      escriuError(["The value of 'stage' attribute in a <cycle> tag is " + stemp + ", which is not valid. The valid attributes are:", svalides])
  temp = ['speed','dispersion', 'divisionshift', 'divisiondispersion','divisiondispersionlimit', 'phase']
  for att in c:
    if 't' not in att.attrib:
      escriuError(["The tag <"+att.tag+"> has not defined the attribute 't'."])
    if att.attrib['t'] not in dglobal['types']:
      escriuError(["The attribute 't' of tag <"+att.tag+"> defines a non valid cell type: '"+att.attrib['t']+"'. Valid types are:",dglobal['types'].keys()])
    if att.tag not in temp:
      escriuError(["The tag <"+att.tag+"> is not a valid tag inside <cycle>. Valid tags are:", temp])
    for s in rangetemp:
      if att.tag == 'phase':
        comprova_propietats_tag(att, {'t': None, 'a0i': float, 'a0f': float, 'reldiv': float, 'prop': float}, int)
        if att.tag not in stages[str(s)]['cycle'][att.attrib['t']]:
          stages[str(s)]['cycle'][att.attrib['t']][att.tag] = dict()
        stages[str(s)]['cycle'][att.attrib['t']][att.tag][int(att.text)] = {'a0i': float(att.attrib['a0i']), 'a0f': float(att.attrib['a0f']), 'reldiv': float(att.attrib['reldiv']), 'prop': float(att.attrib['prop'])}
      elif att.tag == 'speed':
        try:
          stages[str(s)]['cycle'][att.attrib['t']][att.tag] = str(float(att.text))
        except:
          stages[str(s)]['cycle'][att.attrib['t']][att.tag] = att.text
          stages[str(s)]['ConstantSpeed'] = False
      else:
        stages[str(s)]['cycle'][att.attrib['t']][att.tag] = att.text

#Comprova que les fases son completes i correctes
for s in range(1,len(stages)+1):
  for t in dglobal['types']:
    if 'phase' not in stages[str(s)]['cycle'][t]: #stages[str(s)]['cycle'][t][att.tag][int(att.text)]
      stages[str(s)]['cycle'][t]['phase'] = dict()
      stages[str(s)]['cycle'][t]['phase'][1] = {'a0i': 1., 'a0f': 1., 'reldiv': 0., 'prop': 0.5}
      stages[str(s)]['cycle'][t]['phase'][2] = {'a0i': 1., 'a0f': 2., 'reldiv': 0.85, 'prop': 0.5}
    else:
      testprop = 0.
      nphases = len(stages[str(s)]['cycle'][t]['phase'])
      if nphases > 20:
        escriuError(["The mx number of phases for a cycle is 20."])
      for i in range(1, nphases+1):
        if i not in stages[str(s)]['cycle'][t]['phase']:
          escriuError(["The <phase> tags of cell type %s in stage %s are not properly defined. The <phase> tag %s is missing." % (t, str(s), str(i+1))])
        testprop = testprop + stages[str(s)]['cycle'][t]['phase'][i]['prop']
      #if testprop > 1.:
      #  escriuError(["The sum of the 'prop' attribute of <phase> tags of cell type %s in stage %s is greather than 1." % (t, str(s))])
        
"""*********************************************************************************
  Protein Dynamics
*********************************************************************************"""
temp = getTag(eProteins, 'constants', -1)
tempConstants = set()
#Llegim les constants de la dinamica de les proteines.
for c in temp:
  if 'stage' in c.attrib:
    stemp = c.attrib['stage'].strip().lower()
  else:
    escriuError(["A <cycle> tag have not properly set the 'stage' attribute. The defined attributes are:", c.attrib.keys()])
  if stemp == "all":
    rangetemp = range(1,len(stages)+1)
  else:
    if stemp in stages.keys():
      rangetemp = range(int(stemp),int(stemp)+1)
    else:
      svalides = stages.keys()
      svalides.append('all')
      escriuError(["The value of 'stage' attribute in a <cycle> tag is " + stemp + ", which not valid. The valid attributes are:", svalides])
              
  for s in rangetemp:
    for k in c:
      tempConstants = tempConstants.union(set([k.tag]))
      try:
        stages[str(s)]['proteinConstants'][k.tag] = float(k.text.strip())
      except:
        escriuError(["The value of the constant <"+k.tag+">: "+k.text.strip()+" is not a valid number."])
#Les constants que no s'hagin definit les igualem a 0.
for s in stages:
  for k in tempConstants:
    if k not in stages[s]['proteinConstants']:
      stages[s]['proteinConstants'][k] = 0.
#Llegim les funcions que definiexien la dinamica de les proteines.
proteines = dict()
temp = getTag(eProteins, 'protein', -1)
for p in temp:
  if "especie" not in p.attrib:
    escriuError(["A protein tag dos not have defined the 'especie' attribute."])
  if p.attrib['especie'] in proteines:
    escriuError(["There are two proteins with the same name:"])
  proteines[p.attrib['especie']] = dict()
  ptemp = getTag(p, 'dynamics', 1)
  proteines[p.attrib['especie']]['edinamica'] = ptemp.text.replace('%c.', 'c->')
  if 'negval' in ptemp.attrib and ptemp.attrib['negval'] == 'y':
    proteines[p.attrib['especie']]['negval'] = True
  else:
    proteines[p.attrib['especie']]['negval'] = False
  proteines[p.attrib['especie']]['inicial'] = dict()
  for ctype in dglobal['types']:
    proteines[p.attrib['especie']]['inicial'][ctype] = dict()
    proteines[p.attrib['especie']]['inicial'][ctype]['valor'] = 0.
    proteines[p.attrib['especie']]['inicial'][ctype]['stochastic'] = 'no'
    proteines[p.attrib['especie']]['inicial'][ctype]['dispersion'] = 0.
  ictemp = getTag(p, 'iconcentration', -1)
  for ic in ictemp:
    if 't' not in ic.attrib:
      escriuError(["A initial condition <"+ic.tag+"> has not defined the attribute 't'."])
    if ic.attrib['t'] not in proteines[p.attrib['especie']]['inicial']:
      escriuError(["A initial condition <"+ic.tag+"> has defined an invalid attribute 't'. Valid attributes are: ",proteines[p.attrib['especie']]['inicial'].keys()])
    try:
      proteines[p.attrib['especie']]['inicial'][ic.attrib['t']]['valor'] = float(ic.text)
    except:
      escriuError(["A initial condition <"+ic.tag+"> has defined an invalid number: "+ic.text])
    if 'stochastic' in ic.attrib:
      if ic.attrib['stochastic'].lower() == 'y' or  ic.attrib['stochastic'].lower() == 'yes':
        proteines[p.attrib['especie']]['inicial'][ic.attrib['t']]['stochastic'] = 'yes'
        if 'dispersion' in ic.attrib:
          try:
            proteines[p.attrib['especie']]['inicial'][ic.attrib['t']]['dispersion'] = float(ic.attrib['dispersion'])
          except:
            escriuError(["The value of the 'dispersion' attribute of <"+ic.tag+"> tag is not a valida number."])
        else:
          escriuError(["A initial condition <"+ic.tag+"> is defined as stochastic, but 'dispersion' attribute is not deffined."])
    
escriuCodiEspecies(proteines, dglobal, stages)
escriuCodiPoblacio(stages, dictConstants, dglobal)
escriuCodiCelula(dglobal, stages)
#escriuEstructura(stages, proteines, dglobal)

print "The source has been configured properly!"
