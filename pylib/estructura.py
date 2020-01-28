#!/bin/env/python
# -*- coding: utf-8 -*-

def escriu(f, variable, indentation):
  if isinstance(variable, dict):
    for key in variable:
      f.write(indentation + key + ":\n")
      escriu(f, variable[key], indentation + "  ")
  elif isinstance(variable, list):
    f.write(indentation + "[\n")
    for elem in variable:
      escriu(f, elem, indentation + "  ")
      f.write(indentation + ",\n")
    f.write(indentation + "]\n")
  else:
    f.write(indentation + str(variable) + "\n")
     
def escriuEstructura(stages, proteines, defdglobal):
  f = open("estructura.dat","w")

  f.write("**********************************************\n")
  f.write("Global\n")
  f.write("**********************************************\n")
  escriu(f, defdglobal, "")
  f.write("**********************************************\n")
  f.write("Stages\n")
  f.write("**********************************************\n")
  escriu(f, stages, "")
  f.write("**********************************************\n")
  f.write("Proteines\n")
  f.write("**********************************************\n")
  escriu(f, proteines, "")
  
  f.close()