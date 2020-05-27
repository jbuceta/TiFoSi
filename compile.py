#!/bin/env/python
# -*- coding: utf-8 -*-

import subprocess
import os
import sys

try:
  subprocess.check_call('rm ./bin/tifosi -f', shell = True)
except:
  pass

try:
  subprocess.check_call('rm ./src/tifosi -f', shell = True)
except:
  pass

try:
  subprocess.call('rm ./src/*.o -f', shell=True)
except:
  pass

os.chdir("./pylib")
try:
  subprocess.check_call(['python', 'config.py'])
except:
  sys.exit("Error during the configuration process...")

print("Compiling source files...")
os.chdir("../src")
try:
  subprocess.check_call(['make'])
except:
  sys.exit("Error during compilation...")

os.chdir("../")
try:
  subprocess.check_call(['mv', './src/tifosi', './bin/'])
except:
  sys.exit("Error moving the binary...")

print("The process has been successfully completed!""")
