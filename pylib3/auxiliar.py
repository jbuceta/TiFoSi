#!/bin/env/python
# -*- coding: utf-8 -*-

import sys

def stageKeysInOrder(stages):
  keys = [str(s) for s in range(1, len(stages) + 1)]
  if set(keys) != set(stages.keys()) or len(set(keys)) != len(keys):
    print("Error: An internal error was happened in stageKeysInOrder() function")
    sys.exit(1)
  return keys

