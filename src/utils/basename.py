#! /usr/bin/env python
from __future__ import print_function

import sys

if __name__ == "__main__":
  if len(sys.argv) != 2:
    print('wrong number of arguments')
    sys.exit(-1)

  pname = sys.argv[1]

  if pname == "/":
    print("/")
  elif pname == "." or pname == "..":
    print(pname)
  elif "/" not in pname:
    print(pname)
  elif pname[-1] == "/":
    loc = pname[0:-1].rfind("/")
    print(pname[loc+1:-1])
  else:
    loc = pname[0:-2].rfind("/")
    print(pname[loc+1:])

  sys.exit(0)
  
  
