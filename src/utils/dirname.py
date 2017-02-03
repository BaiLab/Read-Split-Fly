#! /usr/bin/env python
from __future__ import print_function

import sys


if __name__ == "__main__":
  if len(sys.argv) != 2:
    print('dirname: wrong number of arguments')
    sys.exit(-1)

  pname = sys.argv[1]

  if pname == "/":
    print("/")
  elif pname == "." or pname == "..":
    print(".")
  elif "/" not in pname:
    print(".")
  elif pname[-1] == "/":
    loc = pname[0:-1].rfind("/")
    if loc == 0:
      print("/")
    elif loc == -1:
      print(".")
    else:
      print(pname[0:loc])
  else:
    loc = pname.rfind("/")
    if loc == 0:
      print("/")
    else:
      print(pname[0:loc])

  sys.exit(0)
  
  
