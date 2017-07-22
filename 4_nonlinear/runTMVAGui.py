#!/usr/bin/env python
import sys
from ROOT import *
argv = sys.argv[1:]
if len(argv) > 0:
    rootfile = argv[0]
else:
    rootfile = "TMVA.root"
gSystem.Load("libTMVAGui.so")

print "\n==> reading Root file", rootfile
TMVA.TMVAGui(rootfile)
gApplication.Run()



