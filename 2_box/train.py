#!/usr/bin/env python
# ----------------------------------------------------------------------------
#  File:        train.py
#  Description: Example of Random Grid Search to find the results of an
#               ensemble cuts. 
#  Created:     10-Jan-2015 Harrison B. Prosper and Sezen Sekmen
#  Updated:     01-Jun-2016 HBP adapt for Bari Lectures and HATS@LPC 2016
# ----------------------------------------------------------------------------
import os, sys, re
from rgsutil import *
from string import *
from ROOT import *
# ----------------------------------------------------------------------------
def main():
    print "="*80
    print "\t=== RGS: Box Cuts ==="
    print "="*80

    # ---------------------------------------------------------------------
    # Load the RGS shared library and check that the various input files
    # exist.
    # ---------------------------------------------------------------------
    if gSystem.Load("libRGS") < 0: error("unable to load libRGS")

    # Name of file containing cut definitions
    # Format of file:
    #   variable-name  cut-type (>, <, <>, |>, |<, ==)
    varfilename = "rgs.cuts"
    if not os.path.exists(varfilename):
        error("unable to open variables file %s" % varfilename)

        
    # Give treename and number of events to use
    treename = "HZZ4LeptonsAnalysisReduced"  # name of Root tree
    start    = 0     # start row 
    numrows  = 10000  # number of events to use
    L        = 300.0 # 1/fb
    maxcuts     = 10000        # maximum number of cut-points to consider
    weightname  = "f_weight"   # name of event weight variable
    sigfilename = '../data/ntuple_4mu_VV.root'
    bkgfilename = '../data/ntuple_4mu_gg.root'
    selection   = "f_massjj>0"
    
    # Check that files exist
    if not os.path.exists(sigfilename):
        error("unable to open signal file %s" % sigfilename)
    if not os.path.exists(bkgfilename):
        error("unable to open background file %s" % bkgfilename)

    # compute weight for files 
    nsig = getEntries(sigfilename, treename)
    wsig = L * float(nsig)/numrows/2.8
    print '=> signal count/weight:     %8d / %8.2f' % (nsig, wsig)
        
    nbkg = getEntries(bkgfilename, treename)
    wbkg = L * float(nbkg)/numrows/2.8    
    print '=> background count/weight: %8d / %8.2f' % (nbkg, wbkg)

    # ---------------------------------------------------------------------
    #  Create RGS object
    #  
    #   The file (cutdatafilename) of cut-points is usually a signal file,
    #   which ideally differs from the signal file on which the RGS
    #   algorithm is run.
    # ---------------------------------------------------------------------
    cutdatafilename = sigfilename
    rgs = RGS(cutdatafilename, start, maxcuts, treename, weightname,
              selection)

    # ---------------------------------------------------------------------
    #  Add signal and background data to RGS object.
    #  Weight each event using the value in the field weightname, if
    #  present.
    #  NB: We asssume all files are of the same format.
    # ---------------------------------------------------------------------
    # The last (optional) argument is a string, which, if given, will be
    # appended to the "count" and "fraction" variables. The "count" variable
    # contains the number of events that pass per cut-point, while "fraction"
    # is count / total, where total is the total number of events per file.
    # If no string is given, the default is to append an integer to the
    # "count" and "fraction" variables, starting at 0, in the order in which
    # the files are added to the RGS object.
    rgs.add(bkgfilename, start, numrows, "_b", wbkg)
    rgs.add(sigfilename, start, numrows, "_s", wsig)

    # ---------------------------------------------------------------------	
    #  Run RGS and write out result
    # ---------------------------------------------------------------------
    rgs.run(varfilename)

    # Write to a root file
    rgsfilename = "%s.root" % nameonly(varfilename)
    rgs.save(rgsfilename)
# ----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print "\tciao!\n"



