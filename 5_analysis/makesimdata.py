#!/usr/bin/env python
# ---------------------------------------------------------------------
#  File:        makesimdata.py
#  Description: make a simulated data set corresponding to 300/fb.
# ---------------------------------------------------------------------
#  Created:     05-Jun-2015 Harrison B. Prosper
# ---------------------------------------------------------------------
import os, sys, re
from maketree import makeTree
from histutil import Ntuple
from random import uniform
from ROOT import *
# ---------------------------------------------------------------------
# Binary search - Slightly modified version of Jessen Havill (Denison U.)
# algorithm.
# http://personal.denison.edu/~havill
#----------------------------------------------------------------------
def binsearch(L, item):
    first = 0
    last = len(L) - 1   
    found = False
    while (first <= last) and not found:
        mid = (first + last) / 2
        if item <= L[mid]:
            last = mid;       
        elif item > L[mid]:
            first = mid + 1
        if first >= last:
            mid = first
            found = True
        if found: return mid
    return -1
# ---------------------------------------------------------------------
def main():
    print "\n\tmakesimdata.py\n"
    
    treename = "HZZ4LeptonsAnalysisReduced"
    # source names   
    srcnames = ['gg', 'VV', 'bkg']
    # ntuple variable names
    varnames = ['D_VVgg_MLP', 'D_VVgg_BDT', 'D_bkg', 'weight']

    # ---------------------------------------------------
    # 1. load data into memory
    # 2. create a bootstrap sample of a given size
    #    by randomly selecting events with a 
    #    probability proportional to event weight
    # 3. write out events to an ntuple
    # ---------------------------------------------------
    records = []
    weight = 0.0
    for name in srcnames:
        filename = 'd_4mu_%s.root' % name
        print 'read %s' % filename
        ntuple = Ntuple(filename, treename)
        for row in ntuple:
            rec = []
            for varname in varnames:
                rec.append(row(varname))
            records.append(rec)
            weight += records[-1][-1] # weight is last column
        
    print "Total weight (300/fb): %8.2f" % weight
    print "\tcompute cdf of weights"
    wcdf = len(records)*[0]
    wcdf[0] = records[0][-1]
    for i in xrange(1,len(records)):
        wcdf[i] = wcdf[i-1] + records[i][-1]
    sumw = wcdf[-1]
    if sumw != weight:
        sys.exit("huh?")

    # randomly select "N" events according to event weight
    N = int(sumw+0.5) # number of events to select
    print "\tselecting %d events" % N
    outrecords = []
    for i in xrange(N):
        w = uniform(0, sumw)
        k = binsearch(wcdf, w)
        if k < 0:
            sys.exit("**error** not found %f *** should not happen!" % w)
        outrecords.append(records[k])
        outrecords[-1][-1] = 1.0

    # write out records to an ntuple
    filename = 'd_4mu_simdata.root'
    makeTree(filename, treename, outrecords)
# ---------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print "\nbye!"


