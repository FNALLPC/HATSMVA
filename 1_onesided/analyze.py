#!/usr/bin/env python
# ---------------------------------------------------------------------
#  File:        analyze.py
#  Description: Analyze the results of RGS and find the best cuts.
# ---------------------------------------------------------------------
#  Created:     10-Jan-2015 Harrison B. Prosper and Sezen Sekmen
# ---------------------------------------------------------------------
import os, sys, re
from string import *
from rgsutil import *
from histutil import *
from time import sleep
from array import array
from ROOT import *
# ---------------------------------------------------------------------
START_ROW=5000
CWD=getCWD()
# ---------------------------------------------------------------------
def plotData():
    xbins = 16
    xmin  =  0.0
    xmax  =  8.0

    ybins =  16
    ymin  =   0.0
    ymax  =5000.0    

    msize = 0.15 # marker size

    fieldx = 'f_deltajj'; varx = '|#Delta#font[12]{#eta_{jj}}|'
    fieldy = 'f_massjj';  vary = '#font[12]{m_{jj}}'
    weightname  = "f_weight"
    treename    = "HZZ4LeptonsAnalysisReduced"
    sigfilename = '../data/ntuple_4mu_VV.root'
    bkgfilename = '../data/ntuple_4mu_gg.root'
    
    cmass = TCanvas("fig_%s_VV_gg" % CWD, "VBF/ggF", 10, 10, 500, 500)    

    # -- signal
    hsig = mkhist2("hsig", varx, vary,
                    xbins, xmin, xmax,
                    ybins, ymin, ymax)
    hsig.Sumw2()
    hsig.SetMarkerSize(msize)
    hsig.GetYaxis().SetTitleOffset(2.10)
    hsig.SetMinimum(0)    
    hsig.SetMarkerSize(msize)
    hsig.SetMarkerColor(kCyan+1)
          
    sntuple = Ntuple(sigfilename, treename, START_ROW)
    total   = 0
    for event in sntuple:
        if not (event.f_massjj > 0): continue
        hsig.Fill(event.f_deltajj, event.f_massjj, event.f_weight)
        if total % 5000 == 0:
            cmass.cd()
            hsig.Draw('p')
            cmass.Update()
        total += 1
    
    # -- background
    hbkg = mkhist2("hbkg", varx, vary,
                    xbins, xmin, xmax,
                    ybins, ymin, ymax)
    hbkg.Sumw2()
    hbkg.SetMarkerSize(msize)
    hbkg.GetYaxis().SetTitleOffset(2.10)
    hbkg.SetMinimum(0)    
    hbkg.SetMarkerSize(msize)
    hbkg.SetMarkerColor(kMagenta+1)     

    bntuple = Ntuple(bkgfilename, treename, START_ROW)
    total  = 0
    for event in bntuple:
        if not (event.f_massjj > 0): continue
        hbkg.Fill(event.f_deltajj, event.f_massjj, event.f_weight)
        total += 1
        if total % 5000 == 0:
            cmass.cd()
            hbkg.Draw('p')
            cmass.Update()
        total += 1
        
    hsig.Scale(1.0/hsig.Integral())
    hbkg.Scale(1.0/hbkg.Integral())
    return (cmass, hsig, hbkg)
# ---------------------------------------------------------------------
def main():
    print "="*80
    print "\t=== RGS: One-Sided Cuts ==="
    print "="*80

    resultsfilename = "rgs.root"
    treename = "RGS"
    print "\n\topen RGS file: %s"  % resultsfilename
    ntuple = Ntuple(resultsfilename, treename)
    
    variables = ntuple.variables()
    for name, count in variables:
        print "\t\t%-30s\t%5d" % (name, count)        
    print "\tnumber of cut-points: ", ntuple.size()

    # -------------------------------------------------------------
    # Plot results of RGS, that is, the fraction of events that
    # pass a given cut-point.
    #  1. Loop over cut points and compute a significance measure
    #     for each cut-point.
    #  2. Find cut-point with highest significance.
    # -------------------------------------------------------------
    # Set up a standard Root graphics style (see histutil.py in the
    # python directory).
    setStyle()

    cmass, hs, hb = plotData()
    
    # Create a 2-D histogram for ROC plot
    msize = 0.30  # marker size for points in ROC plot
    
    xbins =   50  # number of bins in x (background)
    xmin  =  0.0  # lower bound of x
    xmax  =  0.25 # upper bound of y

    ybins =   50
    ymin  =  0.0
    ymax  =  1.0

    color = kBlue+1
    hist  = mkhist2("hroc",
                    "#font[12]{f_{B}}",
                    "#font[12]{f_{S}}",
                    xbins, xmin, xmax,
                    ybins, ymin, ymax,
                    color=color)
    hist.SetMinimum(0)
    hist.SetMarkerSize(msize)


    # loop over all cut-points, compute a significance measure Z
    # for each cut-point, and find the cut-point with the highest
    # significance and the associated cuts.
    print "\tfilling ROC plot..."	
    bestZ = -1 # best Z value
    best  = -1 # row number of with best cuts

    for row, cuts in enumerate(ntuple):
        fb = cuts.fraction_b  #  background fraction
        fs = cuts.fraction_s  #  signal fraction
        b  = cuts.count_b     #  background count
        s  = cuts.count_s     #  signal count
                
        #  Plot fs vs fb
        hist.Fill(fb, fs)
        	
        # Compute measure of significance
        #   Z  = sign(LR) * sqrt(2*|LR|)
        # where LR = log(Poisson(s+b|s+b)/Poisson(s+b|b))
        Z = 0.0
        if b > 1:
            Z = 2*((s+b)*log((s+b)/b)-s)
            absZ = abs(Z)
            if absZ != 0:
                Z = Z*sqrt(absZ)/absZ                    
        if Z > bestZ:
            bestZ = Z
            best  = row

    # -------------------------------------------------------------            
    # get best cut
    # -------------------------------------------------------------
    ntuple.read(best)
    print "\nBest cuts"
    print "Z-value = %-8.2f" % bestZ

    varfilename = "rgs.cuts"
    cutdirs  = getCutDirections(varfilename)
    bestcuts = {}
    for name, cutdir in cutdirs:    
        value = ntuple.get(name)
        bestcuts[name] = value
        print "\t%10s %s %-10.2f" % (name, cutdir, value)
    print
    
    print "Yields and relative efficiencies"
    for name in ['count_s', 'fraction_s',
                 'count_b', 'fraction_b']:    
        value = ntuple(name)
        print "\t%-30s %10.3f" % (name, value)
        if name[0:5] == "fract": print
    
    # -------------------------------------------------------------
    # Save plots
    # -------------------------------------------------------------
    print "\t=== plot ROC ==="
    	
    croc = TCanvas("fig_%s_ROC" % CWD, "ROC", 600, 10, 500, 500)
    croc.cd()
    hist.Draw()

    x = array('d'); x.append(ntuple('fraction_b'))
    y = array('d'); y.append(ntuple('fraction_s'))
    g = TGraph(1, x, y)
    g.SetMarkerSize(1.8)
    g.SetMarkerColor(kRed)
    g.Draw('p')
    croc.Update()

    print "\t=== plot cuts ==="

    xmin = hs.GetXaxis().GetBinLowEdge(1)
    xmax = hs.GetXaxis().GetBinUpEdge(xbins)
    ymin = hs.GetYaxis().GetBinLowEdge(1)
    ymax = hs.GetYaxis().GetBinUpEdge(ybins)    
    hcut = TH2Poly('hcut', '', xmin, xmax, ymin, ymax)
    
    xcut = bestcuts['f_deltajj']
    ycut = bestcuts['f_massjj']
    xname, xdir = cutdirs[0]
    yname, ydir = cutdirs[1]
    
    xlow = xmin
    xupp = xmax
    ylow = ymin
    yupp = ymax
    if xdir == '>':
        xlow = xcut
        if ydir == '>':
            ylow = ycut
        else:
            yupp = ycut
    else:
        xupp = xcut
        if ydir == '>':
            ylow = ycut
        else:
            yupp = ycut

    hcut.AddBin(xlow, ylow, xupp, yupp)
    cmass.cd()
    hs.Draw('p')
    hb.Draw('p same')
    hcut.Draw('same')

    # save plots
    croc.SaveAs(".png")    
    cmass.SaveAs('.png')
    
    sleep(5)
# ---------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print "bye!"


