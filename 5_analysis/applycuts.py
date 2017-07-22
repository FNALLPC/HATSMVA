#!/usr/bin/env python
#------------------------------------------------------------------
# File: applycuts.py
# Description: Apply some cuts and get some counts
# Created: 06-Jun-2016 HBP HATS@LPC 2016
#------------------------------------------------------------------
import os, sys
from histutil import *
from time import sleep
from array import array
from ROOT import *
#------------------------------------------------------------------
def readAndFill(filename, treename, which, c, h, pad):
    useBDT = which == 'BDT'
    print "==> reading %s" % filename
    # open ntuple (see histutil.py for implementation)
    ntuple = Ntuple(filename, treename)
    total  = 0
    weight = 0.0
    passweight = 0.0
    for event in ntuple :
        D = event.D_VVgg_BDT if useBDT else event.D_VVgg_MLP
        h.Fill(D, event.D_bkg, event.weight)

        weight += event.weight
        total  += 1
        # --------------------
        # PLACE YOUR CUTS HERE
        # --------------------
        if not (event.D_bkg > 0.5): continue
        if not (event.D_VVgg_MLP > 0.5): continue
        
        passweight += event.weight
        
        if total % 5000 == 0:
            c.cd(pad)
            h.Draw('lego2')
            c.Update()

    print "==> total (unweighted):            %5d" % total        
    print "==>       (weighted):              %8.2f" % weight
    print "==>       (weighted with cuts):    %8.2f\n" % passweight
    c.cd(pad)
    h.Draw('lego2')
    c.Update()
#------------------------------------------------------------------
def main():
    print 
    print "="*80

    # pick discriminant
    if len(sys.argv) > 1:
        which = sys.argv[1]
    else:
        which = 'MLP'
    isBDT = which == 'BDT'
    
    # pick which variable to use from discriminant ntuple
    if len(sys.argv) > 2:
        whichvar = sys.argv[2]
    else:
        whichvar = 'MLP'        

    # set up a standard graphics style	
    setStyle()

    xbins =  20 
    xmin  =   0.0
    xmax  =   1.0

    ybins =  20
    ymin  =   0.0
    ymax  =   1.0

    msize = 0.15
    
    fieldx = 'D_VVgg_%s' % whichvar
    varx = '#font[12]{D_{VV/gg}}'
    
    fieldy = 'D_bkg';
    vary = '#font[12]{D_{H/bkg}}'
    
    treename    = "HZZ4LeptonsAnalysisReduced"
    filenames   = [('d_4mu_simdata.root', 'simdata', kBlack),
                   ('d_4mu_bkg.root',  'ZZ',   kMagenta+1),
                   ('d_4mu_gg.root',   'gg',   kBlue),
                   ('d_4mu_VV.root',   'VV',   kCyan+1)]
    
    # ---------------------------------------------------------
    # make 2-D plots
    # ---------------------------------------------------------
    #gStyle.SetCanvasPreferGL(True)
    c  = TCanvas("fig_VV_gg_bkg_%s_%s" % (which, whichvar),
                 "", 10, 10, 800, 800)
    c.Divide(2, 2)

    # make some plots
    h = []
    s = []
    k = 0
    for i, (filename, title, color) in enumerate(filenames):
        h.append(mkhist2('h%d' % i, varx, vary,
                            xbins, xmin, xmax,
                            ybins, ymin, ymax))
        h[i].Sumw2()
        h[i].SetMarkerSize(msize)
        h[i].GetXaxis().SetTitleOffset(1.9)
        h[i].GetYaxis().SetTitleOffset(1.9)
        h[i].SetMinimum(0)    
        h[i].SetMarkerSize(msize)
        h[i].SetMarkerColor(color)
        
        readAndFill(filename, treename, which, c, h[i], i+1)
        
        j = i % 2
        s.append(Scribe(0.2+j*0.5, 0.95, 0.07))
        if j == 1: k += 1
        c.cd(i+1)
        s[-1].write(title)
    c.SaveAs('.png')
    sleep(5)
#----------------------------------------------------------------------
main()
