#!/usr/bin/env python
#------------------------------------------------------------------
# Description: plot variables
# HATS@LPC 2016 Harrison B. Prosper
#------------------------------------------------------------------
import os,sys
from string import *
from math import *
from ROOT import *
from histutil import *
from time import sleep
#------------------------------------------------------------------
# potential discriminating variables
VARS = '''
f_deltajj
f_massjj
f_pt4l
f_pfmet
f_Z2mass
f_lept1_pt
'''
VARS = map(strip, split(strip(VARS),'\n'))
#------------------------------------------------------------------
# read data, cache them, and normalize them
def readData(filename, treename, d1=None, d2=None):
    print '\n=> reading file %s' % filename
    ntuple = Ntuple(filename, treename)
    accumulate = d1 == None
    if accumulate:
        d1 = [0.0]*len(VARS)
        d2 = [0.0]*len(VARS)
    data  = []
    total = 0
    weight= 0.0
    for event in ntuple:
        if not (event.f_massjj > 0): continue
        
        if total % 5000 == 0: print total
            
        total  += 1
        weight += event.f_weight
        
        d = [0]*len(VARS)
        for i, var in enumerate(VARS):
            x = eval('event.%s' % var)
            if accumulate: 
                d1[i] += x * event.f_weight
                d2[i] += x*x * event.f_weight
            d[i]   = x
            
        data.append((d, event.f_weight))
    print "unweighted count: %d\tweighted count: %8.2f" % (total, weight)

    for i, var in enumerate(VARS):
        if accumulate:
            d1[i] /= weight
            d2[i] /= weight
            d2[i]  = sqrt(d2[i]-d1[i]*d1[i])
        
    for d, w in data:
        for i, x in enumerate(d):
            d[i] = (x - d1[i])/d2[i]
    return (data, d1, d2)
#------------------------------------------------------------------
# fill 2-D histograms
def fill(canvas, h, data, maxrows=2000):
    for index, (d, w) in enumerate(data):
        ih = 0
        for ii in xrange(len(d)):
            x = d[ii]
            for jj in xrange(ii+1, len(d)):
                y = d[jj]
                h[ih].Fill(x, y)
                ih += 1
        if index % 500 == 0:
            print index
            for ih in xrange(len(h)):
                canvas.cd(ih+1)
                h[ih].Draw('p')
            canvas.Update()
        if index > maxrows: break
#------------------------------------------------------------------    
def main():
    
    # set some standard graphics style (see histutil.py)
    setStyle()

    treename    = "HZZ4LeptonsAnalysisReduced"
    sigfilename = '../data/ntuple_4mu_VV.root'
    bkgfilename = '../data/ntuple_4mu_gg.root'

    sdata,d1,d2 = readData(sigfilename, treename)
    
    bdata,d1,d2 = readData(bkgfilename, treename, d1, d2)

    # create histograms
    hsig = []
    hbkg = []
    msize= 0.08
    for ii in xrange(len(VARS)):
        varx = VARS[ii]
        for jj in xrange(ii+1, len(VARS)):
            vary = VARS[jj]
            hs = mkhist2('hs_%d_%d'  % (ii, jj),
                         varx, vary,
                         50,-2, 2, 50,-2, 2)
            #hs.Sumw2()
            hs.SetMarkerColor(kCyan+1)
            hs.SetMarkerSize(msize)
            hs.SetMinimum(0)
            hs.GetXaxis().CenterTitle()
            hs.GetXaxis().SetTitleSize(0.08)
            hs.GetXaxis().SetTitleOffset(0.5)
            
            hs.GetYaxis().CenterTitle()
            hs.GetYaxis().SetTitleSize(0.08)
            hs.GetYaxis().SetTitleOffset(0.5)
            
            hsig.append(hs)
                                    
            hb = mkhist2('hb_%d_%d'  % (ii, jj),
                         varx, vary,
                         50,-2, 2, 50,-2, 2)
            #hb.Sumw2()
            hb.SetMarkerColor(kMagenta+1)
            hb.SetMarkerSize(msize)
            hb.SetMinimum(0)
            hb.GetXaxis().CenterTitle()
            hb.GetXaxis().SetTitleSize(0.08)
            hb.GetXaxis().SetTitleOffset(0.85)

            hb.GetYaxis().CenterTitle()
            hb.GetYaxis().SetTitleSize(0.08)
            hb.GetYaxis().SetTitleOffset(0.85)                        
            hbkg.append(hb)

    # fill histograms
    canvas = TCanvas('fig_variables', '', 10, 10, 800, 800)
    canvas.Divide(4, 4)

    fill(canvas, hsig, sdata)
    
    fill(canvas, hbkg, bdata)

    # plot histograms
    for ih in xrange(len(hsig)):
        canvas.cd(ih+1)
        hsig[ih].Draw('p')
        hbkg[ih].Draw('p same')
    canvas.Update()
    canvas.SaveAs('.png')
    
    sleep(2)
#------------------------------------------------------------------                
try:
    main()
except KeyboardInterrupt:
    print "\nciao!"
    
    
