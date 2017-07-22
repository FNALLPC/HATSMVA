#!/usr/bin/env python
#------------------------------------------------------------------
#- File: analyze.py
#- Description: plot results of training with TMVA
#  Created: 01-Jun-2013 INFN SOS 2013, Vietri sul Mare, Italy, HBP
#    adapt to CMSDAS 2015 Bari HBP
#    adapt to HATS@LPC 2016 June 9th
#------------------------------------------------------------------
import os, sys
from histutil import *
from time import sleep
from array import array
from ROOT import *
#------------------------------------------------------------------
FIRST_ROW=5000
OPTION='cont1'
#------------------------------------------------------------------
def getTreeWeights(code):
    # hack to get tree weights from TMVA class code
    get = re.compile('(?<=fBoostWeights[.]push_back\().*?(?=\);)', re.DOTALL)
    rec = get.findall(code)
    if len(rec) == 0:
        sys.exit('** cannot get weights from BDT class code')
    weights = []
    for t in rec:
        weights.append(atof(t))
    return weights

def readAndFill(filename, treename, h):
    print "==> reading %s" % filename
    # open ntuple (see histutil.py for implementation)
    ntuple = Ntuple(filename, treename, FIRST_ROW)
    total  = 0
    for event in ntuple :
        if not (event.f_massjj > 0): continue
        h.Fill(event.f_deltajj, event.f_massjj, event.f_weight)
        if total % 5000 == 0:
            print total
        total += 1            
    h.Scale(1.0/h.Integral())
#------------------------------------------------------------------
def readAndFillAgain(filename, treename, reader, which, c, h):
    ntuple = Ntuple(filename, treename, FIRST_ROW)
    total  = 0
    inputvars = vector('double')(2)
    isBDT  = which == 'BDT' 
    for event in ntuple:
        if not (event.f_massjj > 0): continue
        
        inputvars[0] = event.f_deltajj
        inputvars[1] = event.f_massjj
        
        # evaluate discriminant
        D = reader.GetMvaValue(inputvars)
        h.Fill(D, event.f_weight)

        if total % 5000 == 0:
            c.cd()
            h.Draw("hist")
            c.Update()
        total += 1
    h.Scale(1.0/h.Integral())
#------------------------------------------------------------------
def main():
    print "="*80
    # set up a standard graphics style	
    setStyle()

    xbins =  16 
    xmin  =   0.0
    xmax  =   8.0

    ybins =  16
    ymin  =   0.0
    ymax  =5000.0

    msize = 0.15
    
    fieldx = 'f_deltajj'; varx = '|#Delta#font[12]{#eta_{jj}}|'
    fieldy = 'f_massjj';  vary = '#font[12]{m_{jj}}'
    weightname  = "f_weight"
    treename    = "HZZ4LeptonsAnalysisReduced"
    sigfilename = '../data/ntuple_4mu_VV.root'
    bkgfilename = '../data/ntuple_4mu_gg.root'
    
    # pick discriminant
    if len(sys.argv) > 1:
        which = sys.argv[1]
    else:
        which = 'MLP'
    isBDT = which == 'BDT'
    
    # read in trained MLP class and compile it
    codename = 'weights/HATS_%s.class.C' % which
    print "=> compiling code:      %s" % codename
    code = open(codename).read()
    gROOT.ProcessLine(code)

    classname = 'Read%s(inputnames)' % which
    print "=> instantiating class: %s" % classname
    inputnames = vector('string')(2)
    inputnames[0] = fieldx
    inputnames[1] = fieldy
    reader = eval(classname)

    # get tree weights from BDT class code.
    # Assuming the AdaBoost algorithm was used, we
    # need the summed weights in order transform BDT
    # output to a probability:
    # D = 1/(1 + exp(-2*summedalpha*D))
    summedalpha = None
    if isBDT:
        print "=> extracting alphas"
        weights = getTreeWeights(code)
        summedalpha = sum(weights)

    # ---------------------------------------------------------
    # make 2-D surface plot
    # ---------------------------------------------------------

    print "=> plotting"
    c  = TCanvas("fig_VV_gg_%s" % which, "", 10, 10, 800, 800)
    # divide canvas canvas along x-axis
    c.Divide(2, 2)

    # Fill signal histogram
    hsig = mkhist2('hsig', varx, vary,
                   xbins, xmin, xmax,
                   ybins, ymin, ymax)
    hsig.Sumw2()
    hsig.SetMarkerSize(msize)
    hsig.GetYaxis().SetTitleOffset(2.10)
    hsig.SetMinimum(0)    
    hsig.SetMarkerSize(msize)
    hsig.SetMarkerColor(kCyan+1)
          
    readAndFill(sigfilename, treename, hsig)

    # Fill background histogram
    hbkg = mkhist2('hbkg', varx, vary,
                   xbins, xmin, xmax,
                   ybins, ymin, ymax)
    hbkg.Sumw2()
    hbkg.SetMarkerSize(msize)
    hbkg.GetYaxis().SetTitleOffset(2.10)
    hbkg.SetMinimum(0)        
    hbkg.SetMarkerSize(msize)
    hbkg.SetMarkerColor(kMagenta+1)        
    readAndFill(bkgfilename, treename, hbkg)

    # make some plots

    xpos = 0.30
    ypos = 0.85
    tsize= 0.05

    # --- signal

    c.cd(1)
    hsig.SetMinimum(0)
    hsig.Draw('p')
    hsig.Draw('same '+OPTION)
    s1 = Scribe(xpos, ypos, tsize)
    s1.write('VV #rightarrow H #rightarrow ZZ #rightarrow 4l')
    c.Update()

    # --- background

    c.cd(2)
    hbkg.Draw('p')
    hbkg.Draw('same '+OPTION)
    hbkg.Draw('cont1 same')
    s2 = Scribe(xpos, ypos, tsize)
    s2.write('gg #rightarrow H #rightarrow ZZ #rightarrow 4l')
    c.Update()        

    # --- p(S|x) = p(x|S) / [p(x|S) + p(x|B)]

    hD = hsig.Clone('hD')
    hsum = hsig.Clone('hSum')
    hsum.Add(hbkg)
    hD.Divide(hsum)
    hD.SetMinimum(0)
    hD.SetMaximum(1)
    hD.GetYaxis().SetTitleOffset(2.10)

    c.cd(3)
    hD.Draw(OPTION)
    s3 = Scribe(xpos, ypos, tsize)
    s3.write('D(%s, %s) (actual)' % (varx, vary))
    c.Update()

    # ---------------------------------------------------------
    h1 = mkhist2("h1", varx, vary,
                 xbins, xmin, xmax,
                 ybins, ymin, ymax)                 
    h1.SetMinimum(0)
    h1.SetMaximum(1)
    h1.GetYaxis().SetTitleOffset(2.10)
    
    # compute discriminant at a grid of points

    xstep = (xmax-xmin)/xbins
    ystep = (ymax-ymin)/ybins
    inputvars = vector('double')(2)
    for i in xrange(xbins):
        x = xmin + (i+0.5)*xstep
        inputvars[0] = x
        for j in xrange(ybins):
            y = ymin + (j+0.5)*ystep
            inputvars[1] = y
            D = reader.GetMvaValue(inputvars)

            # need to transform BDT output so that
            # we arrive at an apples to apples
            # comparison.
            if isBDT:
                D = 1.0/(1 + exp(-summedalpha*D))            
            h1.Fill(x, y, D)

    # plot MVA approximation to discriminant
    c.cd(4)
    h1.Draw(OPTION)
    s4 = Scribe(xpos, ypos, tsize)
    s4.write('D(%s, %s) (%s)' % (varx, vary, which))
    c.Update()
    c.SaveAs(".png")
    # ---------------------------------------------------------
    # plot distributions of D
    # ---------------------------------------------------------
    c1  = TCanvas("fig_VV_gg_D_%s" % which, "",
                  710, 310, 500, 500)

    xmin =-1 if isBDT else 0
    xmax = 1
    hs = mkhist1("hs", "D(%s, %s)" % (varx, vary), "", 50, xmin, xmax)
    hs.SetFillColor(kCyan+1)
    hs.SetFillStyle(3001)
    readAndFillAgain(sigfilename, treename, reader, which, c1, hs)
    sleep(2)

    hb = mkhist1("hb", "D(%s, %s)" % (varx, vary), "", 50, xmin, xmax)
    hb.SetFillColor(kMagenta+1)
    hb.SetFillStyle(3001)
    readAndFillAgain(bkgfilename, treename, reader, which, c1, hb)

    c1.cd()
    hb.Draw('hist')
    hs.Draw("hist same")
    c1.Update()
    c1.SaveAs(".png")
    sleep(4)
#----------------------------------------------------------------------
main()
