#!/usr/bin/env python
#------------------------------------------------------------------------------
# File: plottrees.py
# Description: plot 2-D trees (experimental code!). Do a git pull to get
#              latest version of histutil
#
# Created: 31-May-2013 INFN SOS 2013, Salerno, Italy, HBP
#  adapt to HATS@LPC 2016 HBP
#-----------------------------------------------------------------------------
import os, sys
from math import *
from string import *
from time import sleep
from array import array
from histutil import *
from ROOT import *
#-----------------------------------------------------------------------------
def fixhist(hb): 
    hb.GetXaxis().CenterTitle()
    hb.GetXaxis().SetTitleSize(0.12)
    hb.GetXaxis().SetTitleOffset(0.55)

    hb.GetYaxis().CenterTitle()
    hb.GetYaxis().SetTitleSize(0.12)
    hb.GetYaxis().SetTitleOffset(0.55)
#----------------------------------------------------------------------------
def main():
    print "="*80
    filename  = 'weights/HATS_BDT.class.C'
    treename  = 'HZZ4LeptonAnalysisReduced'
    weightname= 'f_weight'
        
    # create a Python wrapper of the TMVA BDT class
    bdt = BDT(filename)

    xbins =  16 
    xmin  =   0.0
    xmax  =   8.0

    ybins =  16
    ymin  =   0.0
    ymax  =5000.0

    msize = 0.15
    
    fieldx = 'f_deltajj'; varx = '|#Delta#font[12]{#eta_{jj}}|'
    fieldy = 'f_massjj';  vary = '#font[12]{m_{jj}}'

    # draw individual trees
    pixels = 300
    nx = 2
    ny = 2
    nn = nx*ny

    setStyle()

    c = TCanvas('fig_%dtrees' % nn, 'trees', 10, 10, nx*pixels, ny*pixels)
    c.Divide(nx, ny)
    hist = []
    line = []
    for ii in xrange(nx*ny):
        c.cd(ii+1)
        h = bdt.plot2d(ii, 'h%2.2d' % ii,
                       varx, vary,
                       xmin, xmax, ymin, ymax)        
        hist.append( h )
        h[0].Draw('col')
        h[1].Draw('same')

    c.SaveAs(".png")
    
    # draw 2D plot with increasing numbers of trees
    c1 = TCanvas('fig_forest', 'trees', 820, 10, 600, 600)
    c1.Divide(2,2)
    nx =100
    ny =100
    xstep = (xmax-xmin) / nx
    ystep = (ymax-ymin) / ny
        
    ntrees = [2, 8, 32, 100]
    hh = []
    firstTree = 0
    for ii in xrange(len(ntrees)):
        hname = "h%4.4d" % ntrees[ii]
        lastTree = ntrees[ii]-1
        hh.append(mkhist2(hname, varx, vary,
                          nx, xmin, xmax, ny, ymin, ymax))
        for ix in xrange(nx):
            x = xmin + ix * xstep
            for iy in xrange(ny):
                y = ymin + iy * ystep
                vtuple = (x, y)
                z = bdt(vtuple, firstTree, lastTree)
                hh[-1].SetBinContent(ix+1, iy+1, z);
        c1.cd(ii+1)
        hh[-1].Draw('col')
        addTitle('%5d trees' % ntrees[ii], 0.06)
        c1.Update()
    c1.SaveAs('.png')
    sleep(10)

    ## # make an animiated gif
    ## os.system('rm -rf fig_manytrees.gif')
    ## c1 = TCanvas('fig_manytrees', 'many trees', 610, 310, 2*pixels, 2*pixels)
    ## option = 'col'
    ## hbdt = []
    ## for ii in xrange(100):
    ##     print ii
    ##     hbdt.append(bdt.plot(ii,
    ##                       'hmtree',
    ##                       'SO_{2}', 'alcohol',
    ##                       xmin, xmax, ymin, ymax))
    ##     #fixhist(hbdt[-1])
    ##     hbdt[-1].Draw('col')
    ##     c1.Update()
    ##     c1.Print('fig_manytrees.gif+10')
#-----------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print "\nciao"
