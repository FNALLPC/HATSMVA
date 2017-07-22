#!/usr/bin/env python
#------------------------------------------------------------------------------
# File: maketree.py
# Description: make a ROOT tree
# Created: 22 Sep 2010 Harrison B. Prosper & Sezen Sekmen
#          06 Jun 2016 HBP adapted for HATS@LPC 2016
#------------------------------------------------------------------------------
import os, sys, re
from string import *
from histutil import *
from time import sleep, ctime
from ROOT import *
#------------------------------------------------------------------------------
def nameonly(s):
    import posixpath
    return posixpath.splitext(posixpath.split(s)[1])[0]
#------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------
def getVarnames(code):
    # hack to get variables from TMVA class code
    # get lines from NVar to NSpec
    get = re.compile('(?s)NVar(.*?)NSpec', re.M)
    rec = get.findall(code)
    if len(rec) == 0:
        sys.exit('** cannot get variable names from class code')
    # split at newline (skip first line)
    rec = split(strip(rec[0]), '\n')[1:]

    # first field of each lone is weight
    varnames = vector('string')()
    for t in rec:
        t = split(t)
        varnames.push_back(t[0])
    return varnames
#------------------------------------------------------------------------------
def readData(filename, treename, MLP, BDT, varnames, Lumi, summedalphas):
    # mass4l window (just to check counts)
    lower  = 110 # GeV
    upper  = 136 # GeV
        
    # scale to desired lumi (original ntuples scaled to 2.8/fb)
    if   find(filename, '_data') > 0:
        scale = 1 # don't scale real data
    elif find(filename, '_bkg') > 0:
        # the 1.5 makes the background match the observed count in the
        # sidebands m4l < 110 or m4l > 136 GeV.
        scale = 1.50406* Lumi / 2.8 
    else:
        scale = Lumi / 2.8
     
    print "=> reading file %s, scale weights by %8.1f" % (filename, scale)
    ntuple = Ntuple(filename, treename)

    total   = 0
    i_weight= 0.0
    w_weight= 0.0
    m_weight= 0.0
    inputvars = vector('double')(2)
    records= []
    for event in ntuple:
        w = scale * event.f_weight
        i_weight += w
        
        # impose window cut
        if event.f_mass4l < lower: continue
        if event.f_mass4l > upper: continue
        w_weight += w

        # require at least two jets
        if event.f_massjj <= 0:    continue
        m_weight += w
        
        if total % 5000 == 0:
            print "\t",total

        total  += 1
        
        # evaluate discriminants        
        for i, varname in enumerate(varnames):
            inputvars[i] = eval('event.%s' % varname)
            
        D_MLP = MLP.GetMvaValue(inputvars)
        D_BDT = BDT.GetMvaValue(inputvars)
        D_BDT = 1.0/(1 + exp(-2*summedalphas*D_BDT))
                    
        records.append((D_MLP, D_BDT, event.f_D_bkg, w))

    print
    print "cut flow"
    print "\tno cuts:                 %10.3f (%d)" % (i_weight, len(ntuple))
    print "\twindow cut:              %10.3f" % w_weight
    print "\tmassjj > 0 cut:          %10.3f (%d)\n" % (m_weight, total)
    print
    return records
#------------------------------------------------------------------------------
def makeTree(filename, treename, records, complevel=2):
    print "=> writing to file %s" % filename
    
    # open output root file
    tfile = TFile(filename, "recreate")
    tfile.SetCompressionLevel(complevel)
    ttree = TTree(treename, "%s created: %s" % (treename, ctime()))
    
    # ---------------------------------------------------
    # create branches
    # ---------------------------------------------------
    # variable names in ntuple to be created
    varnames = ['D_VVgg_MLP', 'D_VVgg_BDT', 'D_bkg', 'weight']
    
    # make a struct of variables
    rec = 'struct Bag {'
    for varname in varnames:
        rec += 'double %s;' % varname
    rec += "};"
    
    # compile struct
    gROOT.ProcessLine(rec)
    # and make it visible to Python
    from ROOT import Bag
    bag = Bag()
    
    # OK, actually create branches. note use of AddressOf
    branch = []
    for varname in varnames:
        branch.append( ttree.Branch(varname,
                                    AddressOf(bag, varname),
                                    "%s/D" % varname) )

    # fill tree
    for index, record in enumerate(records):
        if index % 5000 == 0: print '\t',index
        for ii, x in enumerate(record):
            # note use of __setattr__(name, value) to set
            # attributes of struct
            bag.__setattr__(varnames[ii], x)
            
        tfile.cd()
        ttree.Fill()
        
    tfile.Write("", TObject.kOverwrite)
#------------------------------------------------------------------------------
def main():
    print "\n\tmaketree.py\n"
    
    if len(sys.argv) < 2:
        sys.exit('''
Usage:
       maketree.py input-root-file Lumi [300/fb]
        ''')

    filename = sys.argv[1]
    if not os.path.exists(filename):
        sys.exit('** file %s not found' % filename)

    if len(sys.argv) > 2:
        Lumi = atof(sys.argv[2])
    else:
        Lumi = 300.0 # 1/fb
        
    treename = "HZZ4LeptonsAnalysisReduced"
    
    # read and compile MLP class
    codename = '../4_nonlinear/weights/HATS_MLP.class.C'
    if not os.path.exists(codename):
        sys.exit('** file %s NOT found\n'\
                 '** run ../4_nonlinear/train.py to create it\n' % codename)
    print "=> compiling %s" % codename
    code = open(codename).read()
    gROOT.ProcessLine(code)

    # read and compile BDT class
    codename = '../4_nonlinear/weights/HATS_BDT.class.C'
    if not os.path.exists(codename):
        sys.exit('** file %s NOT found\n'\
                 '** run ../4_nonlinear/train.py to create it\n' % codename)    
    print "=> compiling %s" % codename
    BDTcode = open(codename).read()
    gROOT.ProcessLine(BDTcode)
    
    # instantiate discriminants. Annoyingly, we need to pass
    # the names to it first. Let's just extract them from the code
    varnames = getVarnames(BDTcode)
    MLP  = ReadMLP(varnames)
    BDT  = ReadBDT(varnames)

    # load data into memory
    summedalphas = sum(getTreeWeights(BDTcode))
    records  = readData(filename, treename, MLP, BDT, varnames,
                        Lumi, summedalphas)

    filename = '%s.root' % replace(nameonly(filename), 'ntuple', 'd')
    makeTree(filename, treename, records)
    
    print '\ndone!\n'
#------------------------------------------------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print "\nciao!"
    
