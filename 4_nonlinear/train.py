#!/usr/bin/env python
#----------------------------------------------------------------------
# File: train.py
# Description: example of classification with TMVA
# Created: 01-June-2013 INFN SOS 2013, Vietri sul Mare, Italy, HBP
#   adapted for CMSDAS 2015 Bari HBP
# Updated: 01-June-2016 HBP for HAT@LPC June 9th 2016
#----------------------------------------------------------------------
import os, sys, re
from ROOT import *
#----------------------------------------------------------------------
def getTree(filename, treename):
    hfile = TFile(filename)
    if not hfile.IsOpen():
        sys.exit("** can't open file %s" % filename)
    tree = hfile.Get(treename)
    if tree == None:
        sys.exit("** can't find tree %s" % treename)
    return (hfile, tree)
#----------------------------------------------------------------------
#----------------------------------------------------------------------
def main():
    print "\n", "="*80
    print "\tclassification with TMVA"
    print "="*80

    treename = "HZZ4LeptonsAnalysisReduced"

    # get signal and background data for training/testing
    weightname  = "f_weight"  # name of event weight variable
    sigfilename = '../data/ntuple_4mu_VV.root'
    bkgfilename = '../data/ntuple_4mu_gg.root'
    
    sigFile, sigTree = getTree(sigfilename, treename)
    bkgFile, bkgTree = getTree(bkgfilename, treename)
    
    # everything is done via a TMVA factory
    outputFile = TFile("TMVA.root", "recreate")
    factory = TMVA.Factory("HATS", outputFile,
                           "!V:Transformations=I;N;D")

    # define input variables
    factory.AddVariable("f_deltajj", 'D')
    factory.AddVariable("f_massjj", 'D')

    # define from which trees data are to be taken
    factory.AddSignalTree(sigTree)
    factory.AddBackgroundTree(bkgTree)

    # define which variable is to be used for event-by-event weighting
    Lactual =   2.8 # 1/fb
    Lused   = 300.0 # 1/fb
    scale   = Lused / Lactual
    factory.SetWeightExpression("%f*f_weight" % scale)

    # define cuts to be applied to data and sample sizes
    # for training and testing
    cut    = 'f_massjj>0' # so that we have at least two jets
    counts = {'ntrain': 2500,
              'ntest':  1500}
    factory.PrepareTrainingAndTestTree(TCut(cut),
                                       TCut(cut),
                                       "nTrain_Signal=%(ntrain)d:"\
                                       "nTest_Signal=%(ntest)d:"\
                                       "nTrain_Background=%(ntrain)d:"\
                                       "nTest_Background=%(ntest)d:"\
                                       "!V" % counts)

    # define multivariate methods to be run
    # N: no transformation of inputs
    factory.BookMethod( TMVA.Types.kMLP,
                        "MLP",
                        "!H:!V:"\
                        "VarTransform=N:"\
                        "HiddenLayers=4:"\
                        "TrainingMethod=BFGS")

    # 1. Use AdaBoost algorithm
    # 2. Grow a forest of trees
    # 3. Each node must contain a count no less than 1% of the size of
    #    the training sample
    # 4. Consider 100 cuts in each dimension. For this example, this means
    #    for each decision node, consider 200 possible partitions from
    #    which the best cut is chosen.
    # 5. minimum percentage of events/leaf relative to training sample size
    factory.BookMethod( TMVA.Types.kBDT,
                        "BDT",
                        "!V:"\
                        "BoostType=AdaBoost:"\
                        "NTrees=100:"\
                        "MinNodeSize=1.0:"\
                        "nCuts=100")       
  
    factory.TrainAllMethods()  
    factory.TestAllMethods()
    factory.EvaluateAllMethods()
    
    outputFile.Close()
#----------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print "\nciao"
