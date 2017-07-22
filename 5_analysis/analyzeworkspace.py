#!/usr/bin/env python
#-------------------------------------------------------------
# File: analyzeworkspace.py
# Description: Single count analysis with two backgrounds and
#              one signal.
#
# Example:
#   Known                 Unknowns
#   N  = 12               n   = mu * s + b1 + b2
#   B1 =  0.04 +/- 0.04   b1   (ZZ background)
#   B2 =  5.30 +/- 0.53   b2   (gg -> H -> ZZ -> 4mu)
#   S  =  6.70 +/- 0.67   s    (VV -> H -> ZZ -> 4mu)
#
# where mu is the signal strength.
#
# Created: 18-Dec-2015 CMSDAS 2016, LPC Fermilab HBP
#          08-Jun-2016 Adapted to HATS@LPC 2016
#-------------------------------------------------------------
import os,sys,re
from time import sleep
from math import *
from ROOT import *
#-------------------------------------------------------------
def check(o, message):
    if o == None:
        sys.exit(message)
#-------------------------------------------------------------        
def analyzeWorkspace(wsname, wsfilename):

    # Open workspace file
    wsfile = TFile(wsfilename)
    if not wsfile.IsOpen():
        sys.exit("*** can't open file %s" % wsfilename)

    # Get workspace
    wspace = wsfile.Get(wsname)
    check(wspace, "*** can't access workspace %s" % wsname) 

    # Get data
    data = wspace.data('data')
    check(data, "*** can't access workspace data")

    # Get model configuration    
    cfg  = wspace.obj('cfg')
    check(cfg, "*** can't access model configuration cfg")

    #-----------------------------------------------------    
    # Fit model to data
    #-----------------------------------------------------
    results = wspace.pdf('model').fitTo(data, RooFit.Save())
    results.Print()
    
    #-----------------------------------------------------    
    # Compute interval based on profile likelihood
    #-----------------------------------------------------
    # suppress some (apparently) innocuous warnings
    msgservice = RooMsgService.instance()
    msgservice.setGlobalKillBelow(RooFit.FATAL)
        
    print 'compute 68% interval using profile likelihood'
    plc = RooStats.ProfileLikelihoodCalculator(data, cfg)
    CL  = 0.683
    plc.SetConfidenceLevel(CL)
    plcInterval= plc.GetInterval()
    lowerLimit = plcInterval.LowerLimit(wspace.var('mu'))
    upperLimit = plcInterval.UpperLimit(wspace.var('mu'))

    print '\tPL %4.1f%s CL interval = [%5.2f, %5.2f]' % \
      (100*CL, '%', lowerLimit, upperLimit)

    # compute a 95% upper limit on mu by
    # computing a 90% central interval and
    # ignoring the lower limit
    CL = 0.90
    plc.SetConfidenceLevel(CL)
    plcInterval = plc.GetInterval()
    upperLimit = plcInterval.UpperLimit(wspace.var('mu'))

    CL = 0.95
    print '\tPL %4.1f%s upper limit = %5.2f\n' % \
      (100*CL, '%', upperLimit)      

    # plot it
    plcplot = RooStats.LikelihoodIntervalPlot(plcInterval)
    plccanvas = TCanvas('fig_PL', 'plc', 10, 10, 500, 500)
    plcplot.Draw()
    plccanvas.Update()

    # In the frequentist approach, the goal is to reject an hypothesis, 
    # the so-called null (that is, no effect) hypothesis based on a
    # measure called a test statistic with the property that extreme (e.g.,
    # large) values would tend to cast doubt on the veracity of the null
    # hypothesis.
    #
    # The following is a standard test statistic:
    #   q(mu) = -2*log[Lp(mu) / Lp(mu_hat)]
    # where Lp(mu) is the profile likelihood and mu_hat is the
    # maximum likelihood estimate (MLE) of mu (i.e., the best fit
    # value). The null hypothesis in this case is mu=0. In effect, we are
    # comparing the hypothesis that best fits the data, namely, the hypothesis
    # mu = mu_hat to the no effect (i.e. null) hypothesis, mu = 0.
    #
    # The reason q(mu) is useful is that if there is no effect, then the
    # probability density of q(0) is approximately a chi-squared density
    # with (in this case) one degree of freedom. Consequently, Z = sqrt(q(0))
    # is. approximately, the number of standard deviations the result is
    # away from the null hypothesis. The approximation becomes better as the
    # count or sample size grows without limit. In practice, for Poisson
    # problems, the the approximation works well even for small counts,
    # provided that most of the time, mu_hat is not on the boundary of mu.
    #
    # The convention in particle physics is to declare victory if Z >= 5 and,
    # sometimes, state that one has evidence if Z >= 3.
    #
    # In order to compute q(0), one proceeds as follows:
    # 1. get the negative log-profilelikelihood ratio (-log[Lp(mu)/Lp(mu_hat)])
    nlplratio = plcInterval.GetLikelihoodRatio()
    # 2. set the value of mu to zero (the no effect hypothesis)
    wspace.var('mu').setVal(0)
    # 3. and compute q(0)
    q0 = 2*nlplratio.getVal()
    # 4. compute a Z value
    Z  = sqrt(q0)
    print "\tZ-value = %8.1f\n" % Z

    #-----------------------------------------------------    
    # Compute interval based on Bayesian calculator
    #-----------------------------------------------------
    print 'compute interval using Bayesian calculator'
    print "\t\tplease be patient...!"
    bc = RooStats.BayesianCalculator(data, cfg)
    CL  = 0.683
    bc.SetConfidenceLevel(CL)
    bcInterval = bc.GetInterval()
    lowerLimit = bcInterval.LowerLimit()
    upperLimit = bcInterval.UpperLimit()

    print '\tBayes %4.1f%s CL interval = [%5.2f, %5.2f]' % \
      (100*CL, '%', lowerLimit, upperLimit)

    NP = 50 # number of points at which to calculate posterior density
    # calculate posterior density
    bc.SetScanOfPosterior(NP)
    bcplot = bc.GetPosteriorPlot()
    bccanvas = TCanvas('fig_Bayes', 'Bayes', 500, 10, 850, 400)
    bccanvas.Divide(2, 1)
    bccanvas.cd(1)
    bcplot.Draw()
    bccanvas.Update()

    # compute a 95% upper limit on mu
    CL  = 0.950
    bc.SetConfidenceLevel(CL)
    # 0   = upper limit
    # 0.5 = central limits (default)
    # 1   = lower limit
    bc.SetLeftSideTailFraction(0)
    bcInterval = bc.GetInterval()
    upperLimit = bcInterval.UpperLimit()

    print '\tBayes %4.1f%s upper limit = %5.2f\n' % \
      (100*CL, '%', upperLimit)

    print "\t\tplease be patient...computing posterior density again!"      
    # calculate posterior density
    bc.SetScanOfPosterior(NP)
    bcplot2 = bc.GetPosteriorPlot()
    bccanvas.cd(2)
    bcplot2.Draw()
    bccanvas.Update()

    # save canvases
    plccanvas.SaveAs('.png')
    bccanvas.SaveAs('.png')
    
    sleep(5)  
#------------------------------------------------------------------
def main():
    # Suppress all messages except those that matter
    msgservice = RooMsgService.instance()
    msgservice.setGlobalKillBelow(RooFit.WARNING)
    print "="*80

    analyzeWorkspace('HATS@LPC', 'HATSworkspace.root')
#------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print
    print "ciao!"
    print



