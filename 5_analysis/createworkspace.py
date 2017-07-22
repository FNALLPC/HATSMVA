#!/usr/bin/env python
#-------------------------------------------------------------
# File: createworkspace.py
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
def createWorkspace(wsname, wsfilename):
    # The most convenient way to use RooFit/RooStats is to 
    # make a workspace so that we can use its factory method
    # and write the probability model and data to a Root file
    wspace = RooWorkspace(wsname)

    #-----------------------------------------------------
    # Create parameters
    #
    # Use the factory method of the RooWorkspace to create
    # parameters
    #
    # syntax:
    #        <name>[value, min-value, max-value]
    #-----------------------------------------------------
    # observations
    params = [('N',     12,     0,  50),
              
    # background 1 estimate        
              ('B1',     0.04,  0,  1.0),
              ('dB1',    0.04,  0,  1.0),
    # background 2 estimate
              ('B2',     5.3,   0,  15),
              ('dB2',    0.53,  0,   5),
    # signal estimate
              ('S',      6.7,   0,  25),
              ('dS',     0.67,  0,   5),
    # nuisance parameters
              ('b1',     0.04,  1e-3,   1),
              ('b2',     5.3,   1e-3,  10),
              ('s',      6.7,   1e-3,  25),
    # parameter of interest
              ('mu',     1.0,   0,   4)]

    for t in params:
        cmd = '%s[%f, %f, %f]' % t
        wspace.factory(cmd)        
    wspace.var('mu').SetTitle('#mu')

    # fix all background and signal constants (B1 to dS)
    for t in params[1:-4]: # skip first row and last four
        name = t[0]
        wspace.var(name).setConstant()

    #-----------------------------------------------------
    # Create expressions
    #
    # syntax:
    #        expr::<name>("expression", var1, var2, ...)
    #-----------------------------------------------------
    # Given B +/-dB, we define an effective count Q and an
    # effective scale factor q by
    #      Q  = q B
    # sqrt(Q) = q dB,
    # which yield
    #      Q  = (B/dB)^2
    #      q  = B/dB^2
    express = ['Q1("(B1/dB1)^2", B1, dB1)',
               'q1("B1/dB1^2",   B1, dB1)',
               
               'Q2("(B2/dB2)^2", B2, dB2)',
               'q2("B2/dB2^2",   B2, dB2)',
                              
               'Q("(S/dS)^2", S, dS)',
               'q("S/dS^2",   S, dS)',

               'q1b1("q1*b1", q1, b1)',
               'q2b2("q2*b2", q2, b2)',
               'qs("q*s", q, s)',
               
               'n("mu*s + b1 + b2", mu, s, b1, b2)']
        
    for t in express:
        cmd = 'expr::%s' % t
        wspace.factory(cmd)

    print '\neffective counts and scale factors'
    print 'Q1 = %8.2f, q1 = %8.2f' % (wspace.function('Q1').getVal(),
                                      wspace.function('q1').getVal())

    print 'Q2 = %8.2f, q2 = %8.2f' % (wspace.function('Q2').getVal(),
                                      wspace.function('q2').getVal())

    print 'Q  = %8.2f, q  = %8.2f' % (wspace.function('Q').getVal(),
                                      wspace.function('q').getVal())

    #-----------------------------------------------------
    # Create pdfs
    #
    # syntax:
    #        pdf_name::<name>(var1, var2, ...)
    #
    # where the "Roo" prefix is dropped in pdf_name, e.g.
    #-----------------------------------------------------
    pdfs = [('Poisson', 'pN',  '(N, n)'),
            # allow non-integer Q (we are reiterpreting
            # the Poisson as a gamma density in the
            # parameter b1, b2, or s.           
            ('Poisson', 'pQ1', '(Q1, q1b1, 1)'), 
            ('Poisson', 'pQ2', '(Q2, q2b2, 1)'),
            ('Poisson', 'pS',  '(Q,  qs,   1)')]

    prodpdf = ''
    for t in pdfs:
        cmd = '%s::%s%s' % t
        print cmd
        wspace.factory(cmd)
        name = t[1]
        prodpdf += "%s, " % name
    prodpdf = prodpdf[:-2] # remove last ", "
    
    # multiply the pdfs together. use upper case PROD to
    # do this
    wspace.factory('PROD::model(%s)' % prodpdf)

    # create a prior, since one is needed for Bayesian
    # calculations
    wspace.factory('Uniform::prior({mu, s, b1, b2})')

    #-----------------------------------------------------
    # Define a few useful sets. Now we need to decide
    # whether or not to include B and S in the set obs of
    # observations. 
    #-----------------------------------------------------
    sets = [('obs',  'N'),           # observations
            ('poi',  'mu'),          # parameter of interest
            ('nuis', 's,b1,b2')]     # nuisance parameters (leave no spaces)
    for t in sets:
        name, parlist = t
        wspace.defineSet(name, parlist)
    
    #-----------------------------------------------------        
    # Create a dataset
    #-----------------------------------------------------    
    data = RooDataSet('data', 'data', wspace.set('obs'))
    data.add(wspace.set('obs'))
    # import dataset into workspace
    # need last argument to workaround a PyROOT "feature".
    # the last argument ensures the correct version of
    # the import method is called.
    getattr(wspace, 'import')(data, RooCmdArg())
        
    #-----------------------------------------------------
    # Create model configuration. This is needed for the
    # statistical analyses
    #-----------------------------------------------------
    cfg = RooStats.ModelConfig('cfg')
    cfg.SetWorkspace(wspace)
    cfg.SetPdf(wspace.pdf('model'))
    cfg.SetPriorPdf(wspace.pdf('prior'))
    cfg.SetParametersOfInterest(wspace.set('poi'))
    cfg.SetNuisanceParameters(wspace.set('nuis'))

    # import model configuration into workspace
    getattr(wspace, 'import')(cfg)

    wspace.Print()
    
    # write out workspace
    wspace.writeToFile(wsfilename)
#------------------------------------------------------------------
def main():
    # Suppress all messages except those that matter
    msgservice = RooMsgService.instance()
    msgservice.setGlobalKillBelow(RooFit.WARNING)
    print "="*80

    createWorkspace('HATS@LPC', 'HATSworkspace.root')
#------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print
    print "ciao!"
    print



