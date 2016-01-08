#!/user/bin/env python

import os,sys
import array
import math
import ROOT
from scipy.stats.mstats import mquantiles,moment
import numpy

"""
converts an array of estimators to an histogram and fits a gaussian
"""
def fitResolution(vals):
    xmin,xfitmin,xfitmax,xmax=mquantiles(vals,prob=[0.01,0.1,0.9,0.99])
    hresol=ROOT.TH1F('hresol','hresol',100,xmin,xmax)
    for x in vals: hresol.Fill(x)

    hresol.Draw('e1')
    hresol.Fit('gaus','RQ+','',xfitmin,xfitmax)
    fitFunc=hresol.GetFunction('gaus')
    sigma,sigmaUnc=fitFunc.GetParameter(2),fitFunc.GetParError(2)
    return sigma,sigmaUnc,hresol

"""
"""
def main():

    #ROOT configuration
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)

    siData=[ 
        (133,'timesummary_133.root'),
        (211.5,'timesummary_211.root'),
        (285,'timesummary_285.root')
        ]

    enEstProb=[0.1,0.16,0.5,0.68,0.90,0.95,0.99,0.999]
    enEstimators=['calibEn','AoverN']

    #loop over reference channels 
    for ch1,ch0 in [ (1,0),(0,2),(1,2) ]:

        print 'Analysing (ch1,ch0)=(%d,%d)'%(ch1,ch0)

        fOut=ROOT.TFile.Open('timeresol_%dvs%d_summary.root'%(ch1,ch0),'RECREATE')

        #loop over the different data
        for siWidth, fname in siData:
            print '\t si width=',siWidth,' analysed from',fname
    
            fOut.cd()
            outDir=fOut.mkdir('si%d'%siWidth)

            fIn=ROOT.TFile.Open(fname)
            fIn.cd()
            ttree=fIn.Get('ttree')
    
            #loop over different energy estimators
            for est in enEstimators:
                
                #compute the quantiles of the energy
                enEstVals=[]
                for i in xrange(0,ttree.GetEntriesFast()):
                    ttree.GetEntry(i)
                    e=getattr(ttree,est)
                    enEst=e[ ch1 ]
                    if ch0!=2: enEst=0.5*(e[ch1]+e[ch0])
                    if enEst<3 : continue
                    enEstVals.append(enEst)
                enEstQuantiles=mquantiles(enEstVals,prob=enEstProb)
            
                #time estimators to consider
                tEstimators=['t_max','t_max_frac50','t_max_frac30']
                if ch0!=2:
                    tEstimators += ['t_max_m_t_max_frac30', 't_max_frac50_m_t_max_frac30']

                #prepare arrays to store time differences
                tEstVals={}            
                for tEst in tEstimators: 
                    tEstVals[tEst]={}
                    for iq in xrange(0,len(enEstQuantiles)-1):
                        tEstVals[tEst][iq]=[]

                #project the time resolution for the quantile-defined ranges of the energy spectrum
                for i in xrange(0,ttree.GetEntriesFast()):
                    ttree.GetEntry(i)
                    e=getattr(ttree,est)
                    enEst=e[ ch1 ]
                    if ch0!=2: enEst=0.5*(e[ch1]+e[ch0])

                    #map to the quantile
                    estQ=-1
                    for iq in xrange(0,len(enEstQuantiles)-1):
                        if enEst>=enEstQuantiles[iq] and enEst<enEstQuantiles[iq+1]:
                            estQ=iq
                            break
                    if estQ<0 : continue
            
                    #time estimate
                    for tEst in tEstimators:
                        tDiff=0
                        if '_m_' in tEst:
                            tEst1,tEst0=tEst.split('_m_')
                            t1=getattr(ttree,tEst1)
                            t0=getattr(ttree,tEst0)
                            tDiff=(t1[ch1]-t0[ch1])-(t1[ch0]-t0[ch0])
                        else:
                            t=getattr(ttree,tEst)
                            tDiff=t[ch1]-t[ch0]
                        tEstVals[tEst][estQ].append( tDiff )
                        
                #determine bias and resolution
                grBias,grResol,hresol={},{},{}
                for tEst in tEstimators:
                
                    #prepare graphs
                    grBias[tEst]=ROOT.TGraphErrors()
                    grBias[tEst].SetName('bias_%s_%s'%(tEst,est))
                    grResol[tEst]=grBias[tEst].Clone('resol_%s_%s'%(tEst,est))
                    hresol[tEst]=[]

                    for iq in tEstVals[tEst]:
                    
                        #bias is determined from the median of the values
                        tEstMedian=mquantiles(tEstVals[tEst][iq],prob=[0.5])[0]
                        tEstMedianUnc=1.253*numpy.std(tEstVals[tEst][iq])/math.sqrt(len(tEstVals[tEst][iq]))

                        #determine the resolution after correcting for the median bias
                        for idx in xrange(0,len(tEstVals[tEst][iq])): tEstVals[tEst][iq][idx]-=tEstMedian
                        tEstResol,tEstResolUnc,h=fitResolution(tEstVals[tEst][iq])
                        hresol[tEst].append(h)
                        hresol[tEst][-1].SetDirectory(0)
                        hresol[tEst][-1].SetName('resolfit_%s_%s_%d'%(tEst,est,iq))

                        xmin,xmax=enEstQuantiles[iq],enEstQuantiles[iq+1]
                        xcen=0.5*(xmax+xmin)
                        xdiff=0.5*(xmax-xmin)
                        if tEstMedian!=0:
                            if tEstMedianUnc/tEstMedian<0.2:
                                np=grBias[tEst].GetN()
                                grBias[tEst].SetPoint(np,xcen,tEstMedian)
                                grBias[tEst].SetPointError(np,xdiff,tEstMedianUnc)
                        if tEstResol>0:
                            if tEstResolUnc/tEstResol<0.2:
                                np=grResol[tEst].GetN()
                                grResol[tEst].SetPoint(np,xcen,tEstResol)
                                grResol[tEst].SetPointError(np,xdiff,tEstResolUnc)

                #dump results to file                
                outDir.cd()
                for tEst in tEstimators:
                    grBias[tEst].Write()
                    grResol[tEst].Write()
                    for h in hresol[tEst]:
                        h.SetDirectory(outDir)
                        h.Write()
            
            #all done with this input
            fIn.Close()
        
        print 'Summary plots saved in ',fOut.GetName()
        fOut.Close()
    

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
