# /usr/bin/python

import sys
import os
import ROOT

from FTTB215Conditions import CONDITIONS,DIODES

COLORS ={11:1,13:ROOT.kGray+1}
MARKERS={(11,2):22,(11,3):26,(13,2):20,(13,3):24}

def summarizeResultsFor(ntuple,chargeEst,siPadSel,fitType,title,var,varUnc=None,outDir='./'):

    grColl={}
    for i in xrange(0,ntuple.GetEntries()):

        ntuple.GetEntry(i)

        if chargeEst!=ntuple.chargeEstimator : continue
        if siPadSel>0 and ntuple.siPad!=siPadSel:continue
        if fitType != ntuple.fitType : continue
        if ntuple.prob<0.01 : continue
        if ntuple.adc2mipUnc/ntuple.adc2mip>0.2: continue
        #if ntuple.adc2mipUnc/ntuple.adc2mip<0.05 : continue

        #run conditions
        runConds = [x for x in CONDITIONS if x[2]==ntuple.runNumber][0]
        beam     = runConds[0]
        beamEn   = runConds[1]
        diode    = runConds[5]
        siwidth  = DIODES[diode][3]

        varVal=getattr(ntuple,var)
        varValUnc=0
        try:
            varValUnc=getattr(ntuple,varUnc)
        except:
            pass

        if not beam in grColl:
            grColl[beam]=ROOT.TGraphErrors();
            grColl[beam].SetName('%s_beam%d_si%s_fit%d_chargeEst%d' % (var,beam,siPadSel,fitType,chargeEst) )
            grColl[beam].SetMarkerStyle( MARKERS[(beam,ntuple.siPad)] )
            grColl[beam].SetMarkerColor( COLORS[beam] )
            grColl[beam].SetLineColor( COLORS[beam] )
            grColl[beam].SetLineWidth(2)
            grColl[beam].SetFillStyle(0)
                
        np=grColl[beam].GetN()
        grColl[beam].SetPoint(np,siwidth,varVal)    
        grColl[beam].SetPointError(np,0,varValUnc)    
  
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetRightMargin(0.02)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.12)
    drawOpt='ap'
    #leg1=ROOT.TLegend(0.20,0.85,0.24,0.7,'' if 'cen' in var else 'Si1')
    #leg1.SetBorderSize(0)
    #leg1.SetFillStyle(0)
    #leg1.SetTextFont(42)
    #leg1.SetTextSize(0.035)
    #leg2=ROOT.TLegend(0.25,0.85,0.40,0.7, '' if 'cen' in var else 'Si2')
    #leg2.SetBorderSize(0)
    #leg2.SetFillStyle(0)
    #leg2.SetTextFont(42)
    #leg2.SetTextSize(0.035)
    mg=ROOT.TMultiGraph()


    for key in grColl:
        #if key[1]=='1':
        #    leg1.AddEntry(grColl[key],'','p')
        #else:
        #    beamName='e^{-} (50 GeV)'
        #    if key[0]==13: beamName='#mu^{-}(150 GeV)'
        #    leg2.AddEntry(grColl[key],beamName,'p')
        mg.Add(grColl[key])

    mg.Draw('ap')
    mg.GetXaxis().SetTitle('Si depletion thickness [#mum]')
    mg.GetYaxis().SetTitle(title)
    #mg.GetYaxis().SetRangeUser(0.5*grColl[key].GetYaxis().GetXmin(),grColl[key].GetYaxis().GetXmax()*1.5)
    
    #if not 'cen' in var : leg1.Draw()
    #leg2.Draw()
    txt=ROOT.TLatex()
    txt.SetNDC()
    txt.SetTextFont(42)
    txt.SetTextSize(0.035)
    txt.DrawLatex(0.2,0.88,"#splitline{#bf{CMS} #it{work in progress}}{#scale[0.8]{#it{Fast-timing testbeam}}}")

    if var=='adc2mip':
        fitfunc=ROOT.TF1('pol1','[0]*x',0,500)
        mg.Fit(fitfunc,'','0')

        hint = ROOT.TH1D("f95",'',100,100,350)
        ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(hint)
        hint.SetStats(False)
        hint.SetFillColor(ROOT.kCyan-8)
        hint.SetLineColor(ROOT.kBlue+3)
        hint.SetMarkerColor(ROOT.kBlue+3)
        hint.SetFillStyle(3001)
        hint.Draw("e3 same")
        txt.DrawLatex(0.6,0.89,'#splitline{#scale[0.8]{MIP/Si depletion thickness}}{(%3.3f#pm%3.3f) ADC/#mum}'%
                      (fitfunc.GetParameter(0),fitfunc.GetParError(0)))
    
    c.SaveAs(outDir+'/'+var+'_si%d_fit%d_chargeEst%d.png'%(siPadSel,fitType,chargeEst))


def main(argv=sys.argv):
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)
    
    baseDir=argv[1]
    dir_list = next(os.walk(baseDir))[1]
    chain=ROOT.TChain('mipcalib')
    for irun in dir_list:
        chain.Add(os.path.join(baseDir,irun+'/fitsummary.root'))
    print chain.GetEntries()


    for title,var,varUnc in [ ('MIP [ADC]','adc2mip','adc2mipUnc'),
                              ('#xi/MIP','relWidth','relWidthUnc'),
                              ('S(1)/N','mip1frac','mip1fracUnc'),
                              ]:
        for chargeEst in [0,1]:
            for siPadSel in [-1,2,3]:
                for fitType in [1,2,3,4]:
                    summarizeResultsFor(chain,chargeEst,siPadSel,fitType,title,var,varUnc,baseDir)

if __name__ == "__main__":
    main()
