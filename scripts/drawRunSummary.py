# /usr/bin/python

import sys
import os
import ROOT

from FTTB215Conditions import CONDITIONS,DIODES

COLORS ={11:1,13:ROOT.kGray+1}
MARKERS={(11,''):22,(11,'1'):22,(11,'2'):26,(13,''):20,(13,'1'):20,(13,'2'):24}

def summarizeFor(ntuple,title,var,varUnc=None,outDir='./'):

    grColl={}
    subDets=['1','2']
    if 'cen' in var : subDets=['']
    for i in xrange(0,ntuple.GetEntriesFast()):

        ntuple.GetEntry(i)

        runConds = [x for x in CONDITIONS if x[2]==ntuple.run][0]
        beam     = runConds[0]
        beamEn   = runConds[1]
        diode    = runConds[5]
        siwidth  = DIODES[diode][3]

        for sipad in subDets:

            varVal=getattr(ntuple,var+sipad)
            varValUnc=0
            try:
                varValUnc=getattr(ntuple,varUnc+sipad)
            except:
                pass

            if 'mip_Si' in var:
                chi2ndof=getattr(ntuple,'chi2ndof_Si'+sipad)
                if chi2ndof>2 or chi2ndof<0.5:
                    continue
                if varValUnc/varVal < 0.01 :
                    continue

            key=(beam,sipad)
            if not key in grColl:
                grColl[key]=ROOT.TGraphErrors();
                grColl[key].SetName('%s_%d_%s' % (var,beam,sipad) )
                grColl[key].SetMarkerStyle( MARKERS[(beam,sipad)] )
                grColl[key].SetMarkerColor( COLORS[beam] )
                grColl[key].SetLineColor( COLORS[beam] )
                grColl[key].SetLineWidth(2)
                grColl[key].SetFillStyle(0)
                
            np=grColl[key].GetN()
            grColl[key].SetPoint(np,siwidth,varVal)    
            grColl[key].SetPointError(np,0,varValUnc)    
  
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetRightMargin(0.02)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.12)
    drawOpt='ap'
    leg1=ROOT.TLegend(0.20,0.85,0.24,0.7,'' if 'cen' in var else 'Si1')
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(42)
    leg1.SetTextSize(0.035)
    leg2=ROOT.TLegend(0.25,0.85,0.40,0.7, '' if 'cen' in var else 'Si2')
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.SetTextFont(42)
    leg2.SetTextSize(0.035)
    mg=ROOT.TMultiGraph()

    toPlot=[ (11,'1'), (11,'2'),(13,'1'),(13,'2')]
    if 'cen' in var : toPlot=[(11,''),(13,'')]
    for key in toPlot: 
        if key[1]=='1':
            leg1.AddEntry(grColl[key],'','p')
        else:
            beamName='e^{-} (50 GeV)'
            if key[0]==13: beamName='#mu^{-}(150 GeV)'
            leg2.AddEntry(grColl[key],beamName,'p')
        mg.Add(grColl[key])

    mg.Draw('ap')
    mg.GetXaxis().SetTitle('Si depletion thickness [#mum]')
    mg.GetYaxis().SetTitle(title)
    #mg.GetYaxis().SetRangeUser(0.5*grColl[key].GetYaxis().GetXmin(),grColl[key].GetYaxis().GetXmax()*1.5)
    
    if not 'cen' in var : leg1.Draw()
    leg2.Draw()
    txt=ROOT.TLatex()
    txt.SetNDC()
    txt.SetTextFont(42)
    txt.SetTextSize(0.035)
    txt.DrawLatex(0.2,0.88,"#splitline{#bf{CMS} #it{work in progress}}{#scale[0.8]{#it{Fast-timing testbeam}}}")

    if var=='mip_Si' :
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
    
    c.SaveAs(outDir+'/'+var+'.png')


def main(argv=sys.argv):
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)
    
    url=argv[1]
    baseDir=os.path.dirname(url)
    fIn=ROOT.TFile.Open(url)
    ntuple=fIn.Get('summary')
    
    for title,var,varUnc in [ ('#chi^2/ndf','chi2ndof_Si', None),
                              ('MIP [ADC]','mip_Si','mipUnc_Si'),
                              ('<Noise> [ADC]','noise_Si','noiseUnc_Si'),
                              ('#sigma_{noise} [ADC]','noiseSigma_Si','noiseSigmaUnc_Si'),
                              ('<x> [mm]','xcen','xcenUnc'),
                              ('<y> [mm]','ycen','ycenUnc'),
                              ]:
        summarizeFor(ntuple,title,var,varUnc,baseDir)

if __name__ == "__main__":
    main()
