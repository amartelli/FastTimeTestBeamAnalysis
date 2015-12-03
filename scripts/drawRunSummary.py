# /usr/bin/python

import sys
import os
import ROOT
from FTTB215Conditions import CONDITIONS,DIODES

COLORS ={11:ROOT.kBlue+1,13:1}
MARKERS={(11,2):22,(11,3):26,(13,2):20,(13,3):24}

def summarizeResultsFor(ntuple,chargeEst,siPadSel,fitType,title,var,varUnc=None,outDir='./'):

    grColl,shiftGrColl={},{}
    siWidthsCounter={}
    for i in xrange(0,ntuple.GetEntries()):

        ntuple.GetEntry(i)

        if chargeEst!=ntuple.chargeEstimator : continue
        if siPadSel>0 and ntuple.siPad!=siPadSel:continue
        if ntuple.siPad==2 and ntuple.runNumber==3363: continue
        if fitType != ntuple.fitType : continue
        #if ntuple.prob<0.01 : continue
        #if ntuple.adc2mipUnc/ntuple.adc2mip>0.3: continue
        if ntuple.adc2mipUnc/ntuple.adc2mip<0.001 : continue

        #run conditions
        runTag=ntuple.runNumber
        if runTag==0 : runTag=3375 # muons might have been summed
        runConds = [x for x in CONDITIONS if x[2]==runTag][0]
        beam     = runConds[0]
        beamEn   = runConds[1]
        diode    = runConds[5]
        siwidth  = DIODES[diode][3]
        if not siwidth in siWidthsCounter: siWidthsCounter[siwidth]=0
        siWidthsCounter[siwidth]+=1

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
            shiftGrColl[beam]=grColl[beam].Clone()
                
        np=grColl[beam].GetN()
        grColl[beam].SetPoint(np,siwidth,varVal)    
        grColl[beam].SetPointError(np,0,varValUnc)    
        shiftGrColl[beam].SetPoint(np,siwidth+3*(siWidthsCounter[siwidth]-2)-4,varVal)    
        shiftGrColl[beam].SetPointError(np,0,varValUnc)    
  
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetRightMargin(0.02)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.12)
    drawOpt='ap'
  
    mg=ROOT.TMultiGraph()
    shiftmg=ROOT.TMultiGraph()

    for key in grColl:
        mg.Add(grColl[key])
        shiftmg.Add(shiftGrColl[key])
    
    txt=ROOT.TLatex()
    txt.SetNDC()
    txt.SetTextFont(42)
    txt.SetTextSize(0.035)

    toReturn=None
    if var=='adc2mip':
        fitfunc=ROOT.TF1('pol1','[0]*x',0,500)
        mg.Fit(fitfunc,'','0')
        toReturn=fitfunc

        hint = ROOT.TH1D("f95",'',100,120,320)
        ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(hint)
        hint.SetStats(False)
        hint.SetFillColor(ROOT.kCyan-8)
        hint.SetLineColor(ROOT.kBlue+3)
        hint.SetMarkerColor(ROOT.kBlue+3)
        hint.SetFillStyle(3001)
        hint.Draw("e3")
        fitfunc.Draw('same')
        hint.GetXaxis().SetTitle('Si depletion thickness [#mum]')
        hint.GetYaxis().SetTitle(title)
        shiftmg.Draw('p')
        txt.DrawLatex(0.6,0.89,'#splitline{#scale[0.8]{MIP/Si depletion thickness}}{(%3.3f#pm%3.3f) ADC/#mum}'%
                      (fitfunc.GetParameter(0),fitfunc.GetParError(0)))
    else:
        shiftmg.Draw('ap')
        shiftmg.GetXaxis().SetTitle('Si depletion thickness [#mum]')
        shiftmg.GetYaxis().SetTitle(title)
  
    txt.DrawLatex(0.2,0.88,"#splitline{#bf{CMS} #it{work in progress}}{#scale[0.8]{#it{Fast-timing testbeam}}}")
    
    for ext in ['png','pdf','C']:
        c.SaveAs(outDir+'/'+var+'_si%d_fit%d_chargeEst%d.%s'%(siPadSel,fitType,chargeEst,ext))

    return toReturn

"""
"""
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

    allResults={}
    for title,var,varUnc in [ ('MIP [ADC]','adc2mip','adc2mipUnc'),
                              ('#xi/MIP','relWidth','relWidthUnc'),
                              ('S(1)/N','mip1frac','mip1fracUnc'),
                              ]:
        for chargeEst in [0,1]:
            for siPadSel in [-1]: #,2,3]: #[-1,2,3]
                for fitType in [0,1,2]: #[1,2,3,4]
                    result=summarizeResultsFor(chain,chargeEst,siPadSel,fitType,title,var,varUnc,baseDir)
                    if result:
                        allResults[(var,chargeEst,siPadSel,fitType)]=result
    
    #save parameterisations to a file
    fOut=ROOT.TFile.Open( '%s/runsummary.root'%baseDir, 'RECREATE')
    for key in allResults:
        name='%s_chargeEst%d_siPad%d_fit%d' % (key[0], key[1], key[2], key[3])
        allResults[key].SetName(name)
        allResults[key].Write()
    fOut.Close()

if __name__ == "__main__":
    main()
