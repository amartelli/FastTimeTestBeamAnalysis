import ROOT
import os
import sys
import optparse
import pickle
from scipy.stats.mstats import mquantiles

"""
Fits a Landau to the spectra to get MIP calibration
"""
def runMIPcalibration(fileMap):
    
    print 'Starting MIP calibration'

    #loop over the trees to fill histos
    mipCalibHistos={}
    for key in fileMap:
        if not 'mu-' in fileMap[key] : continue
        print '\t analyzing ',fileMap[key]['mu-']
        siWidth,_ = key
        fIn=ROOT.TFile.Open(fileMap[key]['mu-'])
        H4treeSim=fIn.Get('H4treeReco')
        for i in xrange(0,H4treeSim.GetEntriesFast()):
            H4treeSim.GetEntry(i)
            for ich in xrange(0,H4treeSim.maxch):
                chNb=H4treeSim.ch[ich]                
                calibKey=(siWidth,chNb)

                #no point in looking to stuff below 5 keV
                if H4treeSim.charge_integ[ich]<5 : continue

                if not calibKey in mipCalibHistos: 
                    mipCalibHistos[calibKey]=ROOT.TH1F('mipcalib_%d_%d'%(siWidth,chNb),';Energy [keV];Events;',250,0,250)
                    mipCalibHistos[calibKey].SetDirectory(0)
                    mipCalibHistos[calibKey].Sumw2()
                mipCalibHistos[calibKey].Fill(H4treeSim.charge_integ[ich])

    #MIP calibration from Landau fit around max
    mipCalib={}
    mipEvol={}
    for calibKey in mipCalibHistos:
        mipCalibHistos[calibKey].Draw()
        mean=mipCalibHistos[calibKey].GetMean()
        rms=mipCalibHistos[calibKey].GetRMS()
        mipCalibHistos[calibKey].Fit('landau','LMRQ+','',mean-0.7*rms,mean+0.4*rms)

        lan=mipCalibHistos[calibKey].GetFunction('landau')

        if lan: 
            mipCalib[calibKey]=(lan.GetParameter(1),lan.GetParError(1))
        
            siWidth,chNb=calibKey[0],calibKey[1]
    
            if not chNb in mipEvol:
                mipEvol[chNb]=ROOT.TGraphErrors()
                mipEvol[chNb].SetName('mipevol_%d' % chNb)
                mipEvol[chNb].SetTitle('ch = %d' % chNb)
                mipEvol[chNb].SetMarkerStyle(20+len(mipEvol))
                mipEvol[chNb].SetFillStyle(0)
            np=mipEvol[chNb].GetN()
            mipEvol[chNb].SetPoint(np,siWidth,mipCalib[calibKey][0])
            mipEvol[chNb].SetPointError(np,0,mipCalib[calibKey][1])

    #show MIP evolution
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.12)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.05)
    drawOpt='ap'
    for chNb in mipEvol:
        mipEvol[chNb].Draw(drawOpt)
        mipEvol[chNb].GetXaxis().SetTitle('Si width [#mum]')
        mipEvol[chNb].GetYaxis().SetTitle('MIP/keV')
        drawOpt='p'
    leg=c.BuildLegend(0.2,0.95,0.6,0.7,'#bf{Geant4} #it{simulation}')
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    c.Modified()
    c.Update()
    c.SaveAs('mipevol.png')

    return mipCalib,mipCalibHistos

"""
electron spectra analysis
"""
def runElectronAnalysis(fileMap,mipCalib):
    
    #loop over the trees to fill histos
    enMoments={}
    for key in fileMap:
        if not 'e-' in fileMap[key] : continue
        siWidth,x0 = key

        edeps,edeps_full={},{}
        fIn=ROOT.TFile.Open(fileMap[key]['e-'])
        print fileMap[key]['e-']
        H4treeSim=fIn.Get('H4treeReco')
        for i in xrange(0,H4treeSim.GetEntriesFast()):
            H4treeSim.GetEntry(i)
            for ich in xrange(0,H4treeSim.maxch):
                chNb=H4treeSim.ch[ich] 

                calibKey=(siWidth,chNb)
                mip=1.0
                if calibKey in mipCalib: mip=mipCalib[calibKey][0]

                charge_integ      = H4treeSim.charge_integ[ich]/mip
                if charge_integ>0: 
                    if not chNb in edeps: edeps[chNb]=[]
                    edeps[chNb].append(charge_integ)

                charge_integ_full = H4treeSim.charge_integ_full[ich]/mip
                if charge_integ_full>0: 
                    if not chNb in edeps_full: edeps_full[chNb]=[]
                    edeps_full[chNb].append(charge_integ_full)

        for ch in edeps:
          momentKey=(siWidth,x0,ch)
          quantiles=mquantiles(edeps[ch],prob=[0.5,0.9])
          quantiles_full=mquantiles(edeps_full[ch],prob=[0.5,0.9])
          enMoments[momentKey]=(len(edeps[ch]),
                                quantiles[0],
                                quantiles[1],
                                len(edeps_full[ch]),
                                quantiles_full[0],
                                quantiles_full[1])
    return enMoments


"""
do final plots from summary
"""
def showSummaryPlots(url):
    cachefile = open(url,'r')
    enMoments=pickle.load(cachefile)        
    cachefile.close()

    #prepare canvas
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.12)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.05)
    
    #loop over plots
    chMarker={6:20,8:21}
    widthColor={120:1,200:ROOT.kGray,320:ROOT.kCyan-2}
    for idx,name,title in [(0,'fracinwindow','Fraction of events with energy'),
                           (1,'e50','Median energy'),
                           (2,'e90','Energy 90% quantile'),
                           (4,'e50full','Median energy (full pad)'),
                           (5,'e90full','Energy 90% quantile (full pad)')]:

        evolGr={}
        for key in enMoments:
            siWidth,x0,ch=key

            grKey=(siWidth,ch)
            if not grKey in evolGr:
                evolGr[grKey]=ROOT.TGraph()
                evolGr[grKey].SetName('%s_%d_%d' % (name,siWidth,ch) )
                evolGr[grKey].SetTitle('Si width=%d#mum (ch=%d)' % (siWidth,ch) )
                evolGr[grKey].SetFillStyle(0)
                evolGr[grKey].SetMarkerStyle(chMarker[ch])
                evolGr[grKey].SetLineColor(widthColor[siWidth])
                evolGr[grKey].SetMarkerColor(widthColor[siWidth])
            np= evolGr[grKey].GetN()
            val=enMoments[key][idx]
            if name=='fracinwindow' : val=(100.*val)/float(enMoments[key][3])
            evolGr[grKey].SetPoint(np,x0,val)
        
        #show evolution
        drawOpt='ap'
        for grKey in evolGr:
            evolGr[grKey].Draw(drawOpt)
            evolGr[grKey].GetXaxis().SetTitle('Pb X_{0}')
            if name=='fracinwindow':
                evolGr[grKey].GetYaxis().SetTitle('Fraction of events (%)')
            else:
                evolGr[grKey].GetYaxis().SetTitle('E/MIP')
            if drawOpt=='ap':
                evolGr[grKey].GetYaxis().SetRangeUser(0,2*evolGr[grKey].GetYaxis().GetXmax())
            drawOpt='p'
        leg=c.BuildLegend(0.2,0.9,0.6,0.6,'#bf{Geant4} #it{simulation} #scale[0.8]{: %s}' % title)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.03)
        c.Modified()
        c.Update()
        c.SaveAs('%s.png' % name)
        #raw_input()


"""
steer the script
"""
def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    #ROOT.gROOT.SetBatch(True)

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inDir',   dest='inDir',  help='input directory with files',            default=None, type='string')
    parser.add_option(      '--calib',   dest='calib',  help='pickle file with MIP calibration',      default=None, type='string')
    parser.add_option(      '--show',    dest='show',  help='pickle file with final results to show', default=None, type='string')
    (opt, args) = parser.parse_args()

    #map the files to analyze
    fileMap={}
    allFiles = [ os.path.join(opt.inDir,f) for f in os.listdir(opt.inDir) if os.path.isfile(os.path.join(opt.inDir,f)) ]
    for f in allFiles:
        tkns=f.split('_')
        version=int( tkns[ tkns.index('version')+1 ] )
        particle= tkns[ tkns.index('TBSim')+1 ]
        siWidth=version/10
        x0=version-siWidth*10

        #see FTTB215Conditions.py/DetectorConstruction.cc for the depleted thickness @ 600V
        trueSiWidth=140.5
        if siWidth==200: trueSiWidth=215.75
        if siWidth==320 : trueSiWidth=288.5

        key=(trueSiWidth,x0)
        if not key in fileMap: fileMap[key]={}
        fileMap[key][particle]=f

    #run first MIP calibration if needed
    if not opt.calib and not opt.show:
        mipCalib,mipCalibPlots=runMIPcalibration(fileMap)
        opt.calib = '.fasttimetbmipcalib.pck'
        cachefile = open(opt.calib,'w')
        pickle.dump(mipCalib, cachefile,pickle.HIGHEST_PROTOCOL)
        pickle.dump(mipCalibPlots, cachefile, pickle.HIGHEST_PROTOCOL)
        cachefile.close()
    
    #run electron analysis
    if not opt.show:

        #get MIP calibration
        cachefile = open(opt.calib,'r')
        mipCalib=pickle.load(cachefile)
        mipCalibPlots=pickle.load(cachefile)
        cachefile.close()
        
        enMoments=runElectronAnalysis(fileMap,mipCalib)
        opt.show = '.fasttimetbenmoments.pck'
        cachefile = open(opt.show,'w')
        pickle.dump(enMoments, cachefile,pickle.HIGHEST_PROTOCOL)
        cachefile.close()
        

    #show final summary
    showSummaryPlots(opt.show)


    
"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
