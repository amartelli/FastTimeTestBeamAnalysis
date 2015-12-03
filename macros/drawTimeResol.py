#!/user/bin/env python
import array
import ROOT

hdt=ROOT.TH1F('hdt',';#Delta t=t_{2}-t_{1} [ns];Events/ (0.01 ns)',1000,-5,5)
hdt.Sumw2()
hdt.SetMarkerStyle(20)
hdt.GetYaxis().SetTitleOffset(1.2)

probSum = array.array('d', [0.5])
q = array.array('d', [0.0])

#open the file
fIn=ROOT.TFile.Open('timetest.root')
ttree=fIn.Get('ttree')

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c=ROOT.TCanvas('c','c',500,500)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.05)
c.SetTopMargin(0.05)
c.SetBottomMargin(0.12)

fitFunc=ROOT.TF1('gaus','gaus',-0.4,0.4)
resolGr=[]
offsetGr=[]
for t_est,t_est_label in [('t_max',        't_{peak}'),
                          ('t_max_frac30',    't_{30% peak}'),
                          ('t_max_frac50','t_{50% peak}')
                          ]:

    resolGr.append( ROOT.TGraphErrors() )
    resolGr[-1].SetMarkerStyle(20+2*len(resolGr))
    resolGr[-1].SetTitle(t_est_label)
    offsetGr.append( resolGr[-1].Clone() )

    for enmin,enmax in [(0,5),(5,10),(10,20),(20,50),(50,100)]:

        hdt.Reset('ICE')
        hdt.SetDirectory(fIn)
        ttree.Draw('%s[1]-%s[0]>>hdt'% (t_est,t_est),
                   'calibEn[0][0]>=%f && calibEn[0][0]<%f && calibEn[1][0]>=%f && calibEn[1][0]<%f' % (enmin,enmax,enmin,enmax),
                   'goff')
        hdt.GetQuantiles(1,q,probSum)
        t_off=q[0]
        t_offunc=1.253*hdt.GetMeanError()

        hdt.Reset('ICE')
        ttree.Draw('%s[1]-%s[0]-%f>>hdt'% (t_est,t_est,t_off),
                   'calibEn[0][0]>=%f && calibEn[0][0]<%f && calibEn[1][0]>=%f && calibEn[1][0]<%f' % (enmin,enmax,enmin,enmax),
                   'goff')
        hdt.GetQuantiles(1,q,probSum)

        meanEn=0.5*(enmax+enmin)
        h=hdt.Clone('h')
        h.Draw()
        h.Fit(fitFunc,'QR+','',-0.4/ROOT.TMath.Sqrt(meanEn),0.4/ROOT.TMath.Sqrt(meanEn))
        sigma,sigmaUnc=fitFunc.GetParameter(2),fitFunc.GetParError(2)
        
        h.GetYaxis().SetRangeUser(0,1.5*h.GetMaximum())
        h.GetXaxis().SetRangeUser(-3*sigma,3*sigma)
        fitFunc.SetRange(-3*sigma,3*sigma)
        fitFunc.Draw('same')
        np=resolGr[-1].GetN()
        resolGr[-1].SetPoint(np,meanEn,sigma*1e3)
        resolGr[-1].SetPointError(np,0.5*(enmax-enmin),sigmaUnc*1e3)
        offsetGr[-1].SetPoint(np,meanEn,t_off*1e3)
        offsetGr[-1].SetPointError(np,0.5*(enmax-enmin),t_offunc*1e3)

        txt=ROOT.TLatex()
        txt.SetTextFont(42)
        txt.SetTextSize(0.035)
        txt.SetNDC()
        txt.DrawLatex(0.15,0.9,'#bf{CMS}')
        txt.DrawLatex(0.15,0.85,'#scale[0.8]{#it{Fast-timing test beam preliminary}}')
        txt.DrawLatex(0.7,0.9,'#scale[0.8]{%3.0f < #frac{E}{MIP} < %3.0f}'%(enmin,enmax))
        txt.DrawLatex(0.7,0.85,'#scale[0.8]{%s estimator}'%t_est_label)
        txt.DrawLatex(0.7,0.8,'#scale[0.8]{#sigma_{#Deltat}=%3.0f#pm%3.0f ps}' % (1e3*sigma,1e3*sigmaUnc))
        txt.DrawLatex(0.7,0.75,'#scale[0.8]{offset=%3.0f#pm%3.0f ps}' % (t_off*1e3,t_offunc*1e3))
        c.Modified()
        c.Update()
        c.SaveAs('%s_%d.png' % (t_est,np))
        h.Delete()
    

resolFunc=ROOT.TF1('resolfunc','sqrt(pow([0],2)/x+pow([1],2))',0,500)
for grColl,name in [(offsetGr,'toffset'),
                    (resolGr,'tresol')]:
    c.Clear()    
    leg=ROOT.TLegend(0.65,0.88,0.9,0.5)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    txt=ROOT.TLatex()
    txt.SetTextFont(42)
    txt.SetTextSize(0.035)
    txt.SetNDC()
    
    for igr in xrange(0,len(grColl)):
        drawOpt='ap' if igr==0 else 'p'
        grColl[igr].Draw(drawOpt)
        grColl[igr].GetXaxis().SetTitle('Energy [MIP]')
        grColl[igr].GetYaxis().SetTitleOffset(1.5)
        
        if name=='tresol':
            c.SetLogy(True)
            grColl[igr].GetYaxis().SetTitle('#sigma_{#Deltat} [ps]')
            grColl[igr].GetYaxis().SetRangeUser(10,1000)
            grColl[igr].Fit(resolFunc,'Q0')
            if igr==0: resolFunc.Draw('same')
            fitresult_txt='#scale[0.8]{#sigma_{#Deltat}=#frac{%3.0f}{E/MIP}#oplus%3.0f}'%(resolFunc.GetParameter(0),resolFunc.GetParameter(1))    
            leg.AddEntry(resolGr[igr],
                         '#splitline{%s}{%s}'%(resolGr[igr].GetTitle(),fitresult_txt),
                         'p')
        else:
            c.SetLogy(False)
            grColl[igr].GetYaxis().SetTitle('#Deltat bias [ps]')
            grColl[igr].GetYaxis().SetRangeUser(100,300)
            leg.AddEntry(resolGr[igr],
                         resolGr[igr].GetTitle(),
                         'p')

    leg.Draw()
    txt.DrawLatex(0.15,0.9,'#bf{CMS}')
    txt.DrawLatex(0.15,0.85,'#scale[0.8]{#it{Fast-timing test beam preliminary}}')
    c.Modified()
    c.Update()
    c.SaveAs('%s.png'%name)
    raw_input()


