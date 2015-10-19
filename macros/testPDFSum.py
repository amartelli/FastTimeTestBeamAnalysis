import ROOT

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptFit(1111)

c=ROOT.TCanvas('c','c',500,500)

dist,fitDist='Landau','landau'
#dist,fitDist='Gaus','gaus'

mpv=50
sigma=0.2*mpv
h=ROOT.TH1F('totalcharge','totalcharge',400,0,2000)
grMPVLin,grSigmaLin=ROOT.TGraph(),ROOT.TGraph()
grMPVFit,grSigmaFit=ROOT.TGraphErrors(),ROOT.TGraphErrors()

for nmips in xrange(1,11):
    h.Reset('ICE')
    for i in xrange(0,10000):
        charge=0
        for k in xrange(0,nmips):
            charge+=getattr(ROOT.gRandom,dist)(mpv,sigma)
        h.Fill(charge)


    c.Clear()
    h.Fit(fitDist)
    np=grMPVLin.GetN()
    grMPVLin.SetPoint(np,nmips,1)
    grSigmaLin.SetPoint(np,nmips,1)
    grMPVFit.SetPoint(np,nmips,(h.GetFunction(fitDist).GetParameter(1))/(nmips*mpv))
    grMPVFit.SetPointError(np,0,(h.GetFunction(fitDist).GetParError(1))/(nmips*mpv))


    sigmaExp=ROOT.TMath.Sqrt( nmips*(sigma)**2 )
    grSigmaFit.SetPoint(np,nmips,h.GetFunction(fitDist).GetParameter(2)/sigmaExp)
    grSigmaFit.SetPointError(np,0,h.GetFunction(fitDist).GetParError(2)/sigmaExp)
    c.Modified()
    c.Update()
    #raw_input()


testCtr=0
for ytitle,ref,fit in [('MPV/lin. sum',grMPVLin,grMPVFit),
                       ('#Delta #xi/lin. sum',grSigmaLin,grSigmaFit)]:
    c.Clear()
    ref.Draw('al')
    ref.SetFillStyle(0)
    ref.SetLineColor(ROOT.kGray)
    fit.Draw('p')
    fit.SetFillStyle(0)
    fit.SetMarkerStyle(20)
    ref.GetYaxis().SetTitle(ytitle)
    ref.GetXaxis().SetTitle('# MIPs')

    c.Modified()
    c.Update()
    raw_input()
    c.SaveAs('tes%sSum_%d.png'%(dist,testCtr))
    testCtr=testCtr+1
