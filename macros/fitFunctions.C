#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "Math/Math.h"
#include "Math/PdfFuncMathCore.h"

Double_t lanconvgau(Double_t *x, Double_t *par) {

  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location

  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;

  Double_t lanWidth=par[0];


  // MP shift correction
  mpc = par[1] - mpshift * lanWidth;

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,lanWidth) / lanWidth;
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,lanWidth) / lanWidth;
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}

Double_t sigFunc(Double_t *x, Double_t *par) 
{
  float noise=par[0]*TMath::Gaus(x[0],par[1],par[2]);

  par[6]=par[2];
  float mip1=lanconvgau(x,&par[3]);

  //the second landau is assumed to be a sum of 3 random variables Landau distributed as the first
  Float_t relSigma[6]={0.05,0.1,0.15,0.2,0.25,0.3};
  Float_t mpv2expAt3[6]={1.0537841320521606,1.1064830178068685,1.165639071227548,1.2134364957354753,1.2635444715793551,1.3174923805350};
  Float_t sigma2expAt3[6]={1.705351871064611,1.70993653332279,1.7061919950914877,1.687724097141382,1.7021148462819709,1.66699100227524};
  Float_t tryRelSigma(par[4]>0 ? par[3]/par[4] :0.051);
  if(tryRelSigma>0.3) tryRelSigma=0.299;
  par[8]=mpv2expAt3[0]*3*par[4];
  par[7]=sigma2expAt3[0]*TMath::Sqrt(3.0)*par[3];
  for(size_t i=0; i<5; i++)
    {
      if(tryRelSigma<relSigma[i] || tryRelSigma>=relSigma[i+1]) continue;
      par[8]=mpv2expAt3[i]*3*par[4];
      par[7]=sigma2expAt3[i]*TMath::Sqrt(3.0)*par[3];
    }
  
  par[10]=par[2];
  float mip2=lanconvgau(x,&par[7]);  
  return noise+mip1+mip2;
}



Double_t sigFuncCB(Double_t *x, Double_t *par) 
{
  float noise=par[0]*ROOT::Math::crystalball_pdf(x[0],par[3], par[4], par[2], par[1]);

  par[8]=par[2];
  float mip1=lanconvgau(x,&par[5]);

  //the second landau is assumed to be a sum of 3 random variables Landau distributed as the first
  Float_t relSigma[6]     = {0.05,0.1,0.15,0.2,0.25,0.3};
  Float_t mpv2expAt3[6]   = {1.05378, 1.10648, 1.1656, 1.2134, 1.26354, 1.3175};
  Float_t sigma2expAt3[6] = {1.70535, 1.70993, 1.7061, 1.6877, 1.70211, 1.6669};
  Float_t tryRelSigma(par[5]>0 ? par[6]/par[5] : 0.051);
  if(tryRelSigma>0.3) tryRelSigma=0.299;
  par[10]=mpv2expAt3[0]*3*par[6];
  par[9]=sigma2expAt3[0]*TMath::Sqrt(3.0)*par[5];
  for(size_t i=0; i<5; i++)
    {
      if(tryRelSigma<relSigma[i] || tryRelSigma>=relSigma[i+1]) continue;
      par[10]=mpv2expAt3[i]*3*par[6];
      par[9]=sigma2expAt3[i]*TMath::Sqrt(3.0)*par[5];
    }

  par[12]=par[2];
  float mip2=lanconvgau(x,&par[9]);  

  //sum up all contributions
  return noise+mip1+mip2;
}


