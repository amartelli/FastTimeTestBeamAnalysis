#include "interface/WaveformFit.hpp"
#include "Math/Interpolator.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TMath.h"

#include "assert.h"

#define MAX_INTERPOLATOR_POINTS 10000

namespace WaveformFit
{
  TProfile* refWave;
  TProfile* fitWave;
  Waveform* waveData;
  
  //   ROOT::Math::Minimizer* min =0;

  ROOT::Math::Interpolator* y_interpolator=0;
  ROOT::Math::Interpolator* y_err_interpolator=0;

  //   int nPointsUsed;
  //   int nSamplesToInterpolate;

  int xMin,xMax;
  float sampleRMS;
  
  double chi2(const double *par )
  {
    double chisq = 0;
    double delta = 0;
    for (int i=xMin;i<=xMax;++i)
      {
	//	delta=(refWave->GetBinContent(i)-(y_interpolator->Eval(refWave->GetBinCenter(i)+par[1])+par[0]))/(TMath::Sqrt(refWave->GetBinError(i)*refWave->GetBinError(i)+y_err_interpolator->Eval(refWave->GetBinCenter(i)+par[1])*y_err_interpolator->Eval(refWave->GetBinCenter(i)+par[1])));
	delta= (refWave->GetBinContent(i)-(y_interpolator->Eval(refWave->GetBinCenter(i)+par[1])+par[0])); //using no error for the moment in the Chi2 definition
	chisq += delta*delta;
      }
    return chisq;
  }

  double chi2Wave(const double *par )
  {
    double chisq = 0;
    double delta = 0;
    for (int i=xMin;i<=xMax;++i)
      {
	//	delta=(refWave->GetBinContent(i)-(y_interpolator->Eval(refWave->GetBinCenter(i)+par[1])+par[0]))/(TMath::Sqrt(refWave->GetBinError(i)*refWave->GetBinError(i)+y_err_interpolator->Eval(refWave->GetBinCenter(i)+par[1])*y_err_interpolator->Eval(refWave->GetBinCenter(i)+par[1])));
	if ( (*waveData)._times[i]*1.e9-par[1]<=(fitWave->GetXaxis()->GetXmin()+0.25) || (*waveData)._times[i]*1.e9-par[1]>=(fitWave->GetXaxis()->GetXmax()+0.25))
	  {
	    std::cout << "OUT OF BOUNDS " << (*waveData)._times[i]*1.e9-par[1] << "," << fitWave->GetXaxis()->GetXmin() << "," << fitWave->GetXaxis()->GetXmax() << std::endl;
	    chisq += 9999;
	  }
	else
	  {
	    // std::cout  << i << " " << (*waveData)._samples[i] << "," << (*waveData)._times[i] << " " << y_interpolator->Eval((*waveData)._times[i]*1.e9-par[1]) << std::endl;
	    delta= ((*waveData)._samples[i]-(y_interpolator->Eval((*waveData)._times[i]*1.e9-par[1])*par[0]))/sampleRMS; //fit par[0]*ref_shape(t-par[1]) par[0]=amplitude, par[1]=DeltaT
	    chisq += delta*delta;
	  }
      }
    return chisq;
  }

  void residuals(const double *par )
  {
    std::cout << "++++++++++++++++++++++++++++++++++" << std::endl;
    for (int i=xMin;i<=xMax;++i)
      {
	float res = ((*waveData)._samples[i]-(y_interpolator->Eval((*waveData)._times[i]*1.e9-par[1])*par[0]))/sampleRMS; //fit par[0]*ref_shape(t-par[1]) par[0]=amplitude, par[1]=DeltaT
	std::cout << "===> " << (i-xMin) << " " << (*waveData)._samples[i] << " " << res << std::endl;
      }
  }

  void alignWaveform(TProfile* ref_profile, TProfile* fit_profile, ROOT::Math::Minimizer* &minimizer)
  {    
    xMin=1;
    xMax=990;

    refWave=ref_profile;
    fitWave=fit_profile;

    minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    minimizer->SetMaxFunctionCalls(1000000);
    minimizer->SetMaxIterations(100000);
    minimizer->SetTolerance(1e-6);
    minimizer->SetPrintLevel(2);

    if (!y_interpolator)
      y_interpolator=new ROOT::Math::Interpolator(MAX_INTERPOLATOR_POINTS, ROOT::Math::Interpolation::kCSPLINE);
    if (!y_err_interpolator)
      y_err_interpolator=new ROOT::Math::Interpolator(MAX_INTERPOLATOR_POINTS, ROOT::Math::Interpolation::kCSPLINE);

    std::vector<double> x,y,y_err;

    for (int i=1;i<=fit_profile->GetNbinsX();++i)
      {
	x.push_back(fit_profile->GetBinCenter(i));
	y.push_back(fit_profile->GetBinContent(i));
	y_err.push_back(fit_profile->GetBinError(i));
      }

    y_interpolator->SetData(x,y);
    y_err_interpolator->SetData(x,y_err);

    ROOT::Math::Functor f(&chi2,2);
    minimizer->SetFunction(f);

    minimizer->SetVariable(0,"deltaV",0.,1e-2);
    minimizer->SetVariable(1,"deltaT",0.,1e-1);
    minimizer->Minimize();
   
    const double* par=minimizer->X();
    std::cout << "+++++ FIT RESULT: " << par[0] << "," << par[1] << std::endl;
    delete y_interpolator;
    delete y_err_interpolator;
    y_interpolator=0;
    y_err_interpolator=0;
    //     minimizer=min;
  }

  void fitWaveform(Waveform* wave, TProfile* amplitudeProfile, int x1, int x2, const Waveform::max_amplitude_informations& max, const Waveform::baseline_informations& waveRms, ROOT::Math::Minimizer* &minimizer)
  {
    xMin=x1;
    xMax=x2;

    waveData=wave;
    fitWave=amplitudeProfile;

    sampleRMS=waveRms.rms;

    minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    minimizer->SetMaxFunctionCalls(100000);
    minimizer->SetMaxIterations(100);
    minimizer->SetTolerance(1e-3);
    minimizer->SetPrintLevel(0);

    if (!y_interpolator)
      y_interpolator=new ROOT::Math::Interpolator(MAX_INTERPOLATOR_POINTS, ROOT::Math::Interpolation::kCSPLINE);

    std::vector<double> x,y;

    for (int i=1;i<=amplitudeProfile->GetNbinsX();++i)
      {
	x.push_back(amplitudeProfile->GetBinCenter(i));
	y.push_back(amplitudeProfile->GetBinContent(i));
      }

    y_interpolator->SetData(x,y);

    ROOT::Math::Functor f(&chi2Wave,2);
    minimizer->SetFunction(f);

    //    std::cout << "INIITIAL VALUES " << max.max_amplitude << "," << max.time_at_max*1.e9 << "," << max.sample_at_max << std::endl;

    //    minimizer->SetLimitedVariable(0,"amplitude",max.max_amplitude,1e-2,std::max(0.,(double)max.max_amplitude-500),std::min(4000.,(double)max.max_amplitude+500.));
    minimizer->SetLimitedVariable(0,"amplitude",max.max_amplitude,1e-2,-500., 4000.);
    //    minimizer->SetLimitedVariable(1,"deltaT",max.time_at_max*1.e9,1e-3,max.time_at_max*1.e9-0.5,max.time_at_max*1.e9+0.5);
    minimizer->SetFixedVariable(1,"deltaT",0.);

    minimizer->Minimize();

//     if (minimizer->Status()==0)
//       {
// 	const double* par=minimizer->X();
	
// 	// std::cout << "+++++ FIT RESULT: " << par[0] << "," << par[1] << std::endl;
// 	residuals(par);
//       }
    
    delete y_interpolator;
    y_interpolator=0;
  }

  
};

