//#include "interface/Waveform.hpp"
#ifndef WaveForm_h
#define WaveForm_h

#include "TROOT.h"
#include "Math/Interpolator.h"

#include <vector> 
#include <algorithm> 
#include <iostream> 

#include "TGraph.h"
#include "TF1.h"
#include "TMath.h"

using namespace std;

#define WARNING_ERROR 0
#define DEB 0

class Waveform
{
 public:
  std::vector<double> _samples;
  std::vector<double> _times;
  ROOT::Math::Interpolator* _interpolator;
  
  struct max_amplitude_informations
  {
    float max_amplitude;
    float time_at_max;
    int sample_at_max;
  };

  
  struct baseline_informations
  {
    float pedestal;
    float rms;
    float slope;
  };

 Waveform() : _interpolator(0)
    {
    };
  
  //constructor from float array
 Waveform(int nSamples, float* times, float* samples) : _interpolator(0)
    {
      _samples.resize(nSamples);
      _times.resize(nSamples);
      for (int i(0);i<nSamples;++i)
	{
	  _samples[i]=samples[i];
	  _times[i]=times[i];
	}
    };

  //constructor from std::vector<float>
 Waveform(const std::vector<float>& times, const std::vector<float>& samples) : _interpolator(0)
    {
      _samples.resize(samples.size());
      _times.resize(samples.size());
      std::copy(samples.begin(), samples.end(), _samples.begin());
      std::copy(times.begin(), times.end(), _times.begin());
    };

  ~Waveform()
    {
      if (_interpolator)
	delete _interpolator;
    };

  //add samples
  void addSample(const float& sample, float sampleTimeSize=0.2e-9)
  {
    _samples.push_back(sample);
    if (_times.size()>0)
      _times.push_back(_times.back()+sampleTimeSize);
    else
      _times.push_back(0.);
  }
  
  void addTimeAndSample(const float& time, const float& sample)
  {
    _samples.push_back(sample);
    _times.push_back(time);
  }

  void interpolate()
  {
    if  (_interpolator)
      return;

    _interpolator=new ROOT::Math::Interpolator(_samples.size(), ROOT::Math::Interpolation::kCSPLINE);
    _interpolator->SetData(_times,_samples);
  }

  //open a window of nSamples centered keeping into account where is the maximum and some characteristics of the shape
  void find_interesting_samples(int nSamples, const max_amplitude_informations& maxInfos, float riseTime, float fallTime, int& x1,int &x2);

  //Get the max amplitude between x1 and x2 using nSamplesAroundMax samples and a parabolic fit
  max_amplitude_informations max_amplitude(const int& x1, const int& x2, int nSamplesAroundMax=5) const;

  //Get the time at a given fraction of the amplitude for times between x1 and x2 
  float interpolatedValue(const int& i1, int SampleToInterpolate=9) const;

  //Get the time at a given fraction of the amplitude for times between x1 and x2 
  float time_at_frac(const int& i1, const int& i2, const float& frac, const max_amplitude_informations& maxInfos, int SampleToInterpolate=5) const;

  //Get the time at a given fraction of the amplitude for times between x1 and x2 
  float time_at_frac(const float& t1, const float& t2, const float& frac, const max_amplitude_informations& maxInfos, int SampleToInterpolate=5) const;

  //Get the baseline (pedestal and RMS) informations computed between x1 and x2
  baseline_informations baseline(const int& x1, const int& x2, bool fitSlope=false) const;

  //get values at crossing a specif threshold */
  std::vector<float> time_at_threshold(const float& t1, const float& t2, const float& threshold, int SampleToInterpolate=5) const; 

  //get values at crossing a specif threshold */
  std::vector<float> time_at_threshold(const int& x1, const int& x2, const float& threshold, int SampleToInterpolate=5) const; 

  //get value of integrated amplitude between x1 and x2, subtracting pedestal
  float charge_integrated(const int& x1, const int& x2, float pedestal=0) const;
      
  float integral(const int& x1, const int& x2) const
  {
    float integral=0;
    for(int bin=x1; bin<x2; bin++)
      integral += _samples[bin];
    return integral;
  };

  //clear
  void clear()
  {
    _samples.clear();
    _times.clear();
    if (_interpolator)
      delete _interpolator;
  };

  //substract a fixed value from the samples
  void offset(const float& value, const float& slope=0.)
  {
    for (unsigned int i(0);i<_samples.size();++i)
      _samples[i]-=(value+slope*(_times[i]*1.E9)); //slope has to be expressed in ADC/ns
  };

  //rescale all the samples by a given rescale factor
  void rescale(const float& rescale_factor)
  {
    for (unsigned int i(0);i<_samples.size();++i)
      _samples[i]*=rescale_factor;
  };

  //shift all samples by a given amount in time */
  void shift_time(const float& time_offset)
  {
    for (unsigned int i(0);i<_times.size();++i)
      _times[i]+=time_offset;
  } 
    
};

#endif



float Waveform::charge_integrated(const int& x1, const int& x2, float pedestal) const
{

  float return_value = 0;


  if(DEB){
    std::cout << " Waveform::charge_integrated " << std::endl;
    std::cout << " x1 = " << x1 << " x2 = " << x2 << std::endl;
  }

  if (x1<0 || x2> int(_samples.size()-1))
    {
      if (WARNING_ERROR)
	std::cout << "WARNING::charge_integrated::gate is outside samples range" << std::endl;
      return -9999.;
    }
  
  for (int i(x1); i<=x2; i++){
    return_value+=_samples[i]-pedestal;
  }
  
  return return_value;

}

Waveform::max_amplitude_informations Waveform::max_amplitude(const int& x1, const int& x2, int nSamplesAroundMax) const
{
  max_amplitude_informations return_value;  

  if(DEB){
    std::cout << " Waveform::max_amplitude " << std::endl;
    std::cout << " x1 = " << x1 << " x2 = " << x2 << std::endl;
  }


  if (x1<0 || x2>int(_samples.size()-1))
    {
      if (WARNING_ERROR)
	std::cout << "WARNING::Waveform::max_amplitude::gate is outside samples range" << std::endl;
      return return_value;
    }

  int imax=-1;
  float max=-999.;

  for (int i(x1);i<=x2;++i)
    {
      //      std::cout << "++++ " <<  i << "," << _times[i] << "," << _samples[i] << "," << max << "," << imax << std::endl;
      if (_samples[i]>=max)
	{
	  max=_samples[i];
	  imax=i;
	}
    }


  if (imax>-1)
    {
      float x[nSamplesAroundMax];
      float y[nSamplesAroundMax];
      int nSamples=0;
      for (int i(0);i<nSamplesAroundMax;++i)
	{
	  if ( (imax-(nSamplesAroundMax-1)/2+i)>=x1 && (imax-(nSamplesAroundMax-1)/2+i)<=x2)
	    {
	      //	      x[i]=imax-(nSamplesAroundMax-1)/2+i;
	      x[i]=_times[imax-(nSamplesAroundMax-1)/2+i]*1e9;
	      y[i]=_samples[imax-(nSamplesAroundMax-1)/2+i];
	      //	      std::cout <<  imax << "," << imax-(nSamplesAroundMax-1)/2+i << "," << x[i] << "," << y[i] << std::endl;
      
	      ++nSamples;
	    }
	  else
	    {
	      if (WARNING_ERROR)
		std::cout << "WARNING::Waveform::max_amplitude::maximum found too close to gate edges. Increase gate width" << std::endl;
	    }
	}

      if (nSamples>3)
	{
	  //Now fit with parabolic function around maximum value
	  TGraph* graph=new TGraph(nSamples,x,y);
	  graph->Fit("pol2","Q0+");

	  //FIXME Add a check on the FIT status
	  double *par=graph->GetFunction("pol2")->GetParameters();
	  return_value.max_amplitude=par[0]-(par[1]*par[1]/(4*par[2]));
	  return_value.time_at_max=-(par[1]/(2*par[2]))/1.e9;
	  return_value.sample_at_max=imax;

	  delete graph;
	}
      else
	{
	  if (WARNING_ERROR)
	    std::cout << "WARNING::Waveform::max_amplitude::not enough samples to fit fot maximum. Returning unfitted position" << std::endl;
	  return_value.max_amplitude=max;
	  return_value.time_at_max=_times[imax];
	  return_value.sample_at_max=imax;
	}

    }
  return return_value;
};


//Get the residual value between a sample and expected value from interpolating samples
float Waveform::interpolatedValue(const int& i1, int SampleToInterpolate) const
{
  int cfSample=i1;

  //  std::cout <<  "++++ "  << cfSample << std::endl;
  if (cfSample>-1)
    {
      double x[SampleToInterpolate];
      double y[SampleToInterpolate];
      double syy=0,sxy=0,sxx=0,sx=0,sy=0;
      int nSamples=0;
      for (int i(0);i<SampleToInterpolate;++i)
	{
	  if ( (cfSample-(SampleToInterpolate-1)/2+i)>=0 && (cfSample-(SampleToInterpolate-1)/2+i)<=(int)_samples.size() )
	    {
	      //	      x[i]=cfSample-(SampleToInterpolate-1)/2+i;
	      x[i]=(_times[cfSample-(SampleToInterpolate-1)/2+i] - _times[cfSample])*1e9 ;
	      y[i]=_samples[cfSample-(SampleToInterpolate-1)/2+i];
	      //std::cout <<  cfSample << "," << cfSample-(SampleToInterpolate-1)/2+i << "," << x[i] << "," << y[i] << std::endl;
	      syy+=y[i]*y[i];
	      sxy+=x[i]*y[i];
	      sxx+=x[i]*x[i];
	      sy+=y[i];
	      sx+=x[i];
      	      ++nSamples;
	    }
	  else
	    {
	      if (WARNING_ERROR)
		std::cout << "WARNING::Waveform::max_amplitude::maximum found too close to gate edges. Increase gate width" << std::endl;
	    }
	}

      if (nSamples>1)
	{
// 	  //Now fit with parabolic function around maximum value
// 	  TGraph* graph=new TGraph(nSamples,x,y);
// 	  graph->Fit("pol1","Q0+");

// 	  //FIXME Add a check on the FIT status
// 	  double *par=graph->GetFunction("pol1")->GetParameters();

	  double b = (nSamples*sxy-sx*sy)/(double) (nSamples*sxx - sx*sx);
	  double a = (sy - b*sx)/(double) nSamples;
	  return a;
	  
// 	  return -(par[0]/par[1])/1.e9;
// 	  delete graph;
	}
      else
	{
	  if (WARNING_ERROR)
	    std::cout << "WARNING::Waveform::max_amplitude::not enough samples to interpolate. Returning -999." << std::endl;
	  return -999.;
	}
    }

  return -999.;
};


//Get the time at a given fraction of the amplitude for times between x1 and x2 
float Waveform::time_at_frac(const float& t1, const float& t2, const float& frac, const max_amplitude_informations& maxInfos, int SampleToInterpolate) const
{
  //std::cout << "_____ " << t1 << "," << t2 << std::endl;
  int tmin=0;
  int tmax=0;
  for (int i(0);i<(int)_times.size();++i)
    {
      //      std::cout << i << "," << _times[i] << std::endl;
      if (_times[i]<t1)
	tmin=i;
      if (_times[i]<t2)
	tmax=i;
    }

  return time_at_frac(tmin,tmax,frac,maxInfos,SampleToInterpolate);
}


//Get the time at a given fraction of the amplitude for times between x1 and x2 
float Waveform::time_at_frac(const int& x1, const int& x2, const float& frac, const max_amplitude_informations& maxInfos, int SampleToInterpolate) const
{
  int cfSample=maxInfos.sample_at_max;
  if(cfSample>= (int)_samples.size())
    {
      std::cout << "[Waveform::time_at_frac] Sample at max (" << cfSample << ") >= size of samples...bailing out" << std::endl;
      return -999.;
    }
  //  std::cout << "===== " << x1 << "," << x2 << std::endl;

  for(int iSample=(int)maxInfos.sample_at_max; iSample>max(x1,0); iSample--)
    {
      if(_samples[iSample] < maxInfos.max_amplitude*frac) 
	{
	  cfSample = iSample;
 	  break;
	}
    }
  
  //  std::cout <<  "++++ "  << cfSample << std::endl;
  if (cfSample>-1)
    {
      double x[SampleToInterpolate];
      double y[SampleToInterpolate];
      double syy=0,sxy=0,sxx=0,sx=0,sy=0;
      int nSamples=0;
      for (int i(0);i<SampleToInterpolate;++i)
	{
	  if ( (cfSample-(SampleToInterpolate-1)/2+i)>=max(x1,0) && (cfSample-(SampleToInterpolate-1)/2+i)<=min(x2,(int)_samples.size()))
	    {
	      //	      x[i]=cfSample-(SampleToInterpolate-1)/2+i;
	      x[i]=_times[cfSample-(SampleToInterpolate-1)/2+i]*1e9;
	      y[i]=_samples[cfSample-(SampleToInterpolate-1)/2+i]-maxInfos.max_amplitude*frac;
	      //std::cout <<  cfSample << "," << cfSample-(SampleToInterpolate-1)/2+i << "," << x[i] << "," << y[i] << std::endl;
	      syy+=y[i]*y[i];
	      sxy+=x[i]*y[i];
	      sxx+=x[i]*x[i];
	      sy+=y[i];
	      sx+=x[i];
      	      ++nSamples;
	    }
	  else
	    {
	      if (WARNING_ERROR)
		std::cout << "WARNING::Waveform::max_amplitude::maximum found too close to gate edges. Increase gate width" << std::endl;
	    }
	}

      if (nSamples>1)
	{
// 	  //Now fit with parabolic function around maximum value
// 	  TGraph* graph=new TGraph(nSamples,x,y);
// 	  graph->Fit("pol1","Q0+");

// 	  //FIXME Add a check on the FIT status
// 	  double *par=graph->GetFunction("pol1")->GetParameters();

	  double b = (nSamples*sxy-sx*sy)/(double) (nSamples*sxx - sx*sx);
	  double a = (sy - b*sx)/(double) nSamples;
	  return -(a/b)/1.e9;
	  
// 	  return -(par[0]/par[1])/1.e9;
// 	  delete graph;
	}
      else
	{
	  if (WARNING_ERROR)
	    std::cout << "WARNING::Waveform::max_amplitude::not enough samples to fit fot maximum. Returning unfitted position" << std::endl;
	  return _times[cfSample];
	}
    }

  return -999.;
};

//Get the crossing time of the amplitude at a given threshold for times between t1 and t2 
std::vector<float> Waveform::time_at_threshold(const float& t1, const float& t2, const float& threshold, int SampleToInterpolate) const 
{
  //std::cout << "_____ " << t1 << "," << t2 << std::endl;
  int tmin=0;
  int tmax=0;
  for (int i(0);i<(int)_times.size();++i)
    {
      //      std::cout << i << "," << _times[i] << std::endl;
      if (_times[i]<t1)
	tmin=i;
      if (_times[i]<t2)
	tmax=i;
    }

  return time_at_threshold(tmin,tmax,threshold,SampleToInterpolate);
}


//Get all the crossing times of the waveform at a given threshold for samples between x1 and x2 
std::vector<float> Waveform::time_at_threshold(const int& x1, const int& x2, const float& threshold, int SampleToInterpolate) const 
{
  int direction= 1 - 2*( _samples[max(x1,0)]<threshold ); //=-1 if first sample is below threshold, =1 if above
  
  std::vector<int> cfSample;
  
  for(int iSample=max(x1,0); iSample< min( (int) _samples.size() , x2); ++iSample)
    {
      if( (direction==-1 && (_samples[iSample] > threshold ) ) 
	  || (direction==1 && (_samples[iSample] < threshold ) ) ) 
	{
	  cfSample.push_back(iSample);
	  direction*=-1;
	}
    }

  std::vector<float> crossingTimes;

  for(int icross=0;icross<(int)cfSample.size();++icross)
    {
      double x[SampleToInterpolate];
      double y[SampleToInterpolate];
      double syy=0,sxy=0,sxx=0,sx=0,sy=0;

      int nSamples=0;
      for (int i(0);i<SampleToInterpolate;++i)
	{

	  if ( (cfSample[icross]-(SampleToInterpolate-1)/2+i)>=max(x1,0) && (cfSample[icross]-(SampleToInterpolate-1)/2+i)<=min(x2,(int)_samples.size()))
	    {
	      //	      x[i]=cfSample[icross]-(SampleToInterpolate-1)/2+i;
	      x[i]=_times[cfSample[icross]-(SampleToInterpolate-1)/2+i]*1e9;
	      y[i]=_samples[cfSample[icross]-(SampleToInterpolate-1)/2+i]-threshold;
	      //std::cout <<  cfSample[icross] << "," << cfSample[icross]-(SampleToInterpolate-1)/2+i << "," << x[i] << "," << y[i] << std::endl;
	      syy+=y[i]*y[i];
	      sxy+=x[i]*y[i];
	      sxx+=x[i]*x[i];
	      sy+=y[i];
	      sx+=x[i];
      	      ++nSamples;
	    }
	  else
	    {
	      if (WARNING_ERROR)
		std::cout << "WARNING::Waveform::max_amplitude::maximum found too close to gate edges. Increase gate width" << std::endl;
	    }
	}

      if (nSamples>1)
	{
// 	  //Now fit with parabolic function around maximum value
// 	  TGraph* graph=new TGraph(nSamples,x,y);
// 	  graph->Fit("pol1","Q0+");

// 	  //FIXME Add a check on the FIT status
// 	  double *par=graph->GetFunction("pol1")->GetParameters();

// 	  crossingTimes.push_back( -(par[0]/par[1])/1.e9 );
// 	  delete graph;
//        //Use the regression formula instead of the fit
	  double b = (nSamples*sxy-sx*sy)/(double) (nSamples*sxx - sx*sx);
	  double a = (sy - b*sx)/(double) nSamples;
	  crossingTimes.push_back( -(a/b)/1.e9 );
	}
      else
	{
	  if (WARNING_ERROR)
	    std::cout << "WARNING::Waveform::max_amplitude::not enough samples to fit fot maximum. Returning unfitted position" << std::endl;
	  crossingTimes.push_back( _times[cfSample[icross]] );
	}
      
    }
  return crossingTimes;
}

void Waveform::find_interesting_samples(int nSamples, const max_amplitude_informations& maxInfos, float riseTime, float fallTime, int& x1,int &x2)
{
  //std::cout << "_____ " << t1 << "," << t2 << std::endl;
  int tmin=0;
  int tmax=0;
  for (int i(0);i<(int)_times.size();++i)
    {
      //      std::cout << i << "," << _times[i] << std::endl;
      if (_times[i]<maxInfos.time_at_max-2.*riseTime)
	tmin=i;
      if (_times[i]<maxInfos.time_at_max+1.*fallTime)
	tmax=i;
    }

  if ((tmax-tmin+1)==nSamples)
    {
      x1=tmin;
      x2=tmax;
    }
  else if ((tmax-tmin+1)<nSamples)
    {
      if ((int)_samples.size()<nSamples)
	{
	  x1=0;
	  x2=(int)_samples.size()-1;
	}
      else if (tmin+nSamples-1<(int)_samples.size())
	{
	  x1=tmin;
	  x2=(tmin+nSamples-1);
	}
      else
	{
	  x2=(int)_samples.size();
	  x1=x2-nSamples+1;
	}
    }
  else if ((tmax-tmin+1)>nSamples)
    {
      x1=tmin;
      x2=(tmin+nSamples-1);
    }
}
//Get the baseline (pedestal and RMS) informations computed between x1 and x2
Waveform::baseline_informations Waveform::baseline(const int& x1, const int& x2, bool fitSlope) const
{

  baseline_informations return_value;

  if (x1<0 || x2>((int)_samples.size()-1))
    {
      if (WARNING_ERROR)
	std::cout << "WARNING::Waveform::baseline::gate is outside samples range" << std::endl;
      return return_value;
    }
  
  if ((x2-x1)<2)
    {
      if (WARNING_ERROR)
	std::cout << "WARNING::Waveform::baseline::you need >2 samples to get pedestal & rms" << std::endl;
      return return_value;
    }

  double mean=0;
  double rms=0;
  double x[x2-x1+1];
  double y[x2-x1+1];
  double syy=0,sxy=0,sxx=0,sx=0,sy=0;
	  
  for (int i(x1);i<=x2;++i)
    {
      mean+=_samples[i];
      rms+=_samples[i]*_samples[i];
      x[i-x1]=_times[i]*1e9;
      y[i-x1]=_samples[i];
      syy+=y[i-x1]*y[i-x1];
      sxy+=x[i-x1]*y[i-x1];
      sxx+=x[i-x1]*x[i-x1];
      sy+=y[i-x1];
      sx+=x[i-x1];
    }
  
  mean=mean/(double)(x2-x1+1);
  rms=TMath::Sqrt((x2-x1+1)/(x2-x1))*TMath::Sqrt( rms/(x2-x1+1) - mean*mean );
  
  if (!fitSlope)
    {
      return_value.pedestal = mean;
      return_value.rms = rms;
      return_value.slope = 0;
    }
  else
    {
      double nSamples = x2-x1+1;
      double b = (nSamples*sxy-sx*sy)/(nSamples*sxx - sx*sx);
      double a = (sy - b*sx)/nSamples;

      return_value.pedestal = a;
      return_value.rms = rms;
      return_value.slope = b;

      // std::cout << "BASELINE INFOS SLOPE: " 
      // 		<< "\tPED " << return_value.pedestal 
      // 		<< "\tRMS " << return_value.rms 
      // 		<< "\tSLOPE " << return_value.slope
      // 		<< std::endl;
    }
      
  return return_value;

};
