#ifndef WAVEFORM_FIT_H
#define WAVEFORM_FIT_H

#include "TProfile.h"
#include "interface/Waveform.hpp"
#include "Math/Minimizer.h"

#include "TH1F.h"
#include "TF1.h"
#include "TVirtualFitter.h"

#include "interface/histoFuncT.h"
#include "TError.h"
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>


namespace WaveformFit
{
  void alignWaveform(TProfile* ref_profile, TProfile* fit_profile, ROOT::Math::Minimizer* &minimizer);
  void fitWaveform(Waveform* wave, TProfile* amplitudeProfile,  int x1, int x2, const Waveform::max_amplitude_informations& max, const Waveform::baseline_informations& rms, ROOT::Math::Minimizer* &minimizer);
  void fitTemplate(TH1F* waveToFit, TH1F* waveTemplate, TF1** f_template, float binWidth, float offset, int tStart, int tStop, int x1, int x2, const Waveform::max_amplitude_informations& max, int& status);

}
#endif
