#ifndef histoFuncT_h
#define histoFuncT_h

#include "TH1.h"
#include "TF1.h"
#include <cstdio>

class histoFuncT
{
 public:


  //! ctor                                                                                                                                            
  histoFuncT(TH1F* histo, int& tStart, int& tStop);

  //! dtor                                                                                                                                            
  ~histoFuncT();

  //! operator()                                                                                                                                      
  double operator()(double* x, double* par);


 private:

  TH1F* histo_p;
  int tStart_p;
  int tStop_p;
};

#endif
