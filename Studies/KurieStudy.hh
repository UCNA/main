#ifndef KURIE_STUDY_HH
#define KURIE_STUDY_HH

#include "FloatErr.hh"
#include <TH2F.h>

float_err testKurie(); 
float_err makeKurieFit(TH1F* hIn);
float_err makeKurieFitsforRun(std::string runnum);

#endif
