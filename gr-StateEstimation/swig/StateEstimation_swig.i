/* -*- c++ -*- */

#define STATEESTIMATION_API

%include "gnuradio.i"			// the common stuff

//load generated python docstrings
%include "StateEstimation_swig_doc.i"

%{
#include "StateEstimation/se_ff.h"
%}


%include "StateEstimation/se_ff.h"
GR_SWIG_BLOCK_MAGIC2(StateEstimation, se_ff);
