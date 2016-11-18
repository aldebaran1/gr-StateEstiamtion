/* -*- c++ -*- */
/* 
 * Copyright 2015 <+YOU OR YOUR COMPANY+>.
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifndef INCLUDED_STATEESTIMATION_SE_FF_IMPL_H
#define INCLUDED_STATEESTIMATION_SE_FF_IMPL_H

#include <StateEstimation/se_ff.h>
#include <armadillo>

//using namespace arma;

namespace gr {
  namespace StateEstimation {

    class se_ff_impl : public se_ff
    {
    private:
      int vector_size;
      //double phi;
      double sample_rate;
      //arma::vec nwls(arma::mat &z, arma::mat &t, arma::vec &X);
      arma::mat multiplyMatrix(arma::mat A, arma::mat B);
      arma::mat transpose(arma::mat M);
      arma::mat cholesky(arma::mat A);
      arma::mat forward_solve(arma::mat G, arma::mat g);
      arma::mat backward_solve(arma::mat G, arma::mat u);
      arma::vec h_x(arma::mat t, arma::vec x);
      arma::vec h_dx(arma::mat t, arma::vec x);
      arma::vec h_dw(arma::mat t, arma::vec x);
      arma::vec h_dphi(arma::mat t, arma::vec x);
    public:
      se_ff_impl(int size, int smp_rate);
      ~se_ff_impl();

      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      // Where all the action really happens
      int general_work(int noutput_items,
		       gr_vector_int &ninput_items,
		       gr_vector_const_void_star &input_items,
		       gr_vector_void_star &output_items);
    };

  } // namespace StateEstimation
} // namespace gr

#endif /* INCLUDED_STATEESTIMATION_SE_FF_IMPL_H */

