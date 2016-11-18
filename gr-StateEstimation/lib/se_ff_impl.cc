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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gr_io_signature.h>
#include "se_ff_impl.h"
#include <vector>
#include <math.h>
#include <armadillo>
#include <typeinfo>
#include <iostream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <ctime>

#include <boost/date_time/posix_time/posix_time.hpp>

using namespace arma;
using namespace std;

namespace gr {
  namespace StateEstimation {

    se_ff::sptr
    se_ff::make(int size, int smp_rate)
    {
      return gnuradio::get_initial_sptr
        (new se_ff_impl(size, smp_rate));
    }

    /*
     * The private constructor
     */
    se_ff_impl::se_ff_impl(int size, int smp_rate)
      : gr_block("se_ff",
		      gr_make_io_signature(1, 1, sizeof (float)*size),
		      gr_make_io_signature(1, 1, sizeof (float)*4))
    {
        vector_size = size;
        sample_rate = smp_rate;
    }

    /*
     * Our virtual destructor.
     */
    se_ff_impl::~se_ff_impl()
    {
    }

    void
    se_ff_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      /* <+forecast+> e.g. ninput_items_required[0] = noutput_items */
      ninput_items_required[0] = noutput_items; // + history() - 1;
    }

    int
    se_ff_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        float *in0 = (float *) input_items[0];
        //const float *in0 = (const float *) &((const float *)input_items[0])[history()-1]; //begin in0 with new data

        float *out0 = (float *) output_items[0];
        //float *out1 = (float *) output_items[1]; 
        //float *out2 = (float *) output_items[2];        
        const double pi = 3.14159265;
        double delta_time = (1.0 / sample_rate);

        // Input vector!
        std::vector<float> v;
        for (int i = 0;i < vector_size;i++){
          v.push_back(in0[i]);}
        vec input_data;
        input_data = conv_to<mat>::from(v);
        //Time vector
        std::vector<double> t;
        for(int i = 0; i < input_data.size(); i++){
          t.push_back(i*delta_time);
        }
        mat time = conv_to<mat>::from(t);
        //Get F_estimate
        /*
        std::vector<double> a;
        int co = 1;
        for(int i = 0; i < 50; i++){
          if ( ((input_data(i) > 0) && (input_data(i+1) < 0)) || ((input_data(i) < 0 ) && (input_data(i+1) > 0))){
            a.push_back(i); co++;
          }
        }        
        double f_est = 1 / (delta_time * (a[1] - a[0]) * 2);
        */
        //Starting vector
        mat finals(3,1); finals.fill(0);
        mat temp(3,1); temp.fill(0);
        std::string line;
        std::ifstream infile;
        infile.open("/home/sebastijan/wls/bound_ele");
        for (int i = 0;i < 3; i++){
          getline(infile,line);
          finals(i,0) = atof(line.c_str());
          temp(i,0) = atof(line.c_str());
        }
        finals(1,0) = finals(1,0)*2*pi;

        //Inputs
        mat inp, inpt;
        int N = vector_size;

        //NWLS algoritem
        //finals = nwls(input_data, time, finals);
        int count = 0;
        float error = 0.00001;
        mat delta_x(3,1);
        delta_x(0,0) = 10.0;
        delta_x(1,0) = 10.0;
        delta_x(2,0) = 10.0;
        int A = 0;
        //while (std::abs(delta_x.max()) > error){
        while (std::abs(delta_x(2,0)) > error){
          count++;
          vec h = h_x(time, finals);
          vec hdx = h_dx(time, finals);
          vec hdw = h_dw(time, finals);
          vec hdphi = h_dphi(time, finals);
          mat H = join_rows(hdx, hdw);
          H = join_rows(H, hdphi);
          mat differ = input_data - h;
          mat Ht = transpose(H);
          mat g = multiplyMatrix(Ht, differ);
          mat G = multiplyMatrix(Ht, H);
          mat L = cholesky(G);
          mat Lt = transpose(L);
          mat u = forward_solve(L, g);          
          delta_x = backward_solve(Lt, u);
          finals = finals + delta_x;
          if (count > 6){
            if ( (std::abs(delta_x(0,0)) > error) && 
                 (std::abs(delta_x(2,0)) > error) && 
                 (std::abs(delta_x(1,0)) > 0.001) ){
              A = 1;
            }
            break;
          }
        }
        double P = 10*log10(pow((finals(0,0)*sqrt(2)),2)/0.1);

        //phase correction
        double Tau = sample_rate / (finals(1,0)/2/pi);
        double skip_mistake = fmod(vector_size, Tau) * ( (2*pi) / Tau);
        double phi_corr = finals(2,0) + skip_mistake;
        //BLEF
        double dif = std::abs(finals(0,0)) - std::abs(temp(0,0));
        //TIME STAMP
        // Get current time from the clock, using microseconds resolution
        const boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();
        // Get the time offset in current day
        const boost::posix_time::time_duration td = now.time_of_day();
        const long stamp = td.total_milliseconds();                            
        cout << "a = " << finals(0,0) << " |f = " << finals(1,0)/2/pi << " |phi = " << finals(2,0) << " |P= " << P << "           " << A << " | " << delta_x(0,0) << " | " << delta_x(1,0) << " | " << delta_x(2,0) << " | " << dif << endl;

        //out0[0] = std::abs(finals(0));
        //out1[0] = std::abs(finals(1)/2/pi);
        //out2[0] = fmod(finals(2), 2 * pi);

        double f_diff = std::abs(temp(1,0) - finals(1,0));
        
        if ( (std::abs(delta_x(0,0)) < error) && 
             (std::abs(delta_x(2,0)) < error) && 
             (std::abs(delta_x(1,0)) < 0.0001) ){
          FILE *fw;
          //cout << "BINGO!^" << endl;
          //FILE *fwp;
          //fwp = fopen("/home/satprosi-station/Dropbox/wls/Astra3B_converter_512_160304_1", "a");
          fw = fopen("/home/sebastijan/wls/bound_ele", "w");
          fprintf(fw, "%12.8f\n%12.8f\n%12.8f\n", std::abs(finals(0,0)), std::abs(finals(1,0)/2/pi), finals(2,0)); //fmod(finals(2)+skip_mistake, 2 * pi)
          //fprintf(fwp, "%d, %12.5f, %12.3f, %12.8f, %12.5f, ||, %12.8f\n", stamp, std::abs(finals(0,0)), std::abs(finals(1,0)/2/pi), finals(2,0), P, delta_x(1,0)); //fmod(finals(2)+skip_mistake, 2 * pi)
          fclose(fw);
          //fclose(fwp);
          out0[0] = std::abs(finals(0,0));
          out0[1] = std::abs(finals(1,0)/2/pi);
          out0[2] = std::abs(finals(2,0));
          out0[3] = 1;
        }
        else{
          out0[0] = std::abs(finals(0,0));
          out0[1] = std::abs(finals(1,0)/2/pi);
          out0[2] = std::abs(finals(2,0));
          out0[3] = 0; 
        }
        //cout << out0[0] << "," << out0[1] << "," << out0[2] << endl;

        consume_each (noutput_items); // Tell runtime system how many output items we produced.        
        return noutput_items;
    }

    mat se_ff_impl::transpose(mat M){
      int n = M.n_rows;
      int m = M.n_cols;
      mat Mt(m,n); Mt.fill(0);
      for (int i = 0; i < m; i++){
        for(int j = 0; j < n; j++)
          Mt(i,j) = M(j,i);
      }
      return Mt;
    }

    mat se_ff_impl::multiplyMatrix(mat A, mat B){
      int n_row = A.n_rows;
      int multiply = A.n_cols;
      int n_col = B.n_cols;
      mat C(n_row, n_col); C.fill(0);
      for (int row = 0; row < n_row; row++){
        for (int col = 0; col < n_col; col++){
          for (int inner = 0; inner < multiply; inner++){
            C(row, col) += A(row, inner) * B(inner, col);
          }
        }
      }
      return C;
    }

    mat se_ff_impl::cholesky(mat A){
      int n = A.n_rows;
      mat L(n,n); L.fill(0);
      for (int i = 0; i < n; i++){
        for (int j = 0; j < (i+1); j++){
          double s = 0;
          for (int k = 0; k < j; k++){
            s += L(i,k) * L(j,k);
          }
            L(i,j) = (i == j) ? sqrt(A(i,i) - s) : ((A(i,j) - s) / L(j,j));
        }
      }
      return L;
    }

    mat se_ff_impl::forward_solve(mat G, mat g){
      int n = g.n_rows;
      mat u(n,1); u.fill(0);
      for (int i = 0; i < n; i++){
        u(i,0) = g(i,0);
        for (int j = 0; j <= (i-1); j++){
          u(i,0) = u(i,0) - G(i,j) * u(j,0);
        }
        u(i,0) = u(i,0) / G(i,i);
      }
      return u;
    }

    mat se_ff_impl::backward_solve(mat G, mat u){
      int n = u.n_rows;
      mat x(n,1); x.fill(0);
      for (int i = n-1; i >= 0; i--){
        x(i,0) = u(i,0);
        for (int j = (i+1); j <= n-1; j++){
          x(i,0) = x(i,0) - G(i,j) * x(j,0);
        }
        x(i,0) = x(i,0) / G(i,i);
      }
      return x;
    }
    
    vec se_ff_impl::h_x(mat t, vec x){
      vec hx(t.size()); hx.fill(0);
      for (int i = 0; i < t.size(); i++){
        hx(i) = x(0) * cos(x(1) * t(i) + x(2));
      } return hx;
    }

    vec se_ff_impl::h_dx(mat t, vec x){
      vec hdx(vector_size); hdx.fill(0);
      for (int i = 0; i < vector_size; i++){
        hdx(i) = cos(x(1) * t(i) + x(2));
      } return hdx;
    }

    vec se_ff_impl::h_dw(mat t, vec x){
      vec hdw(t.size()); hdw.fill(0);
      for (int i = 0; i < t.size(); i++){
        hdw(i) = (-x(0) * t(i)) * sin(x(1) * t(i) + x(2));
      } return hdw;
    }

    vec se_ff_impl::h_dphi(mat t, vec x){
      vec hdphi(t.size()); hdphi.fill(0);
      for (int i = 0; i < t.size(); i++){
        hdphi(i) = -x(0) * sin(x(1) * t(i) + x(2));
      } return hdphi;
    }
    
  } /* namespace StateEstimation */
} /* namespace gr */

