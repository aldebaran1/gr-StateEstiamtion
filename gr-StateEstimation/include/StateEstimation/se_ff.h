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


#ifndef INCLUDED_STATEESTIMATION_SE_FF_H
#define INCLUDED_STATEESTIMATION_SE_FF_H

#include <StateEstimation/api.h>
#include <gr_block.h>

namespace gr {
  namespace StateEstimation {

    /*!
     * \brief <+description of block+>
     * \ingroup StateEstimation
     *
     */
    class STATEESTIMATION_API se_ff : virtual public gr_block
    {
    public:
       typedef boost::shared_ptr<se_ff> sptr;

       /*!
        * \brief Return a shared_ptr to a new instance of StateEstimation::se_ff.
        *
        * To avoid accidental use of raw pointers, StateEstimation::se_ff's
        * constructor is in a private implementation
        * class. StateEstimation::se_ff::make is the public interface for
        * creating new instances.
        */
       static sptr make(int size, int smp_rate);
    };

  } // namespace StateEstimation
} // namespace gr

#endif /* INCLUDED_STATEESTIMATION_SE_FF_H */

