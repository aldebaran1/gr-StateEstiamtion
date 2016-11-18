#!/usr/bin/env python
# 
# Copyright 2015 <+YOU OR YOUR COMPANY+>.
# 
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
# 
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this software; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
# 

from gnuradio import gr, gr_unittest
import StateEstimation_swig as StateEstimation

class qa_se_ff (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def test_001_se_ff (self):
    	#input data
    	size = 512
    	signal = []
    	f = open("/home/satprosi-station/gr-sewls/python/input_test_data.txt", 'r')
    	for i in range(1, size+1):
		line = f.readline()
		signal.append(float(line))
	f.close()
	input_signal = tuple(signal)
	#src = gr.vector_source_f(input_signal, False, size)

	pwr = -6.8811
	sample_rate = 400000
	#print signal
	#print len(signal)
	#print size(input_signal)
	src = gr.vector_source_f(input_signal, False, size)	
	#src1 = gr.vector_source_f((0.1,))	
	bea = StateEstimation.se_ff(size, sample_rate)
	dst = gr.vector_sink_f()
	self.tb.connect(src, (bea,0))
	#self.tb.connect(src1, (bea,1))
	self.tb.connect((bea,0), dst)
	self.tb.run()
	result_data = dst.data()
	expected_result = (pwr,)
	self.assertFloatTuplesAlmostEqual(expected_result, result_data, 4)


if __name__ == '__main__':
    gr_unittest.run(qa_se_ff, "qa_se_ff.xml")
