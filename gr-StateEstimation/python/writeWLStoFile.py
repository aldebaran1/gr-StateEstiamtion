#!/usr/bin/env python
# 
# Copyright 2016 <+YOU OR YOUR COMPANY+>.
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

import numpy
from gnuradio import gr
from datetime import datetime
import os.path

class writeWLStoFile(gr.sync_block):
    """
    docstring for block writeWLStoFile
    """

    def check_for_empty_file(self, fname):
        file_name_tmp = fname
        counter = 0
        while (os.path.isfile(file_name_tmp)):
            counter += 1
            file_name_tmp = fname + '_' + str(counter)
        return file_name_tmp

    def __init__(self, path, prefix):
        gr.sync_block.__init__(self, name="writeWLStoFile", in_sig=[(numpy.float32, 4)], out_sig=None)
        local_time = datetime.now()
        self.file_path = path + prefix
        file_name = self.file_path + str(local_time.year)[2:4] + '%(mon)02d%(day)02d' %{"mon":local_time.month, "day":local_time.day}
        file_name_use = self.check_for_empty_file(file_name)
        self.current_day = local_time.day
        self.f = open(file_name_use, 'w')


    def work(self, input_items, output_items):
        in0 = input_items[0]
        #print in0[0]
        #print in0[0][0]
        #print in0[0][2]
        
        local_time = datetime.now()
        file_name_tmp = self.file_path + str(local_time.year)[2:4] + '%(mon)02d%(day)02d' %{"mon":local_time.month, "day":local_time.day}
        if (self.current_day != local_time.day):
            self.current_day = local_time.day
            self.f = open(file_name_tmp, 'w')
        flag = 1
        time_in_centisec = local_time.hour*60*60*1000 + local_time.minute*60*1000 + local_time.second*1000 + round(local_time.microsecond/1000)
        if (int(in0[0][3]) == 1):
            self.f.write(str(time_in_centisec)[0:10] + ',' 
                + str(in0[0][0]) + ','
                + str(in0[0][1]) + ','
                + str(in0[0][2]) + ','
                + str(int(in0[0][3])) + '\n')

        return len(input_items [0])

    def close(self):
        self.f.close()
        print('Closing file')

