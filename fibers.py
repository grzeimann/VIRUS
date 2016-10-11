#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fiber Class
-----------
To be used in conjuction with IFU reduction code, Panacea


"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np

__all__ = ["Fibers"]

class Fibers:
    def __init__(self, D, order=3):
        ''' 
        Initialize class
        ----------------
        :param D:
            The number of columns in the direction of the trace.
        :param order:
            The order of the polynomial used to fit the trace.                
        :init D:
            Same as described above.
        :init order:
            Same as described above.
        :init flag:
            Flag value for undefined values.
        :init trace_x:
            Columns values for trace.
        :init trace_y:
            Row values for the trace.
        :init trace: 
            Polynomial fitted values for the trace.
        :init polyvals:
            Polynomial coefficients for 
            

        '''
        self.D = D
        self.order = order
        self.flag = -99999
        self.trace_x = self.flag * np.ones((D,),dtype = np.int)
        self.trace_y = np.zeros((D,))
        self.trace = np.zeros((D,))
        self.polyvals = np.zeros((order,))
        self.norm = np.zeros((D,))
        self.deadfiber = False
        
    def fitpoly(self):
        sel = self.x != self.flag
        self.polyvals = np.polyfit(self.x[sel] / 1032.,self.y[sel],self.order)
        
    def evalpoly(self):
        self.trace = np.polyval(self.polyvals, np.arange(len(self.x)) / 1032.)