#!/usr/bin/env python

__author__ = "Satya Prasad"
__email__ = "szm888@psu.edu"

import numpy as np

def initialize():
	global flag_update, global_acceleration, aligning_moment
	flag_update = True
	global_acceleration = np.zeros(5)
	aligning_moment = np.zeros(2)