#!/usr/bin/env python

__author__ = "Satya Prasad"
__email__ = "szm888@psu.edu"

import numpy as np
from kinematic import Kinematic
import globals
import csv

kinematic_model = Kinematic()
globals.initialize()

# Define a python dictionary that specifies the physical values for a vehicle.
# For convenience, we ask that you call this dictionary 'vehicle'.
vehicle = {'m': 1600.0, # mass (kg)
		   'Izz': 2500.0, # mass moment of inertia (kg m^2)
		   'Iw': 1.2, # mass moment of inertia of a wheel (kg m^2)
		   'Re': 0.32, # effective radius of a wheel (m)
		   'a': 1.3, # length from front axle to CG (m)
		   'L': 2.6, # wheelbase (m)
		   'b': 1.3, # length from rear axle to CG (m)
		   'd': 1.5, # track width (m)
		   'h_cg': 0.42, # height of the cg (m)
		   'Ca': np.array([95000.0, 110000.0]), # wheel cornering stiffnesses
		   'Cx': np.array([65000.0, 65000.0]), # longitudinal stiffnesses
		   'contact_patch_length': 0.15, # [m]
		   'friction_ratio': 1} # [No units]
		   
# road properties
road_properties = {'grade': 0.0, 'bank_angle': 0.0}
friction_coefficient = np.array([0.9, 0.9])

# load transfer conditions
type_of_transfer = 'longitudinal';

# RK4
def fcn_VD_RungeKutta(input_function, initial_states, initial_time, time_interval):
	final_time = initial_time+time_interval # find the final time

	k1 = input_function(initial_time, initial_states)
	# print(k1)
	k2 = input_function(initial_time+time_interval/2, initial_states+(time_interval/2)*k1)
	# print(k2)
	k3 = input_function(initial_time+time_interval/2, initial_states+(time_interval/2)*k2)
	k4 = input_function(final_time, initial_states+time_interval*k3)
	final_states = initial_states + (time_interval/6)*(k1+2*k2+2*k3+k4)

	return final_states

# Items used to define steering input
steering_amplitude = 2*np.pi/180 # 2 degrees of steering amplitude for input sinewave
Period = 3.0 # Units are seconds. A typical lane change is about 3 to 4 seconds based on experimental highway measurements
wheel_torque = np.array([0, 5]) # wheel torque [Nm]

# Define items used to determine how long to run sim
TotalTime = 1.5*Period # This is how long the simulation will run.
deltaT = 0.01
N_timeSteps = int(np.floor(TotalTime/deltaT)+1) # This is the number of time steps we should have

# variables to store outputs of simulation
Sim_time = np.empty((N_timeSteps,1))
Sim_time[:] = np.NaN
States = np.empty((N_timeSteps,7))
States[:] = np.NaN
Pose = np.empty((N_timeSteps,3))
Pose[:] = np.NaN

input_states = np.array([25.0, 0.0, 0.0, 25.0/vehicle['Re'], 25.0/vehicle['Re'], 0.0, 0.0, 0.0])
in_states = np.array([25.0, 0.0, 0.0, 0.0])
input_time = 0.0
for counter in range(N_timeSteps):
	Sim_time[counter]   = input_time
	States[counter,0:5] = input_states[0:5]
	Pose[counter,:]     = input_states[5:8]

	## Steering Input
	steering_angle = 2*np.pi/180 # front steering angle [rad]

	globals.flag_update = True # set it to to true before every call to RK4 method
	y = fcn_VD_RungeKutta(lambda t_start,y_in: kinematic_model.fcn_VD_KinematicModel(t_start,y_in,steering_angle,wheel_torque,\
																	                 vehicle,road_properties,friction_coefficient,type_of_transfer),\
						  					   in_states,input_time,deltaT)
	U = y[0]
	V = 0
	r = kinematic_model.fcn_VD_kYawRate(U,steering_angle,vehicle)
	omega = np.array([U, U])/vehicle['Re']
	input_states = np.concatenate(([U, V, r], omega, y[1:4]))
	in_states = y
	States[counter,5:7] = globals.global_acceleration[0:2]

	input_time = input_time+deltaT

	with open('States.csv','a') as fd1:
		writer = csv.writer(fd1)
		writer.writerow(States[counter,:])

	with open('Pose.csv','a') as fd2:
		writer = csv.writer(fd2)
		writer.writerow(Pose[counter,:])

	with open('Sim_time.csv','a') as fd3:
		writer = csv.writer(fd3)
		writer.writerow(Sim_time[counter,:])

	with open('moment.csv','a') as fd4:
		writer = csv.writer(fd4)
		writer.writerow(globals.aligning_moment)