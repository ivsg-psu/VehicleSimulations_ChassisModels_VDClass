#!/usr/bin/env python
# Vehicle Model: Kinematic

__author__ = "Satya Prasad"
__email__ = "szm888@psu.edu"

import numpy as np
import globals

class Kinematic:

	def __init__(self):
		self.g = 9.81 # [m/s^2]
	'''
		==========================================================================
		# Purpose: Calculates longitudinal acceleration.
		#
		# INPUTS:
		#	U: Longitudinal velocity [m/s]
		#	normal_force: A 1x2 vector of normal forces.
		#	[Front, Rear]
		#	steering_angle: Front steering angle. [rad]
		#	wheel_torque: A 1x2 vector of wheel torque.
		#	[Front, Rear]
		#	vehicle: Dictionary containing vehicle properties.
		#	road_properties: Dictionary containing road properties.
		#	friction_coefficient: A 1x2 vector of friction coefficients.
		#	[Front, Rear]
		#
		# OUTPUTS:
		#	dUdt: Longitudinal acceleration [m/s^2]
		==========================================================================
	'''
	def fcn_VD_LongitudinalModel(self,U,normal_force,steering_angle,wheel_torque,\
								 vehicle,road_properties,friction_coefficient):
		tractive_force = wheel_torque/vehicle['Re'] # both driving and brake force
		rolling_resistance_force = 0.015*normal_force # rolling resistance
		gravity_force = vehicle['m']*self.g*np.cos(road_properties['bank_angle'])*np.sin(road_properties['grade'])
		drag_force = 0.0 # drag force

		dUdt = (sum(tractive_force-rolling_resistance_force)-gravity_force-drag_force)/vehicle['m']

		if 0.0 < U and 0.0 < dUdt:
			dUdt = 0.0

		return dUdt
	'''
		==========================================================================
		# Purpose: Calculates normal force on both wheels.
		#
		# INPUTS:
		#	acceleration: A 1x2 vector of accelerations. [m/s^2]
		#	[Longitudinal; Lateral] Note: Lateral acceleration is not used.
		#	vehicle: Dictionary containing vehicle properties.
		#	road_properties: Dictionary containing road properties.
		#	type_of_transfer: To decide type of load transfer.
		#		'longitudinal': Only longitudinal weight transfer
		#		'default': No weight transfer. Any string argument will give this result.
		#
		# OUTPUTS:
		#	normal_force: A 1x2 vector of normal forces.
		#	[Front, Rear]
		==========================================================================
	'''
	def fcn_VD_kNormalForce(self,acceleration,vehicle,road_properties,type_of_transfer):
		g = 9.81 # [m/s^2]
		Fz_front = vehicle['m']*g*np.cos(road_properties['bank_angle'])*np.cos(road_properties['grade'])*\
				   (vehicle['b']-vehicle['h_cg']*np.tan(road_properties['grade']))/vehicle['L']
		Fz_rear  = vehicle['m']*g*np.cos(road_properties['bank_angle'])*np.cos(road_properties['grade'])*\
				   (vehicle['a']+vehicle['h_cg']*np.tan(road_properties['grade']))/vehicle['L']
		normal_force = np.array([Fz_front, Fz_rear]) # load transfer due to grade and bank

		if (type_of_transfer == 'longitudinal'): # Only Longitudinal Transfer
			Fz_longitudinal = vehicle['m']*acceleration[0]*(vehicle['h_cg']/vehicle['L'])*np.array([-1, 1])
			normal_force    = normal_force + Fz_longitudinal

		return normal_force
	'''
		==========================================================================
		# Purpose: Calculates yaw rate.
		#
		# INPUTS:
		#	U: Longitudinal velocity [m/s]
		#	steering_angle: Front steering angle. [rad]
		#	vehicle: Dictionary containing vehicle properties.
		#
		# OUTPUTS:
		#	r: Yaw rate [rad/s]
		==========================================================================
	'''
	def fcn_VD_kYawRate(self,U,steering_angle,vehicle):
		if 0.0 != steering_angle:
			turn_radius = vehicle['d']/steering_angle # turn radius
			r = U/turn_radius # yaw rate
		else:
			r = 0.0

		return r
	'''
		==========================================================================
		# Purpose: Calculates velocites in global coordinates.
		#
		# INPUTS:
		#	y: A 1X3 vector of global pose [X, Y, Phi] OR [East, North, Phi]
		#	U: Longitudinal velocity [m/s]
		#	r: Yaw rate [rad/s]
		#
		# OUTPUTS:
		#	DposeDt: A 1X3 vector of velocities in global coordinates
		==========================================================================
	'''
	def fcn_VD_Body2GlobalCoordinates(self,t,y,U,r):
		psi = y[2]

		dXdt    = U*np.cos(psi)
		dYdt    = U*np.sin(psi)
		dPhidt  = r
		DposeDt = np.array([dXdt, dYdt, dPhidt])

		return DposeDt
	'''
		==========================================================================
		# Sets-up a 5-DoF Single-Track vehicle model using Brush tire model.
		#
		# INPUTS:
		#	t: A number indicating time corresponding to y.
		#	y: A 1x4 vector of velocities and pose. [U, X, Y, Phi]
		# 	steering_angle: Front steering angle.
		#	wheel_torque: A 1x2 vector of wheel torque. [Front, Rear]
		#	vehicle: Dictionary containing vehicle properties.
		#	road_properties: Dictionary containing road properties.
		#	friction_coefficient: A 1x2 vector of friction coefficients.
		#	[Front, Rear]
		#	type_of_transfer: To decide type of load transfer.
		#		'longitudinal': Only longitudinal weight transfer
		#		'default': No weight transfer. Any string argument will give this result.
		#
		# OUTPUTS:
		#	dydt: A 1x4 vector of accelerations and velocities.
		==========================================================================
	'''
	def fcn_VD_KinematicModel(self,t,y,steering_angle,wheel_torque,\
						      vehicle,road_properties,friction_coefficient,type_of_transfer):
		global delayed_acceleration
		U = y[0] # longitudinal velocity [m/s]
		pose  = y[1:4] # position and orientation of the vehicle
		
		## Normal Forces
		if globals.flag_update:
			# load transfer is calculated based on the acceleration in the previous time-step 't-deltaT'
			delayed_acceleration = globals.global_acceleration
		normal_force = self.fcn_VD_kNormalForce(delayed_acceleration[0:2],vehicle,road_properties,type_of_transfer)
		
		## Longitudinal Model
		dUdt = self.fcn_VD_LongitudinalModel(U,normal_force,steering_angle,wheel_torque,\
											 vehicle,road_properties,friction_coefficient)
		if globals.flag_update:
			## Aligning Moment
			globals.aligning_moment = np.array([2, 2])
			# 'global_acceleration' will be updated only in the first call to the
			# function in RK4 method i.e., while computing k1 in 'fcn_VD_RungeKutta'
			globals.global_acceleration = np.array([dUdt, 0, 0, 0, 0])
			globals.flag_update = False

		## Body2Global coordinates
		r = self.fcn_VD_kYawRate(U,steering_angle,vehicle) # yaw rate
		DposeDt = self.fcn_VD_Body2GlobalCoordinates(t,pose,U,r)
		
		## Write to output
		dydt = np.concatenate(([dUdt], DposeDt))

		return dydt