#!/usr/bin/env python

# Vehicle Model: Bicycle
# Tire Model: Brush

__author__ = "Satya Prasad"
__email__ = "szm888@psu.edu"

import numpy as np
import globals

class Bicycle:

	def __init__(self):
		self.g = 9.81 # [m/s^2]
		self.eps = 0.00000000000001
	'''
		==========================================================================
		# Purpose: Transforms variables from wheel to body coordinates
		#
		# INPUTS:
		#	wheel_variable: A 2x2 matrix of wheel coordinates.
		#	steering_angle: A 1x2 vector of steering angles. [rad]
		#
		# OUTPUTS:
		#	body_variable: A 2x2 matrix of body coordinates.
		==========================================================================
	'''
	def fcn_VD_stWheel2BodyCoordinates(self,wheel_variable,steering_angle):
		body_variable = np.empty((2,2))
		body_variable[:] = np.NaN

		body_variable[0,:] = np.multiply(wheel_variable[0,:], np.cos(steering_angle))-\
							 np.multiply(wheel_variable[1,:], np.sin(steering_angle))
		body_variable[1,:] = np.multiply(wheel_variable[0,:], np.sin(steering_angle))+\
							 np.multiply(wheel_variable[1,:], np.cos(steering_angle))

		return body_variable
	'''
		==========================================================================
		# Purpose: Computes velocities for front and rear wheels.
		#
		# INPUTS:
		#	U: Longitudinal velocity [m/s]
		#	V: Lateral velocity [m/s]
		#	r: Yaw rate [rad/s]
		#	vehicle: Dictionary containing vehicle properties.
		#
		# OUTPUTS:
		#	wheel_velocity: A 2x2 vector of wheel velocities [m/s]
		#	First row is along vehicle (X-Direction) and second row is 
		#	perpendicular to vehicle (Y-Direction) [Front, Rear]
		==========================================================================
	'''
	def fcn_VD_stWheelVelocity(self,U,V,r,vehicle):
		# Velocity of wheel along vehicle body (X-Direction) and perpendicular to vehicle body (Y-Direction)
		wheel_velocity = np.array([[U, U], [V+vehicle['a']*r, V-vehicle['b']*r]])

		return wheel_velocity
	'''
		==========================================================================
		# Purpose: Computes slip-angles for front and rear wheels.
		#
		# INPUTS:
		#	U: Longitudinal velocity [m/s]
		#	V: Lateral velocity [m/s]
		#	r: Yaw rate [rad/s]
		#	steering_angle: A 1x2 vector of steering angles [rad] [Front, Rear]
		#	vehicle: Dictionary containing vehicle properties.
		#
		# OUTPUTS:
		#	slip_angles: A 1x2 vector containing slip-angles of front and rear wheels [rad]
		==========================================================================
	'''
	def fcn_VD_stSlipAngle(self,U,V,r,steering_angle,vehicle):
		wheel_velocity = self.fcn_VD_stWheelVelocity(U,V,r,vehicle) # calculate wheel velocities
		slip_angles = np.arctan(wheel_velocity[1,:]/wheel_velocity[0,:])-steering_angle # Slip angle

		return slip_angles
	'''
		==========================================================================
		# Purpose: Computes Wheel Slip/Longitudinal Slip for front and rear wheels.
		#
		# INPUTS:
		#	U: Longitudinal velocity [m/s]
		#	V: Lateral velocity [m/s]
		#	r: Yaw rate [rad/s]
		#	wheel_angular_velocity: A 1x2 vector of wheel angular velocities [rad/s]
		#	[Front, Rear]
		#	steering_angle: A 1x2 vector of steering angles [rad] [Front, Rear]
		#	vehicle: Dictionary containing vehicle properties.
		#
		# OUTPUTS:
		#	wheel_slip: A 1x2 vector of wheel slip [No Units]
		#	[Front, Rear]
		==========================================================================
	'''
	def fcn_VD_stWheelSlip(self,U,V,r,wheel_angular_velocity,steering_angle,vehicle):
		wheel_slip = np.empty((1,2))
		wheel_slip[:] = np.NaN

		wheel_velocity = self.fcn_VD_stWheelVelocity(U,V,r,vehicle) # calculate wheel velocities
		wheel_longitudinal_velocity = np.multiply(wheel_velocity[0,:], np.cos(steering_angle))+\
									  np.multiply(wheel_velocity[1,:], np.sin(steering_angle)) # velocity of wheel along the wheel direction
		wheel_longitudinal_velocity[0==wheel_longitudinal_velocity] = self.eps

		# wheel slip is calculated only when the tire is moving in the direction it is pointed
		temp_wheel_slip = (vehicle['Re']*wheel_angular_velocity-wheel_longitudinal_velocity)/wheel_longitudinal_velocity
		wheel_slip = np.where(0<=wheel_longitudinal_velocity, temp_wheel_slip, 0)

		# Saturation
		wheel_slip[1.0<wheel_slip]  = 1.0 # maximum wheel slip
		wheel_slip[-1.0>wheel_slip] = -1.0 # minimum wheel slip

		return wheel_slip
	'''
		==========================================================================
		#	Computes Combined Slip for front and rear wheels.
		#
		# INPUTS:
		#	slip_angle: A 1x2 vector of slip-angles. [rad]
		#	[Front, Rear]
		#	wheel_slip: A 1x2 vector of wheel slip.
		#	[Front, Rear]
		#
		# OUTPUTS:
		#	combined_slip: A 2x2 matrix of combined slip. Row-1 is X and Row-2 is Y.
		#	[Front, Rear]
		==========================================================================
	'''
	def fcn_VD_stCombinedSlip(self,slip_angle,wheel_slip):
		combined_slip = np.empty((2,2)) # Initialize a variable
		combined_slip[:] = np.NaN

		denominator = wheel_slip+1 # denominator that's used in calculating combined slip
		denominator[0==denominator] = self.eps # to avoid division by zero

		combined_slip[0,:] = wheel_slip/denominator # sigma-x
		combined_slip[1,:] = np.tan(slip_angle)/denominator # sigma-y

		return combined_slip
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
	def fcn_VD_stNormalForce(self,acceleration,vehicle,road_properties,type_of_transfer):
		Fz_front = vehicle['m']*self.g*np.cos(road_properties['bank_angle'])*np.cos(road_properties['grade'])*\
				   (vehicle['b']-vehicle['h_cg']*np.tan(road_properties['grade']))/vehicle['L']
		Fz_rear  = vehicle['m']*self.g*np.cos(road_properties['bank_angle'])*np.cos(road_properties['grade'])*\
				   (vehicle['a']+vehicle['h_cg']*np.tan(road_properties['grade']))/vehicle['L']
		normal_force = np.array([Fz_front, Fz_rear]) # load transfer due to grade and bank

		if (type_of_transfer == 'longitudinal'): # Only Longitudinal Transfer
			Fz_longitudinal = vehicle['m']*acceleration[0]*(vehicle['h_cg']/vehicle['L'])*np.array([-1, 1])
			normal_force    = normal_force + Fz_longitudinal

		return normal_force
	'''
		==========================================================================
		# Computes tire forces using Combined-Slip Brush Model.
		# Ref: Tire Modeling and Friction Estimation by Jacob Svendenius
		#
		# INPUTS:
		# 	slip_angle: A 1x2 vector of slip-angles. [rad]
		#	[Front, Rear]
		#	wheel_slip: A 1x2 vector of wheel slip.
		#	[Front, Rear]
		#	normal_force: A 1x2 vector of normal forces. [N]
		#	[Front, Rear]
		#	friction_coefficient: A 1x2 vector of friction coefficients.
		#	[Front, Rear]
		#	vehicle: Dictionary containing vehicle properties
		#
		# OUTPUTS:
		#	tire_force: A 2x2 matrix of tire forces. Row-1 is X-direction and
		#	Row-2 is Y-direction.
		#	[Front, Rear]
		==========================================================================
	'''
	def fcn_VD_stTireForceBrush(self,slip_angle,wheel_slip,normal_force,friction_coefficient,vehicle):
		tire_force = np.empty((2,2)) # Initialize a variable
		tire_force[:] = np.NaN

		combined_slip = self.fcn_VD_stCombinedSlip(slip_angle,wheel_slip) # calculate combined slip

		sigma_x = combined_slip[0,:]
		sigma_y = combined_slip[1,:]
		sigma   = np.sqrt(np.square(sigma_x) + np.square(sigma_y))
		sigma[0==sigma] = self.eps # To avoid division by zero

		# eq 4.10, 4.13 on pg 44
		Psi = np.sqrt(np.square(np.multiply(vehicle['Cx'],sigma_x)) + np.square(np.multiply(vehicle['Ca'],sigma_y)))/\
			  (3*np.multiply(friction_coefficient,normal_force)) # use peak-friction here
		Psi[1<=Psi] = 1 # account for pure slip

		# Adhesion forces:
		# eq 4.14 on pg 44
		Fax = np.multiply(np.multiply(vehicle['Cx'],sigma_x),np.square(1-Psi))
		Fay = -np.multiply(np.multiply(vehicle['Ca'],sigma_y),np.square(1-Psi))

		# Slide forces:
		# eq 4.17 on pg 46
		Fsz = np.multiply(np.multiply(normal_force,np.square(Psi)),(3-2*Psi))
		# eq 4.18 on pg 46
		Fsx = np.multiply(np.multiply(sigma_x/sigma,friction_coefficient),Fsz) # use sliding-friction here
		Fsy = -np.multiply(np.multiply(sigma_y/sigma,friction_coefficient),Fsz) # use sliding-friction here

		tire_force[0,:] = Fax+Fsx
		tire_force[1,:] = Fay+Fsy

		return tire_force
	'''
		# Calculates aligning moment on both wheels.
		# Ref: https://www.tandfonline.com/doi/full/10.1080/00423114.2019.1580377
		#
		# INPUTS:
		# 	slip_angle: A 1x2 vector of slip-angles. [rad]
		#	[Front, Rear]
		#	normal_force: A 1x2 vector of normal forces.
		#	[Front, Rear]
		#	friction_coefficient: A 1x2 vector of friction coefficients.
		#	[Front, Rear]
		#	vehicle: Dictionary containing vehicle properties
		#
		# OUTPUTS:
		#	aligning_moment: A 1x2 vector of aligning moment.
		#	[Front, Rear]
	'''
	def fcn_VD_stAligningMomentBrush(self,slip_angle,normal_force,friction_coefficient,vehicle):
		temp_variable = np.multiply(vehicle['Ca'],np.tan(slip_angle))

		Xi1 = (vehicle['contact_patch_length']/3)*temp_variable
		Xi2 = ((2-vehicle['friction_ratio'])*vehicle['contact_patch_length']/3)*(temp_variable**2)
		Xi3 = ((1-2*vehicle['friction_ratio']/3)*vehicle['contact_patch_length']/3)*(temp_variable**3)
		Xi4 = ((4/27-vehicle['friction_ratio']/9)*vehicle['contact_patch_length']/3)*(temp_variable**4)
		# equation 12 on pg 9
		temp_variable2  = np.multiply(friction_coefficient,normal_force)
		aligning_moment = Xi1-\
						  Xi2/temp_variable2+\
						  Xi3/(temp_variable2**2)-\
						  Xi4/(temp_variable2**3)
		# check for saturation
		saturation_test = (3*temp_variable2)/vehicle['Ca']-abs(np.tan(slip_angle))
		aligning_moment = np.where(0>=saturation_test, 0, aligning_moment)

		return aligning_moment
	'''
		==========================================================================
		# Purpose: Calculates accelerations.
		#
		# INPUTS:
		#	wheel_force: A 2x2 matrix of tire forces. Row-1 is X-direction and
		#	Row-2 is Y-direction.
		#	[Front, Rear]
		#	wheel_torque: A 1x2 vector of wheel torque.
		#	[Front, Rear]
		#	steering_angle: A 1x2 vector of steering angles. [rad]
		#	[Front, Rear]
		#	vehicle: Dictionary containing vehicle properties.
		#	road_properties: Dictionary containing road properties.
		#
		# OUTPUTS:
		#	acceleration: A 1X5 vector of accelerations.
		#	[longitudinal acceleration, lateral acceleration, yaw acceleration,
		#	wheel angular acceleration[Front, Rear]]
		==========================================================================
	'''
	def fcn_VD_st5dofForceEquation(self,wheel_force,wheel_torque,steering_angle,vehicle,road_properties):
		body_force = self.fcn_VD_stWheel2BodyCoordinates(wheel_force,steering_angle) # find body forces
		# Force and torque balance on the vehicle body
		ax = sum(body_force[0,:])/vehicle['m']-self.g*np.cos(road_properties['bank_angle'])*np.sin(road_properties['grade']) # longitudinal acceleration
		ay = sum(body_force[1,:])/vehicle['m']+self.g*np.sin(road_properties['bank_angle']) # lateral acceleration
		drdt = ((vehicle['a']*body_force[1,0])-(vehicle['b']*body_force[1,1]))/vehicle['Izz'] # yaw acceleration
		# Torque balance on vehicle's wheels
		domegadt = (wheel_torque-vehicle['Re']*wheel_force[0,:])/vehicle['Iw'] # angular acceleration of wheel

		acceleration = np.array([ax, ay, drdt, domegadt[0], domegadt[1]])

		return acceleration
	'''
		==========================================================================
		# Purpose: Differential equation for longitudinal velocity, lateral velocity, 
		#			yaw rate, and angular velocities of wheels.
		#
		# INPUTS:
		#	y: A 1X5 vector of velocities in body coordinates.
		#	[longitudinal velocity, lateral velocity, yaw rate, wheel angular
		#	velocity[Front, Rear]]
		#	acceleration: A 1X5 vector of accelerations.
		#	[longitudinal acceleration, lateral acceleration, yaw acceleration,
		#	wheel angular acceleration[Front, Rear]]
		#
		# OUTPUTS:
		#	DvelDt: A 1X5 vector of accelerations.
		#	[longitudinal acceleration, lateral acceleration, yaw acceleration,
		#	wheel angular acceleration[Front, Rear]]
		==========================================================================
	'''
	def fcn_VD_5dofStateEquation(self,t,y,acceleration):
		dUdt = acceleration[0]+y[2]*y[1] # Vxdot
		dVdt = acceleration[1]-y[2]*y[0] # Vydot
		drdt = acceleration[2] # yaw acceleration

		DvelDt = np.array([dUdt, dVdt, drdt, acceleration[3], acceleration[4]])

		return DvelDt
	'''
		==========================================================================
		# Purpose: Calculates velocites in global coordinates.
		#
		# INPUTS:
		#	y: A 1X3 vector of global pose [X, Y, Phi] OR [East, North, Phi]
		#	U: Longitudinal velocity [m/s]
		#	V: Lateral velocity [m/s]
		#	r: Yaw rate [rad/s]
		#
		# OUTPUTS:
		#	DposeDt: A 1X3 vector of velocities in global coordinates
		==========================================================================
	'''
	def fcn_VD_Body2GlobalCoordinates(self,t,y,U,V,r):
		psi = y[2]

		dXdt    = U*np.cos(psi)-V*np.sin(psi)
		dYdt    = U*np.sin(psi)+V*np.cos(psi)
		dPhidt  = r
		DposeDt = np.array([dXdt, dYdt, dPhidt])

		return DposeDt
	'''
		==========================================================================
		# Sets-up a 5-DoF Single-Track vehicle model using Brush tire model.
		#
		# INPUTS:
		#	t: A number indicating time corresponding to y.
		#	y: A 1x10 vector of velocities and pose. [U, V, r, omega, X, Y, Phi]
		# 	steering_angle: A 1x2 vector of steering angle. [Front, Rear]
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
		#	dydt: A 1x8 vector of accelerations and velocities.
		==========================================================================
	'''
	def fcn_VD_st5dofModel(self,t,y,steering_angle,wheel_torque,\
						   vehicle,road_properties,friction_coefficient,type_of_transfer):
		global delayed_acceleration
		U = y[0] # longitudinal velocity [m/s]
		V = y[1] # lateral velocity [m/s]
		r = y[2] # yaw rate [rad/s]
		omega = y[3:5] # angular velocity of wheels [rad/s]
		pose  = y[5:8] # position and orientation of the vehicle
		
		## Slips
		# Slip Angle/Lateral Slip
		slip_angle = self.fcn_VD_stSlipAngle(U,V,r,steering_angle,vehicle)
		# Wheel Slip/Longitudinal Slip
		wheel_slip = self.fcn_VD_stWheelSlip(U,V,r,omega,steering_angle,vehicle)

		## Normal Forces
		if globals.flag_update:
			# load transfer is calculated based on the acceleration in the previous time-step 't-deltaT'
			delayed_acceleration = globals.global_acceleration
		normal_force = self.fcn_VD_stNormalForce(delayed_acceleration[0:2],vehicle,road_properties,type_of_transfer)
		
		## Tire Forces
		tire_force = self.fcn_VD_stTireForceBrush(slip_angle,wheel_slip,normal_force,friction_coefficient,vehicle)
		
		## Newtonian Dynamics
		acceleration = self.fcn_VD_st5dofForceEquation(tire_force,wheel_torque,steering_angle,vehicle,road_properties)
		DvelDt = self.fcn_VD_5dofStateEquation(t,y[0:5],acceleration)
		if globals.flag_update:
			## Aligning Moment
			globals.aligning_moment = self.fcn_VD_stAligningMomentBrush(slip_angle,normal_force,friction_coefficient,vehicle)
			# 'global_acceleration' will be updated only in the first call to the
			# function in RK4 method i.e., while computing k1 in 'fcn_VD_RungeKutta'
			globals.global_acceleration = acceleration
			globals.flag_update = False

		## Body2Global coordinates
		DposeDt = self.fcn_VD_Body2GlobalCoordinates(t,pose,U,V,r)
		
		## Write to output
		dydt = np.concatenate((DvelDt, DposeDt))

		return dydt