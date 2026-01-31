function tire_force = fcn_VD_stTireForceBrush(slip_angle,wheel_slip,...
                        normal_force,friction_coefficient,vehicle)
%% fcn_VD_stTireForceBrush
%   This function computes tire forces using Combined-Slip Brush Model.
%   Ref: Tire Modeling and Friction Estimation by Jacob Svendenius
%   It uses single-track/bicycle vehicle model and brush tire model.
%
% FORMAT: 
%
%   tire_force = fcn_VD_stTireForceBrush(slip_angle,wheel_slip,...
%                   normal_force,friction_coefficient,vehicle)
%
% INPUTS:
%
%   slip_angle: A 2x1 vector of slip-angles. [rad]
%   [Front; Rear]
%   wheel_slip: A 2x1 vector of wheel slip.
%   [Front; Rear]
%   normal_force: A 2x1 vector of normal forces. [N]
%   [Front; Rear]
%   friction_coefficient: A 2x1 vector of friction coefficients.
%   [Front; Rear]
%   vehicle: MATLAB structure containing vehicle properties
%
% OUTPUTS:
%
%   tire_force: A 2x2 matrix of tire forces. Column-1 is X-direction and
%   Column-2 is Y-direction.
%   [Front; Rear]
%
% This function was written on 2021/08/13 by Satya Prasad
% Questions or comments? szm888@psu.edu

flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1, 'STARTING function: %s, in file: %s\n', st(1).name, st(1).file);
end

%% Check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _       
%  |_   _|                 | |      
%    | |  _ __  _ __  _   _| |_ ___ 
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |                  
%              |_| 
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_check_inputs
    % Are there the right number of inputs?
    if 5~=nargin
        error('Incorrect number of input arguments.')
    end
    

	% Check the slip_angle input to be sure it has 2 col, 1+ rows
	fcn_DebugTools_checkInputsToFunctions(slip_angle, '2column_of_numbers',[1 2]);

	% Check the normal_force input to be sure it has 2 col, 1+ rows
	fcn_DebugTools_checkInputsToFunctions(normal_force, '2column_of_numbers',[1 2]);

	% Check the friction_coefficient input to be sure it has 2 col, 1+ rows
	fcn_DebugTools_checkInputsToFunctions(friction_coefficient, '2column_of_numbers',[1 2]);

	% Check the wheel_slip input to be sure it has 2 col, 1+ rows
	fcn_DebugTools_checkInputsToFunctions(wheel_slip, '2column_of_numbers',[1 2]);

end

%% Calculate Tire Forces using Brush Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
combined_slip = fcn_VD_stCombinedSlip(slip_angle,wheel_slip); % calculate combined slip

sigma_x = combined_slip(:,1);
sigma_y = combined_slip(:,2);
sigma   = sqrt(sigma_x.^2 + sigma_y.^2);
sigma(0==sigma) = eps; % To avoid division by zero

% eq 4.10, 4.13 on pg 44
Psi = sqrt((vehicle.Cx.*sigma_x).^2 + (vehicle.Ca.*sigma_y).^2)./...
        abs(3*friction_coefficient.*normal_force); % use peak-friction here
Psi(1<=Psi) = 1; % account for pure slip

% Adhesion forces:
% eq 4.14 on pg 44
Fax = vehicle.Cx.*sigma_x.*((1-Psi).^2);
Fay = -vehicle.Ca.*sigma_y.*((1-Psi).^2);

% Slide forces:
% eq 4.17 on pg 46
Fsz = normal_force.*(Psi.^2).*(3-2*Psi);
% eq 4.18 on pg 46
Fsx = (sigma_x./sigma).*friction_coefficient.*Fsz; % use sliding-friction here
Fsy = -(sigma_y./sigma).*friction_coefficient.*Fsz; % use sliding-friction here

tire_force = [Fax+Fsx, Fay+Fsy];

%% Any debugging?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _                 
%  |  __ \     | |                
%  | |  | | ___| |__  _   _  __ _ 
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_do_debug
    fprintf(1, 'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end

end