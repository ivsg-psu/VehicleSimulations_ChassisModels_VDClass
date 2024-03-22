function aligning_moment = fcn_VD_dtAligningMomentBrush(slip_angle,...
                            normal_force,friction_coefficient,vehicle)
%% fcn_VD_dtAligningMomentBrush
%   This function calculates aligning moment on all the four wheels.
%   Reference: https://www.tandfonline.com/doi/full/10.1080/00423114.2019.1580377
%   It uses double-track vehicle model and brush tire model.
%
% FORMAT:
%
%   aligning_moment = fcn_VD_dtAligningMomentBrush(slip_angle,normal_force,...
%   friction_coefficient,vehicle)
%
% INPUTS:
%
%   slip_angle: A 4x1 vector of slip-angles.
%   normal_force: A 4x1 vector of normal forces.
%   friction_coefficient: A 4x1 vector of friction coefficients.
%   vehicle: MATLAB structure containing vehicle properties.
%
% OUTPUTS:
%
%   aligning_moment: A 4x1 vector of aligning moment.
%   [Front Left; Front Right; Rear Left; Rear Right]
%
% This function was written on 2021/05/18 by Satya Prasad
% Questions or comments? szm888@psu.edu
%

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
    if 4~=nargin
        error('Incorrect number of input arguments.')
    end
    
    % Check the inputs
    fcn_VD_checkInputsToFunctions(slip_angle,'vector4');
    fcn_VD_checkInputsToFunctions(normal_force,'vector4');
    fcn_VD_checkInputsToFunctions(friction_coefficient,'vector4');
end

%% Calculate Aligning Moment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% equation 13a, 13b on pg 9
Xi1 = (vehicle.contact_patch_length*vehicle.Ca.*tan(slip_angle))/3;
Xi2 = ((2-vehicle.friction_ratio)*vehicle.contact_patch_length*(vehicle.Ca.^2).*(tan(slip_angle).^2))/3;
Xi3 = ((1-2*vehicle.friction_ratio/3)*vehicle.contact_patch_length*(vehicle.Ca.^3).*(tan(slip_angle).^3))/3;
Xi4 = ((4/27-vehicle.friction_ratio/9)*vehicle.contact_patch_length*(vehicle.Ca.^4).*(tan(slip_angle).^4))/3;
% equation 12 on pg 9
aligning_moment = Xi1-...
                  Xi2./(friction_coefficient.*normal_force)+...
                  Xi3./((friction_coefficient.*normal_force).^2)-...
                  Xi4./((friction_coefficient.*normal_force).^3);
% check for saturation
saturation_test = (3*friction_coefficient.*normal_force)./vehicle.Ca-abs(tan(slip_angle));
aligning_moment(0>=saturation_test) = 0;

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