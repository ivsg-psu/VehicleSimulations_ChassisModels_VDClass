function friction_coefficient = fcn_VD_estimateFrictionCoefficient(...
    slip_angle, normal_force, aligning_moment, ...
    vehicle)
%% fcn_VD_estimateFrictionCoefficient
%   This function calculates friction coefficient at all four wheels.
%   Reference: https://www.tandfonline.com/doi/full/10.1080/00423114.2019.1580377
%
%   Note: Straight line condition gives friction coefficient equals to NaN.
%
% FORMAT:
%
%   friction_coefficient = fcn_VD_estimateFrictionCoefficient(...
%   slip_angle, normal_force, aligning_moment, vehicle)
%
% INPUTS:
%
%   slip_angle: A 4x1 vector of slip-angles.
%   normal_force: A 4x1 vector of normal forces.
%   aligning_moment: A 4x1 vector of aligning moment.
%   vehicle: MATLAB structure containing vehicle properties.
%
% OUTPUTS:
%
%   friction_coefficient: A 4x1 vector of friction coefficient.
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
    fcn_VD_checkInputsToFunctions(aligning_moment,'vector4');
end

%% Estimate Friction Coefficient
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

diff_Mz_Xil = aligning_moment-Xi1; % Difference between aligning moment and Xi1
diff_Mz_Xil(0==diff_Mz_Xil) = eps; % To avoid division by zero

% equation 20 on pg 10
p = (-3*(aligning_moment-Xi1).*Xi3 - Xi2.^2)./(3*diff_Mz_Xil.^2);
% equation 21 on pg 10
q = (2*Xi2.^3 + 9*(aligning_moment-Xi1).*Xi2.*Xi3 + 27*((aligning_moment-Xi1).^2).*Xi4)./(27*diff_Mz_Xil.^3);

indices_negative_p = 0>p; % indices/tires with negative 'p'
indices_positive_p = 0<p; % indices/tires with positive 'p'

% p, q, Xi2, difference between Mz & Xi1 values corresponding to indices/tires with -ve p
pNegative   = p(indices_negative_p);
qNegative   = q(indices_negative_p);
Xi2Negative = Xi2(indices_negative_p);
diff_Mz_XilNegative = diff_Mz_Xil(indices_negative_p);

% p, q, Xi2, difference between Mz & Xi1 values corresponding to indices/tires with non -ve p
pPositive   = p(indices_positive_p);
pPositive(0==pPositive) = eps; % To avoid division by zero
qPositive   = q(indices_positive_p);
Xi2Positive = Xi2(indices_positive_p);
diff_Mz_XilPositive = diff_Mz_Xil(indices_positive_p);

friction_force = nan(4,1); % initialize output variable
% write frictional force to the output
if ~isempty(pNegative)
    % equation 19 on pg 10
    friction_force(indices_negative_p) = -2*sign(qNegative).*sqrt(-pNegative/3).*...
        cosh(acosh(-1.5*(abs(qNegative)./pNegative).*sqrt(-3./pNegative))/3) - ...
        Xi2Negative./(3*diff_Mz_XilNegative);
end
if ~isempty(pPositive)
    % equation 19 on pg 10
    friction_force(indices_positive_p) = -2*sqrt(pPositive/3).*...
        sinh(asinh(1.5*(qPositive./pPositive).*sqrt(3./pPositive))/3) - ...
        Xi2Positive./(3*diff_Mz_XilPositive);
end
friction_coefficient = friction_force./normal_force; % calculate friction coefficient

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