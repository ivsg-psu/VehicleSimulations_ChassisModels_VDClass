%%%%%%%%%%%%%% S-Function sfcn_estimateFrictionForce_Bicycle %%%%%%%%%%%%%%
% Purpose:
%   Estimates frictional force at tire locations based on the paper
%   Rapid Road Friction Estimation using Independent Left/Right Steering Torque Measurements
%   Link: https://www.tandfonline.com/doi/full/10.1080/00423114.2019.1580377
% 
%   Note: Straight line condition gives friction force equals to NaN.
%   Note: Change the sample time depending on the simulation.
%
% INPUTS:
%   INPUT 1: Aligning moment, eq. 12 on page 9
%   INPUT 2: Xi1 in eq. 13a on page 9
%   INPUT 3: Xi2 in eq. 13a on page 9
%   INPUT 4: Xi3 in eq. 13b on page 9
%   INPUT 5: Xi3 in eq. 13b on page 9
% 
% OUTPUTS:
%   OUTPUT 1: estimate of frictional force
% 
% NOTE: coefficient of friction in the exception is same as true value as
% the computations involved are quivalent to mu*N/N
% 
% Author: Satya Prasad
% Create Date: 2020-07-26
% ======== to do list ============
% 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sfcn_estimateFrictionForce_Bicycle(block)
% The setup method is used to set up the basic attributes of the
% S-function such as ports, parameters, etc. Do not add any other
% calls to the main body of the function.
setup(block);

%% Function: setup ===================================================
%%   C MEX counterpart: mdlInitializeSizes
function setup(block)
% Register number of input and output ports
block.NumInputPorts  = 5;
block.NumOutputPorts = 1;

% Setup functional port properties to dynamically inherited
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
for i = 1:5
    block.InputPort(i).Dimensions = 2;
end

% Override output port properties
block.OutputPort(1).Dimensions = 2;

% Set block sample time to be inherited
block.SampleTimes = [-1, 0];

% Set the block simStateCompliance to default (i.e., same as a built-in block)
block.SimStateCompliance = 'DefaultSimState';

% Register methods
block.RegBlockMethod('Outputs', @Output);   % Required

function Output(block)
%% Output function with vectorization
% i = 1,2,3,4 represent different tires of a vehicle
% i = 1 -> Front left tire
% i = 2 -> Front right tire
% i = 3 -> Rear left tire
% i = 4 -> Rear right tire

Mz  = block.InputPort(1).Data;   % Aligning moment on all tires
Xi1 = block.InputPort(2).Data;   % Xi1 corresponding to all tires
Xi2 = block.InputPort(3).Data;   % Xi2 corresponding to all tires
Xi3 = block.InputPort(4).Data;   % Xi3 corresponding to all tires
Xi4 = block.InputPort(5).Data;   % Xi4 corresponding to all tires
diff_Mz_Xil = Mz - Xi1;          % Difference between Mz and Xi1
diff_Mz_Xil( diff_Mz_Xil == 0 ) = eps;  % To avoid divide by zero

% equation 20 on pg 10
p = (-3*(Mz-Xi1).*Xi3 - Xi2.^2)./(3*diff_Mz_Xil.^2);
% equation 21 on pg 10
q = (2*Xi2.^3 + 9*(Mz-Xi1).*Xi2.*Xi3 + 27*((Mz-Xi1).^2).*Xi4)./(27*diff_Mz_Xil.^3);

% indices/tires with negative 'p'
indicesNegativeP = p < 0;
indices_positive_p = 0<p; % indices/tires with positive 'p'

% p, q, Xi2, difference between Mz & Xi1 values corresponding to indices/tires with -ve p
pNegative   = p(indicesNegativeP);
qNegative   = q(indicesNegativeP);
Xi2Negative = Xi2(indicesNegativeP);
diff_Mz_XilNegative = diff_Mz_Xil(indicesNegativeP);

% p, q, Xi2, difference between Mz & Xi1 values corresponding to indices/tires with +ve p
pPositive   = p(indices_positive_p);
qPositive   = q(indices_positive_p);
Xi2Positive = Xi2(indices_positive_p);
diff_Mz_XilPositive = diff_Mz_Xil(indices_positive_p);

% initialize output variable
output_temp = NaN(2,1);

if ~isempty(pNegative)
    % equation 19 on pg 10
    frictionForceNegative = -2*sign(qNegative).*sqrt(-pNegative/3)...
        .*cosh(acosh(-1.5*(abs(qNegative)./pNegative).*sqrt(-3./pNegative))/3) ...
        - Xi2Negative./(3*diff_Mz_XilNegative);
    % write frictional force to the output temporary variable
    output_temp(indicesNegativeP) = frictionForceNegative;
    
end

if ~isempty(pPositive)
    % equation 19 on pg 10
    frictionForcePositive = -2*sqrt(pPositive/3)...
        .*sinh(asinh(1.5*(qPositive./pPositive).*sqrt(3./pPositive))/3) ...
        - Xi2Positive./(3*diff_Mz_XilPositive);
    % write frictional force to the output temporary variable
    output_temp(indices_positive_p) = frictionForcePositive;
    
end

block.OutputPort(1).Data = output_temp;
