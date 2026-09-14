% Solve for equivalent rotation and translation matrices
deltaTranslation = [2; 3];
deltaRotation = 10*pi/180;

% Fill in values for 2D XY rotation-translation
roll = 0; % in rad
pitch = 0; % in rad
yaw = deltaRotation; % in rad

% Perform trig calculations
cosR = cos(roll);
sinR = sin(roll);
cosP = cos(pitch);
sinP = sin(pitch);
cosY = cos(yaw);
sinY = sin(yaw);

Rx = [
    1 0 0;
    0 cosR -sinR;
    0 sinR  cosR;
    ];

Ry = [
    cosP 0 sinP;
    0    1    0;
    sinP 0 cosP;
    ];

Rz = [
    cosY -sinY 0;
    sinY  cosY 0;
    0     0    1;
    ];
Rotation = Rx*Ry*Rz;

tx = deltaTranslation(1);
ty = deltaTranslation(2);
tz = 0;
translation = [tx; ty; tz];

Tmatrix = [...
    Rotation translation;
    zeros(1,3) 1;
    ]

Ttemp = se3([yaw, pitch, roll],"eul",'ZYX',translation')

% Create the translation matrix
translation_matrix = makehgtform('translate',translation);
% Create the z-rotate matrix
rotation_matrix_z = makehgtform('zrotate',yaw);
% Create the y-rotate matrix
rotation_matrix_y = makehgtform('yrotate',pitch);
% Create the x-rotate matrix
rotation_matrix_x = makehgtform('xrotate',roll);
% Compute the transformation matrix
Transformation_Matrix = translation_matrix*rotation_matrix_z*...
                                rotation_matrix_y*rotation_matrix_x