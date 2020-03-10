% Complete Kinematics (matrix notation)
function q_dot = QuaternionsKinematicsMat(q,rot_vel,n)

% - Kinematic equation
N = (1/2)*[
    0 rot_vel(3) -rot_vel(2) rot_vel(1);
    -rot_vel(3) 0 rot_vel(1) rot_vel(2);
    rot_vel(2) -rot_vel(1) 0 rot_vel(3);
    -rot_vel(1) -rot_vel(2) -rot_vel(3) 0;
    ]; % rotational velocity made for quaternions
f11 = N; 

% - Gravity gradient term:
f12 = 2*n*[0 0 1 0;
           0 0 0 1;
           -1 0 0 0;
           0 -1 0 0;
           ]*q; % Gravity gradient term

% - Sum of previous terms gives total angular velocity of the Spacecraft:
q_dot = f11*q + f12;

end