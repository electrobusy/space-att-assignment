function [x_k1,y,x_dot] = QuaternionsModel(x_k0,T_d,T_c,J,n,dt)

q = x_k0(1:4);
rot_vel = x_k0(5:7);

u = [T_d; T_c];

% -->  "State" Equation:
% - Kinematic equation
N = (1/2)*[
    0 rot_vel(3) -rot_vel(2) rot_vel(1);
    -rot_vel(3) 0 rot_vel(1) rot_vel(2);
    rot_vel(2) -rot_vel(1) 0 rot_vel(3);
    -rot_vel(1) -rot_vel(2) -rot_vel(3) 0;
    ]; % rotational velocity made for quaternions
f11 = N; % Kinematics equation
f12 = 2*n*[0 0 1 0;
           0 0 0 1;
           -1 0 0 0;
           0 -1 0 0;
           ]*q; % Gravity gradient term
g1 = [zeros(4,3) zeros(4,3)];

% - Dynamics equation
rot_vel_mat = [
              0 -rot_vel(3) rot_vel(2);
              rot_vel(3) 0 -rot_vel(1);
              -rot_vel(2) rot_vel(1) 0;
              ]; % rotational velocity made for the external product
f21 = -J^-1*rot_vel_mat*J; % Dynamics equation
a3 = [
     2*(q(1)*q(3) - q(2)*q(4));
     2*(q(2)*q(3) + q(1)*q(4));
     (1 - 2*(q(1)^2 + q(2)^2));
     ];
a3_mat = [
         0 -a3(3) a3(2);
         a3(3) 0 -a3(1);
         -a3(2) a3(1) 0;
         ];
f22 = -J^-1*3*n^3*a3_mat*J*a3; % Gravity gradient term
g2 = [J^-1 J^-1];

% - Final equation:
x_dot = [
        f11 zeros(4,3); 
        zeros(3,4) f21
        ]*x_k0 + [g1; g2]*u + [f12; f22];

% - Numerical Euler Integration to obtain the next state:
x_k1 = x_k0 + dt*x_dot;

% --> "Output" Equation:
H = [eye(4) zeros(4,3)];
y = H*x_k0;

end