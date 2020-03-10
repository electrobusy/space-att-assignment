function [x_k1,y,x_dot] = EulerModel(x_k0,T_d,T_c,J,n,dt)

att = x_k0(1:3);
rot_vel = x_k0(4:6);

u = [T_d; T_c];

% -->  "State" Equation:
% - Kinematic equation
N = [
    1 sin(att(1))*tan(att(2)) cos(att(1))*tan(att(2));
    0 cos(att(1)) -sin(att(1));
    0 sin(att(1))/cos(att(2)) cos(att(1))/cos(att(2));
    ];
f11 = N; % Kinematics equation
f12 = n*[
        sin(att(3))/cos(att(2));
        cos(att(3));      
        tan(att(2))*sin(att(3));
        ]; % Gravity gradient term
g1 = [zeros(3) zeros(3)];

% - Dynamics equation
rot_vel_mat = [
            0 -rot_vel(3) rot_vel(2);
            rot_vel(3) 0 -rot_vel(1);
            -rot_vel(2) rot_vel(1) 0;
            ]; % rotational velocity made for the external product
f21 = -J^-1*rot_vel_mat*J; % Dynamics equation
a3 = [
     -sin(att(2));
     sin(att(1))*cos(att(2));
     cos(att(1))*cos(att(2));
     ];
a3_mat = [
         0 -a3(3) a3(2);
         a3(3) 0 -a3(1);
         -a3(2) a3(1) 0;
         ]; 
f22 = J^-1*3*n^2*a3_mat*J*a3; % Gravity gradient term
g2 = [J^-1 J^-1];

% - Final equation:
x_dot = [
        zeros(3) f11; 
        zeros(3) f21;
        ]*x_k0 + [g1; g2]*u + [f12; f22];

% - Numerical Euler Integration to obtain the next state:
x_k1 = x_k0 + dt*x_dot;

% --> "Output" Equation:
H = [eye(3) zeros(3)];
y = H*x_k0;

end