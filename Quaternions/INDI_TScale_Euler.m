function [x_k1,x_k1_dot,e_o_E,u_in,rot_vel_ref] = INDI_TScale_Euler(att,rot_vel,att_ref,rot_vel_ref_0,x_k0_dot,T_c,T_d,J,n,dt,Kp_o,Kp_i)

x_k1 = zeros(6,1);
x_k1_dot = zeros(6,1);

% Define N (as a function handdler, since it will be used more the once to
% compute things:
N = @(x) ([
    1 sin(x(1))*tan(x(2)) cos(x(1))*tan(x(2));
    0 cos(x(1)) -sin(x(1));
    0 sin(x(1))/cos(x(2)) cos(x(1))/cos(x(2));
    ]);

% - Kinematic equation
f11 = N(att); % Kinematics equation
f12 = n/cos(att(2))*[
                    sin(att(3));
                    cos(att(2))*cos(att(3));      
                    sin(att(2))*sin(att(3));
                   ]; % Gravity gradient term
% g1 = [zeros(3,1)];

% - Dynamics equation - Inner Loop
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
g2 = J^-1*T_d;

% NDI Control Loops:
% 1 - Computation of the outer loop error:
e_o_E = att_ref - att;

% --------
% [NOT USING DERIVATIVE COMPONENT FOR CONTROL - rely only on proportional control]
% 2 - Computation of the derivative of the error:
% de_o_E = - x_k0_dot(1:3);
% --------

% 3 - Virtual Control choice for the outer loop:
v_o = Kp_o.*e_o_E;

% 4 - Calculation of rot_vel_ref as an NDI input 
N_inv = [
        1 0 -sin(att(2));
        0 cos(att(1)) sin(att(1))*cos(att(2));
        0 -sin(att(1)) cos(att(1))*cos(att(2));
        ];

rot_vel_ref = N_inv*(v_o - f12); % u_out -> input for the outter loop
      
% 5 - Computation of the inner loop error:
e_i_E = rot_vel_ref - rot_vel;

% --------
% [NOT USING DERIVATIVE COMPONENT FOR CONTROL - rely only on proportional control]
% 6 - Computation of the derivative of the error:
% de_i_E = - x_k0_dot(4:6);
% --------

% 7 - Virtual Control choice for the inner loop:
d_rot_vel_ref = (rot_vel_ref - rot_vel_ref_0)/dt; % Euler derivative 
                                                 % of reference rotational velocity
v_i = Kp_i.*e_i_E;% + d_rot_vel_ref; % NOTE: d(omega_ref)/dt ~= 0, because its value
                                   % depend on the outer loop.

% 8 - NDI inner loop controller (for the dynamics)
M = J^-1;
du_in = M^-1*(v_i - x_k0_dot(4:6)); % -> incremental input for the inner loop
u_in = T_c + du_in; % INDI controller

% 9 - Apply controller in the Dynamics equation:
l = f21*rot_vel + f22 + g2;
x_k1_dot(4:6) = l + J^-1*u_in;

% 10 - Numerical Euler Integration to obtain the next state for omega: 
x_k1(4:6) = rot_vel + dt*x_k1_dot(4:6);

% 11 - Apply kinematics equation to obtain the Euler rates:
x_k1_dot(1:3) = N(att)*x_k1(4:6);

% 12 - Numerical Euler integration to obtain the Euler angles:
x_k1(1:3) = att + dt*x_k1_dot(1:3);

end