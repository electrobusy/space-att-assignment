function [x_k1,x_k1_dot,e_o_Q,u_in,rot_vel_ref] = INDI_TScale_Quaternions(q,rot_vel,q_ref,rot_vel_ref_0,x_k0_dot,T_c,T_d,J,n,dt,Kp_o,Kp_i)

x_k1 = zeros(7,1);
x_k1_dot = zeros(7,1);

% Define N (as a function handdler, since it will be used more the once to
% compute things:
N = @(x)[
    0 x(3) -x(2) x(1);
    -x(3) 0 x(1) x(2);
    x(2) -x(1) 0 x(3);
    -x(1) -x(2) -x(3) 0;
    ]; % rotational velocity made for quaternions

% - Kinematic equation
f11 = N(rot_vel); % Kinematics equation
f12 = (1/2)*n*[
               0 0 1 0;
               0 0 0 1;
              -1 0 0 0;
               0 -1 0 0;
               ]; % Gravity gradient term
% g1 = [zeros(4,1)];

% - Dynamics equation - Inner Loop
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
f22 = -J^-1*3*n^2*a3_mat*J*a3; % Gravity gradient term
g2 = J^-1*T_d;

% NDI Control Loops:
% 1 - Computation of the outer loop error:
Q_com = [
        q_ref(4) q_ref(3) -q_ref(2) -q_ref(1);
        -q_ref(3) q_ref(4) q_ref(1) -q_ref(2);
        q_ref(2) -q_ref(1) q_ref(4) -q_ref(3);
        q_ref(1) q_ref(2) q_ref(3) q_ref(4);
        ];
e_o_Q = Q_com*q;

% --------
% [NOT USING DERIVATIVE COMPONENT FOR CONTROL - rely only on proportional control]
% 2 - Computation of the derivative of the error:
% de_o_Q = Q_com*N(rot_vel)*q;
% --------

% 3 - Virtual Control choice for the outer loop:
v_o = Kp_o.*e_o_Q(1:3);

% 4 - Calculation of rot_vel_ref as an NDI input 
N_inv = 2*[
        q(4) q(3) -q(2) -q(1);
        -q(3) q(4) q(1) -q(2);
        q(2) -q(1) q(4) -q(3);
        q(1) q(2) q(3) q(4);
        ];
b_o = f12*q;
rot_vel_ref = N_inv(1:3,1:3)*(v_o - b_o(1:3)); % u_o -> input for the outter loop
      
% 5 - Computation of the inner loop error:
e_i_Q = rot_vel_ref - rot_vel;

% --------
% [NOT USING DERIVATIVE COMPONENT FOR CONTROL - rely only on proportional control]
% 6 - Computation of the derivative of the error:
% de_i_Q = - x_k0_dot(5:7);
% --------

% 7 - Virtual Control choice for the inner loop:
d_rot_vel_ref = (rot_vel_ref - rot_vel_ref_0)/dt; % Euler derivative 
                                                  % of reference rotational velocity
v_i = Kp_i.*e_i_Q + d_rot_vel_ref; % NOTE: d(omega_ref)/dt ~= 0, because its value
                                   % depend on the outer loop.

% 8 - NDI inner loop controller (for the dynamics)
M = J^-1;
du_in = M^-1*(v_i - x_k0_dot(5:7)); % -> incremental input for the inner loop
u_in = T_c + du_in; % INDI controller

% 9 - Apply controller in the Dynamics equation:
l = f21*rot_vel + f22 + g2;
x_k1_dot(5:7) = l + J^-1*u_in;

% 10 - Numerical Euler Integration to obtain the next state for omega: 
x_k1(5:7) = rot_vel + dt*x_k1_dot(5:7);

% 11 - Apply kinematics equation to obtain the Quaternions:
x_k1_dot(1:4) = (N(x_k1(5:7)) + f12)*q;

% 12 - Numerical Euler integration to obtain the Euler angles:
x_k1(1:4) = q + dt*x_k1_dot(1:4);

end