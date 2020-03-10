% TuDelft - Faculty of Aerospace Engineering
% Spacecraft Attitude Dynamics & Control 
% Spacecraft Attitude Dynamics & Control: Exercise
% Rohan Camlesh Chotalal -> Student Number: 4746317
% ASSIGNMENT: Assignment 2

clear all;
close all;

% -- Figure numbers chosen to compare all methods
% -- Quaternions:
fig1 = 100;
% -- Quaternion error (first 3):
fig2 = 101;
% -- Angular velocities:
fig3 = 102;
% -- Torques:
fig4 = 103;

%% CONVERSION PARAMETERS

deg2rad = pi/180;

%% DATA

h = 700e3; % [m] - orbit height
Re = 6371e3; % [m] - Radius of the Earth
Me = 5.972e24; % [kg] - Mass of the Earth
G = 6.67408e-11; % [m^3 kg^-1 s^-2] - Gravitational Constant

mu_e = G*Me;

n = sqrt(mu_e/(Re + h)^3); % [rad/s] - orbital velocity at height h

I = [
    124.531 0 0;
    0 124.586 0;
    0 0 0.704;
    ]; % [Kg m^2] - Satellite ineria property

T_d = [0.0001 0.0001 0.0001]'; % [N m] - Disturbances torques

dt = 0.1; % [s] - sample time
T = 1500; % [s] - simulation duration

att_ini = [30 30 30]'*deg2rad; % [rad] - initial attitude (psi, theta, phi)
q_ini = att2q(att_ini); % - initial quaternions

% Reference comands:
t1 = 0:dt:99.9;
att_ref_t1 = 0*ones(3,length(t1));
t2 = 100:dt:500;
att_ref_t2 = 70*deg2rad*ones(3,length(t2));
t3 = 500.1:dt:900;
att_ref_t3 = -70*deg2rad*ones(3,length(t3));
t4 = 900.1:dt:1500;
att_ref_t4 = 0*deg2rad*ones(3,length(t4)); 

t = [t1 t2 t3 t4]; % concatenate all time intervals
att_ref_t = [att_ref_t1 att_ref_t2 att_ref_t3 att_ref_t4]; % concatenate all ref. theta intervals
q_ref_t = zeros(4,size(att_ref_t,2));
for i = 1:size(att_ref_t,2)
    q_ref_t(:,i) = att2q(att_ref_t(:,i)); % compute ref. q 
end

% Plot of the reference quaternions:
fig_num = 1;
figure(fig_num);
suptitle('Quaternions - Reference commands for quaternions control');

% - q1
subplot(3,1,1);
plot(t,q_ref_t(1,:));
xlabel('t [sec]');
ylabel('q_{ref_1} [-]');

% - q2
subplot(3,1,2);
plot(t,q_ref_t(2,:));
xlabel('t [sec]');
ylabel('q_{ref_2} [-]');

% - q3
subplot(3,1,3);
plot(t,q_ref_t(3,:));
xlabel('t [sec]');
ylabel('q_{ref_3} [-]');

% - q4
fig_num = fig_num + 1;
figure(fig_num);
plot(t,q_ref_t(4,:));
xlabel('t [sec]');
ylabel('q_{ref_4} [-]');
title('Quaternions - Reference 4th quaternion command');

%% Ex_1: Model Implementation

% --> Euler Angles: 
% Implemented in function EulerKinematics.m (or EulerKinematicsMat.m) & 
% EulerDynamics.m (or EulerDynamicsMat.m)
% --> Quaternions:
% Implemented in function QuaternionsKinematicsMat.m & 
% QuaternionsDynamicsMat.m 

%% Ex_2, Ex_3 and Ex_4: Design of linear PD controller including simulations

% Variable Initialization:
rot_vel_ini = pi/6*[1 1 1]'; % [Ask the initial angular velocity]

x_Q = zeros(7,length(t));

x_Q_dot = zeros(7,length(t));

T_c_Q = zeros(3,length(t));

e_Q = zeros(4,length(t));
de_Q = zeros(4,length(t));

x_Q(:,1) = [q_ini; rot_vel_ini];

% Proportional and Derivative gains for:
% -- Quaternion model
% Kp_Q = [I(1,1); I(2,2); I(3,3)*0.5];
% Kd_Q = [I(1,1)*6; I(2,2)*6; I(3,3)];

Kp_Q = [200; 200; 0.5];
Kd_Q = [700; 500; 1];

for i = 1:length(t)-1
    % --> Current Quaternions and Angular Velocity 
    q = x_Q(1:4,i); % current Quaternions Angles
    rot_vel = x_Q(5:7,i); % current Angular Velocity
    
    % --> Kinematics Equation:
    x_Q_dot(1:4,i) = QuaternionsKinematicsMat(q,rot_vel,n);
    
    % --> Computation of the error:
    e_Q(:,i) = q_ref_t(:,i) - q; % e = q_ref - q
    
    % --> Computation of the derivative of the error:
    d_q = x_Q_dot(1:4,i);
    de_Q(:,i) = - d_q; % d_e = - de_q
    
    % ------- 2nd way to compute the quaternion errors -------
%     Q_com = [
%         q_ref_t(4) q_ref_t(3) -q_ref_t(2) -q_ref_t(1);
%         -q_ref_t(3) q_ref_t(4) q_ref_t(1) -q_ref_t(2);
%         q_ref_t(2) -q_ref_t(1) q_ref_t(4) -q_ref_t(3);
%         q_ref_t(1) q_ref_t(2) q_ref_t(3) q_ref_t(4);
%         ];
%     e_Q(:,i) = Q_com*q;
%     de_Q(:,i) = Q_com*d_q;
    % ------- 2nd way to compute the quaternion errors -------
    T_c_Q(:,i) = Kp_Q.*e_Q(1:3,i) + Kd_Q.*de_Q(1:3,i); 
    
    % --> Dynamics Equation:
    x_Q_dot(5:7,i) = QuaternionsDynamicsMat(q,rot_vel,I,n,T_d,T_c_Q(:,i));
    
    % --> Numerical Euler Integration to obtain the next state:
    x_Q(:,i+1) = x_Q(:,i) + dt*x_Q_dot(:,i);
    
    % --> Quaternions Normalization:
    % - This does not work - yields complex numbers sometimes
    %x_Q(1,i+1) = sqrt(1 - x_Q(2,i+1)^2 - x_Q(3,i+1)^2 - x_Q(4,i+1)^2);
    %x_Q(2,i+1) = sqrt(1 - x_Q(3,i+1)^2 - x_Q(4,i+1)^2 - x_Q(1,i+1)^2);
    %x_Q(3,i+1) = sqrt(1 - x_Q(4,i+1)^2 - x_Q(1,i+1)^2 - x_Q(2,i+1)^2);
    %x_Q(4,i+1) = sqrt(1 - x_Q(1,i+1)^2 - x_Q(2,i+1)^2 - x_Q(3,i+1)^2);
    % - This is the proper way, and is used in linear algebra to normalize vectors:
    x_Q(1,i+1) = x_Q(1,i+1)/sum(x_Q(1:4,i+1).^2);
    x_Q(2,i+1) = x_Q(2,i+1)/sum(x_Q(1:4,i+1).^2);
    x_Q(3,i+1) = x_Q(3,i+1)/sum(x_Q(1:4,i+1).^2);
    x_Q(4,i+1) = x_Q(4,i+1)/sum(x_Q(1:4,i+1).^2);
end

% ---- Plots:
% 1) Quaternions:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - PD -> q_1, q_2 and q_3');
subplot(3,1,1);
plot(t,x_Q(1,:),t,q_ref_t(1,:));
xlabel('t [sec]');
ylabel('q_1 [-]');
legend('Attitude','Command');

subplot(3,1,2);
plot(t,x_Q(2,:),t,q_ref_t(2,:));
xlabel('t [sec]');
ylabel('q_2 [-]');

subplot(3,1,3);
plot(t,x_Q(3,:),t,q_ref_t(3,:));
xlabel('t [sec]');
ylabel('q_3 [-]');

fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - PD -> q_4');
plot(t,x_Q(4,:),t,q_ref_t(4,:));
xlabel('t [sec]');
ylabel('q_4 [-]');
legend('Attitude','Command');

% --> Angular Velocity:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - PD -> Angular velocities');

subplot(3,1,1);
plot(t,x_Q(5,:));
xlabel('t [sec]');
ylabel('\omega_1 [rad/s]');

subplot(3,1,2);
plot(t,x_Q(6,:));
xlabel('t [sec]');
ylabel('\omega_2 [rad/s]');

subplot(3,1,3);
plot(t,x_Q(7,:));
xlabel('t [sec]');
ylabel('\omega_3 [rad/s]');

% --> System input (Control torque):
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - PD - Control torques (inputs)');

subplot(3,1,1);
plot(t,T_c_Q(1,:));
xlabel('t [sec]');
ylabel('T_c_1 [N m]');

subplot(3,1,2);
plot(t,T_c_Q(2,:));
xlabel('t [sec]');
ylabel('T_c_2 [N m]');

subplot(3,1,3);
plot(t,T_c_Q(3,:));
xlabel('t [sec]');
ylabel('T_c_3 [N m]');
 
% --> Quaternion error:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - PD - Quaternions error');

subplot(4,1,1);
plot(t,e_Q(1,:));
xlabel('t [sec]');
ylabel('e_1 [-]');

subplot(4,1,2);
plot(t,e_Q(2,:));
xlabel('t [sec]');
ylabel('e_2 [-]');

subplot(4,1,3);
plot(t,e_Q(3,:));
xlabel('t [sec]');
ylabel('e_3 [-]');

subplot(4,1,4);
plot(t,e_Q(4,:));
xlabel('t [sec]');
ylabel('e_4 [-]');

% --> Derivative of the quaternion error:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - PD - Derivative of the quaternion error');

subplot(4,1,1);
plot(t,de_Q(1,:));
xlabel('t [sec]');
ylabel('de_1/dt [s^{-1}]');

subplot(4,1,2);
plot(t,de_Q(2,:));
xlabel('t [sec]');
ylabel('de_2/dt [s^{-1}]');

subplot(4,1,3);
plot(t,de_Q(3,:));
xlabel('t [sec]');
ylabel('de_3/dt [s^{-1}]');

subplot(4,1,4);
plot(t,de_Q(4,:));
xlabel('t [sec]');
ylabel('de_4/dt [s^{-1}]');

% % ---- For comparison:
% -- q_1 q_2 and q_3:
figure(fig1);
suptitle('Quaternions -> q_1, q_2 and q_3');
subplot(3,1,1);
plot(t,q_ref_t(1,:),t,x_Q(1,:));
xlabel('t [sec]');
ylabel('q_1 [-]');
hold on;
subplot(3,1,2);
plot(t,q_ref_t(2,:),t,x_Q(2,:));
xlabel('t [sec]');
ylabel('q_2 [-]');
hold on;
subplot(3,1,3);
plot(t,q_ref_t(3,:),t,x_Q(3,:));
xlabel('t [sec]');
ylabel('q_3 [-]');
hold on;

% -- Quaternion error (first 3):
figure(fig2);
suptitle('Quaternions - Quaternions error');
subplot(3,1,1);
plot(t,e_Q(1,:));
xlabel('t [sec]');
ylabel('e_1 [-]');
hold on;
subplot(3,1,2);
plot(t,e_Q(2,:));
xlabel('t [sec]');
ylabel('e_2 [-]');
hold on;
subplot(3,1,3);
plot(t,e_Q(3,:));
xlabel('t [sec]');
ylabel('e_3 [-]');
hold on;

% -- Angular velocities:
figure(fig3);
suptitle('Quaternions -> Angular velocities');
subplot(3,1,1);
plot(t,x_Q(5,:));
xlabel('t [sec]');
ylabel('\omega_1 [rad/s]');
hold on;
subplot(3,1,2);
plot(t,x_Q(6,:));
xlabel('t [sec]');
ylabel('\omega_2 [rad/s]');
hold on;
subplot(3,1,3);
plot(t,x_Q(7,:));
xlabel('t [sec]');
ylabel('\omega_3 [rad/s]');
hold on;

% -- Torques:
figure(fig4);
suptitle('Quaternions - Control torques (inputs)');
subplot(3,1,1);
plot(t,T_c_Q(1,:));
xlabel('t [sec]');
ylabel('T_c_1 [N m]');
hold on;
subplot(3,1,2);
plot(t,T_c_Q(2,:));
xlabel('t [sec]');
ylabel('T_c_2 [N m]');
hold on;
subplot(3,1,3);
plot(t,T_c_Q(3,:));
xlabel('t [sec]');
ylabel('T_c_3 [N m]');
hold on;

%% Ex_5 to Ex_8: Design and Simulation of attitude control system of the 3 NDI variants
%% 1) NDI Approach

% Reinitialize variables:
x_Q = zeros(7,length(t));

x_Q_dot = zeros(7,length(t));

u = zeros(3,length(t));

x_Q(:,1) = [q_ini; rot_vel_ini];

% Proportional and Derivative gains for:
% -- Quaternions model
% Kp_Q = -[1 1 1]';
% Kd_Q = -[1.5 3 4.5]';

Kp_Q = -[0.5 0.5 0.5]';
Kd_Q = -[1.8 2 1.8]';

for i = 1:length(t)-1
    % --> Current Quaternions Angles and Angular Velocity 
    q = x_Q(1:4,i); % current Quaternions Angles
    rot_vel = x_Q(5:7,i); % current Angular Velocity
    
    % --> Kinematics Equation:
    x_Q_dot(1:4,i) = QuaternionsKinematicsMat(q,rot_vel,n);
    
    % --> Computation of the error:
    Q_c_mat = [
              q_ref_t(4,i) q_ref_t(3,i) -q_ref_t(2,i) -q_ref_t(1,i);
              -q_ref_t(3,i) q_ref_t(4,i) q_ref_t(1,i) -q_ref_t(2,i);
              q_ref_t(2,i) -q_ref_t(1,i) q_ref_t(4,i) -q_ref_t(3,i);
              q_ref_t(1,i) q_ref_t(2,i) q_ref_t(3,i) q_ref_t(4,i);
              ];
    e_Q(:,i) = Q_c_mat*q; % e = Q_c_mat*q
    
    % --> Computation of the derivative of the error:
    d_q = x_Q_dot(1:4,i);
    de_Q(:,i) = Q_c_mat*d_q; % d_e = Q_c_mat*d_q
    
    % --> Virtual Control choice:
    v = Kp_Q.*e_Q(1:3,i) + Kd_Q.*de_Q(1:3,i); % NOTE: d^2(q_ref)/dt^2 = 0
    
    % --> NDI control:
    u(:,i) = NDI_Quaternions(q,rot_vel,T_d,I,n,v);
    
    % --> Dynamics Equation:
    x_Q_dot(5:7,i) = QuaternionsDynamicsMat(q,rot_vel,I,n,T_d,u(:,i));
    
    % --> Numerical Euler Integration to obtain the next state:
    x_Q(:,i+1) = x_Q(:,i) + dt*x_Q_dot(:,i);     
    
    % --> Quaternions Normalization:
%     x_Q(1,i) = sqrt(1 - x_Q(2,i)^2 - x_Q(3,i)^2 - x_Q(4,i)^2);
%     x_Q(2,i) = sqrt(1 - x_Q(3,i)^2 - x_Q(4,i)^2 - x_Q(1,i)^2);
%     x_Q(3,i) = sqrt(1 - x_Q(4,i)^2 - x_Q(1,i)^2 - x_Q(2,i)^2);
%     x_Q(4,i) = sqrt(1 - x_Q(1,i)^2 - x_Q(2,i)^2 - x_Q(3,i)^2);
    x_Q(1,i+1) = x_Q(1,i+1)/sum(x_Q(1:4,i+1).^2);
    x_Q(2,i+1) = x_Q(2,i+1)/sum(x_Q(1:4,i+1).^2);
    x_Q(3,i+1) = x_Q(3,i+1)/sum(x_Q(1:4,i+1).^2);
    x_Q(4,i+1) = x_Q(4,i+1)/sum(x_Q(1:4,i+1).^2);
end

% ---- Plots:
% 1) Quaternions:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - NDI -> q_1, q_2 and q_3');
subplot(3,1,1);
plot(t,x_Q(1,:),t,q_ref_t(1,:));
xlabel('t [sec]');
ylabel('q_1 [-]');
legend('Attitude','Command');

subplot(3,1,2);
plot(t,x_Q(2,:),t,q_ref_t(2,:));
xlabel('t [sec]');
ylabel('q_2 [-]');

subplot(3,1,3);
plot(t,x_Q(3,:),t,q_ref_t(3,:));
xlabel('t [sec]');
ylabel('q_3 [-]');

fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - NDI -> q_4');
plot(t,x_Q(4,:),t,q_ref_t(4,:));
xlabel('t [sec]');
ylabel('q_4 [-]');
legend('Attitude','Command');

% --> Angular Velocity:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - NDI -> Angular velocities');

subplot(3,1,1);
plot(t,x_Q(5,:));
xlabel('t [sec]');
ylabel('\omega_1 [rad/s]');

subplot(3,1,2);
plot(t,x_Q(6,:));
xlabel('t [sec]');
ylabel('\omega_2 [rad/s]');

subplot(3,1,3);
plot(t,x_Q(7,:));
xlabel('t [sec]');
ylabel('\omega_3 [rad/s]');
 
% --> System input (Control torque):
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - NDI - Control torques (inputs)');

subplot(3,1,1);
plot(t,u(1,:));
xlabel('t [sec]');
ylabel('T_c_1 [N m]');

subplot(3,1,2);
plot(t,u(2,:));
xlabel('t [sec]');
ylabel('T_c_2 [N m]');

subplot(3,1,3);
plot(t,u(3,:));
xlabel('t [sec]');
ylabel('T_c_3 [N m]');

% --> Quaternion error:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - NDI - Quaternions error');

subplot(4,1,1);
plot(t,e_Q(1,:));
xlabel('t [sec]');
ylabel('e_1 [-]');

subplot(4,1,2);
plot(t,e_Q(2,:));
xlabel('t [sec]');
ylabel('e_2 [-]');

subplot(4,1,3);
plot(t,e_Q(3,:));
xlabel('t [sec]');
ylabel('e_3 [-]');

subplot(4,1,4);
plot(t,e_Q(4,:));
xlabel('t [sec]');
ylabel('e_4 [-]');

% --> Derivative of the quaternion error:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - NDI - Derivative of the quaternion error');

subplot(4,1,1);
plot(t,de_Q(1,:));
xlabel('t [sec]');
ylabel('de_1/dt [s^{-1}]');

subplot(4,1,2);
plot(t,de_Q(2,:));
xlabel('t [sec]');
ylabel('de_2/dt [s^{-1}]');

subplot(4,1,3);
plot(t,de_Q(3,:));
xlabel('t [sec]');
ylabel('de_3/dt [s^{-1}]');

subplot(4,1,4);
plot(t,de_Q(4,:));
xlabel('t [sec]');
ylabel('de_4/dt [s^{-1}]');

% % ---- For comparison:
% -- q_1 q_2 and q_3:
figure(fig1);
subplot(3,1,1);
plot(t,x_Q(1,:));
xlabel('t [sec]');
ylabel('q_1 [-]');
hold on;
subplot(3,1,2);
plot(t,x_Q(2,:));
xlabel('t [sec]');
ylabel('q_2 [-]');
hold on;
subplot(3,1,3);
plot(t,x_Q(3,:));
xlabel('t [sec]');
ylabel('q_3 [-]');
hold on;

% -- Quaternion error (first 3):
figure(fig2);
subplot(3,1,1);
plot(t,e_Q(1,:));
xlabel('t [sec]');
ylabel('e_1 [-]');
hold on;
subplot(3,1,2);
plot(t,e_Q(2,:));
xlabel('t [sec]');
ylabel('e_2 [-]');
hold on;
subplot(3,1,3);
plot(t,e_Q(3,:));
xlabel('t [sec]');
ylabel('e_3 [-]');
hold on;

% -- Angular velocities:
figure(fig3);
subplot(3,1,1);
plot(t,x_Q(5,:));
xlabel('t [sec]');
ylabel('\omega_1 [rad/s]');
hold on;
subplot(3,1,2);
plot(t,x_Q(6,:));
xlabel('t [sec]');
ylabel('\omega_2 [rad/s]');
hold on;
subplot(3,1,3);
plot(t,x_Q(7,:));
xlabel('t [sec]');
ylabel('\omega_3 [rad/s]');
hold on;

% -- Torques:
figure(fig4);
subplot(3,1,1);
plot(t,u(1,:));
xlabel('t [sec]');
ylabel('T_c_1 [N m]');
hold on;
subplot(3,1,2);
plot(t,u(2,:));
xlabel('t [sec]');
ylabel('T_c_2 [N m]');
hold on;
subplot(3,1,3);
plot(t,u(3,:));
xlabel('t [sec]');
ylabel('T_c_3 [N m]');
hold on;

%% 2) NDI Approach with time-scale separation

% Reinitialize variables:
x_Q = zeros(7,length(t));

x_Q_dot = zeros(7,length(t));

u = zeros(3,length(t));

rot_vel_ref = zeros(3,length(t));

x_Q(:,1) = [q_ini; rot_vel_ini];

% Proportional and Derivative gains for:
% -- Quaternions model
Kp_Q_o = -[0.9 0.8 1.0]'; % outer loop gain

Kp_Q_i = [0.4 0.6 0.7]'; % inner loop gain

for i = 1:length(t)-1
    % --> Current Quaternions and Angular Velocity 
    q = x_Q(1:4,i); % current Quaternions
    rot_vel = x_Q(5:7,i); % current Angular Velocity
    
    % --> NDI Time-Scale separation:
    [x_Q(:,i+1),x_Q_dot(:,i),e_Q(:,i),u(:,i),rot_vel_ref(:,i+1)] = NDI_TScale_Quaternions(q,rot_vel,q_ref_t(:,i),rot_vel_ref(:,i),T_d,I,n,dt,Kp_Q_o,Kp_Q_i); 
    
    % --> Quaternions Normalization:
%     x_Q(1,i) = sqrt(1 - x_Q(2,i)^2 - x_Q(3,i)^2 - x_Q(4,i)^2);
%     x_Q(2,i) = sqrt(1 - x_Q(3,i)^2 - x_Q(4,i)^2 - x_Q(1,i)^2);
%     x_Q(3,i) = sqrt(1 - x_Q(4,i)^2 - x_Q(1,i)^2 - x_Q(2,i)^2);
%     x_Q(4,i) = sqrt(1 - x_Q(1,i)^2 - x_Q(2,i)^2 - x_Q(3,i)^2);
    x_Q(1,i+1) = x_Q(1,i+1)/sum(x_Q(1:4,i+1).^2);
    x_Q(2,i+1) = x_Q(2,i+1)/sum(x_Q(1:4,i+1).^2);
    x_Q(3,i+1) = x_Q(3,i+1)/sum(x_Q(1:4,i+1).^2);
    x_Q(4,i+1) = x_Q(4,i+1)/sum(x_Q(1:4,i+1).^2);
end

% ---- Plots:
% 1) Quaternions:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - NDI + Time-Scale -> q_1, q_2 and q_3');
subplot(3,1,1);
plot(t,x_Q(1,:),t,q_ref_t(1,:));
xlabel('t [sec]');
ylabel('q_1 [-]');
legend('Attitude','Command');

subplot(3,1,2);
plot(t,x_Q(2,:),t,q_ref_t(2,:));
xlabel('t [sec]');
ylabel('q_2 [-]');

subplot(3,1,3);
plot(t,x_Q(3,:),t,q_ref_t(3,:));
xlabel('t [sec]');
ylabel('q_3 [-]');

fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - NDI + Time-Scale -> q_4');
plot(t,x_Q(4,:),t,q_ref_t(4,:));
xlabel('t [sec]');
ylabel('q_4 [-]');
legend('Attitude','Command');

% --> Angular Velocity:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - NDI + Time-Scale -> Angular velocities');

subplot(3,1,1);
plot(t,x_Q(5,:));
xlabel('t [sec]');
ylabel('\omega_1 [rad/s]');

subplot(3,1,2);
plot(t,x_Q(6,:));
xlabel('t [sec]');
ylabel('\omega_2 [rad/s]');

subplot(3,1,3);
plot(t,x_Q(7,:));
xlabel('t [sec]');
ylabel('\omega_3 [rad/s]');

% --> System input (Control torque):
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - NDI + Time-Scale - Control torques (inputs)');

subplot(3,1,1);
plot(t,u(1,:));
xlabel('t [sec]');
ylabel('T_c_1 [N m]');

subplot(3,1,2);
plot(t,u(2,:));
xlabel('t [sec]');
ylabel('T_c_2 [N m]');

subplot(3,1,3);
plot(t,u(3,:));
xlabel('t [sec]');
ylabel('T_c_3 [N m]');

% --> Quaternion error:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - NDI + Time-Scale - Quaternions error');

subplot(4,1,1);
plot(t,e_Q(1,:));
xlabel('t [sec]');
ylabel('e_1 [-]');

subplot(4,1,2);
plot(t,e_Q(2,:));
xlabel('t [sec]');
ylabel('e_2 [-]');

subplot(4,1,3);
plot(t,e_Q(3,:));
xlabel('t [sec]');
ylabel('e_3 [-]');

subplot(4,1,4);
plot(t,e_Q(4,:));
xlabel('t [sec]');
ylabel('e_4 [-]');

% --> Derivative of the quaternion error:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - NDI + Time-Scale - Derivative of the quaternion error');

subplot(4,1,1);
plot(t,de_Q(1,:));
xlabel('t [sec]');
ylabel('de_1/dt [s^{-1}]');

subplot(4,1,2);
plot(t,de_Q(2,:));
xlabel('t [sec]');
ylabel('de_2/dt [s^{-1}]');

subplot(4,1,3);
plot(t,de_Q(3,:));
xlabel('t [sec]');
ylabel('de_3/dt [s^{-1}]');

subplot(4,1,4);
plot(t,de_Q(4,:));
xlabel('t [sec]');
ylabel('de_4/dt [s^{-1}]');

% % ---- For comparison:
% -- q_1 q_2 and q_3:
figure(fig1);
subplot(3,1,1);
plot(t,x_Q(1,:));
xlabel('t [sec]');
ylabel('q_1 [-]');
hold on;
subplot(3,1,2);
plot(t,x_Q(2,:));
xlabel('t [sec]');
ylabel('q_2 [-]');
hold on;
subplot(3,1,3);
plot(t,x_Q(3,:));
xlabel('t [sec]');
ylabel('q_3 [-]');
hold on;

% -- Quaternion error (first 3):
figure(fig2);
subplot(3,1,1);
plot(t,e_Q(1,:));
xlabel('t [sec]');
ylabel('e_1 [-]');
hold on;
subplot(3,1,2);
plot(t,e_Q(2,:));
xlabel('t [sec]');
ylabel('e_2 [-]');
hold on;
subplot(3,1,3);
plot(t,e_Q(3,:));
xlabel('t [sec]');
ylabel('e_3 [-]');
hold on;

% -- Angular velocities:
figure(fig3);
subplot(3,1,1);
plot(t,x_Q(5,:));
xlabel('t [sec]');
ylabel('\omega_1 [rad/s]');
hold on;
subplot(3,1,2);
plot(t,x_Q(6,:));
xlabel('t [sec]');
ylabel('\omega_2 [rad/s]');
hold on;
subplot(3,1,3);
plot(t,x_Q(7,:));
xlabel('t [sec]');
ylabel('\omega_3 [rad/s]');
hold on;

% -- Torques:
figure(fig4);
subplot(3,1,1);
plot(t,u(1,:));
xlabel('t [sec]');
ylabel('T_c_1 [N m]');
hold on;
subplot(3,1,2);
plot(t,u(2,:));
xlabel('t [sec]');
ylabel('T_c_2 [N m]');
hold on;
subplot(3,1,3);
plot(t,u(3,:));
xlabel('t [sec]');
ylabel('T_c_3 [N m]');
hold on;

%% 3) INDI Approach (with time-scale separation)

% Reinitialize variables:
x_Q = zeros(7,length(t));

x_Q_dot = zeros(7,length(t));

u = zeros(3,length(t));

rot_vel_ref = zeros(3,length(t));

x_Q(:,1) = [q_ini; rot_vel_ini];

% Proportional and Derivative gains for:
% -- Quaternions model
Kp_Q_o = -[1 1 1]'; % outer loop gain

Kp_Q_i = [1 1 1]'; % inner loop gain

for i = 1:length(t)-1
    % --> Current Euler Angles and Angular Velocity 
    q = x_Q(1:4,i); % current Euler Angles
    rot_vel = x_Q(5:7,i); % current Angular Velocity
    
    % --> NDI Time-Scale separation:
    [x_Q(:,i+1),x_Q_dot(:,i+1),e_Q(:,i),u(:,i+1),rot_vel_ref(:,i+1)] = INDI_TScale_Quaternions(q,rot_vel,q_ref_t(:,i),rot_vel_ref(:,i),x_Q_dot(:,i),u(:,i),T_d,I,n,dt,Kp_Q_o,Kp_Q_i);    
    
    % --> Quaternions Normalization:
%     x_Q(1,i) = sqrt(1 - x_Q(2,i)^2 - x_Q(3,i)^2 - x_Q(4,i)^2);
%     x_Q(2,i) = sqrt(1 - x_Q(3,i)^2 - x_Q(4,i)^2 - x_Q(1,i)^2);
%     x_Q(3,i) = sqrt(1 - x_Q(4,i)^2 - x_Q(1,i)^2 - x_Q(2,i)^2);
%     x_Q(4,i) = sqrt(1 - x_Q(1,i)^2 - x_Q(2,i)^2 - x_Q(3,i)^2);
    x_Q(1,i+1) = x_Q(1,i+1)/sum(x_Q(1:4,i+1).^2);
    x_Q(2,i+1) = x_Q(2,i+1)/sum(x_Q(1:4,i+1).^2);
    x_Q(3,i+1) = x_Q(3,i+1)/sum(x_Q(1:4,i+1).^2);
    x_Q(4,i+1) = x_Q(4,i+1)/sum(x_Q(1:4,i+1).^2);
end

% ---- Plots:
% 1) Quaternions:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - INDI + Time-Scale -> q_1, q_2 and q_3');
subplot(3,1,1);
plot(t,x_Q(1,:),t,q_ref_t(1,:));
xlabel('t [sec]');
ylabel('q_1 [-]');
legend('Attitude','Command');

subplot(3,1,2);
plot(t,x_Q(2,:),t,q_ref_t(2,:));
xlabel('t [sec]');
ylabel('q_2 [-]');

subplot(3,1,3);
plot(t,x_Q(3,:),t,q_ref_t(3,:));
xlabel('t [sec]');
ylabel('q_3 [-]');

fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - INDI + Time-Scale -> q_4');
plot(t,x_Q(4,:),t,q_ref_t(4,:));
xlabel('t [sec]');
ylabel('q_4 [-]');
legend('Attitude','Command');

% --> Angular Velocity:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - INDI + Time-Scale -> Angular velocities');

subplot(3,1,1);
plot(t,x_Q(5,:));
xlabel('t [sec]');
ylabel('\omega_1 [rad/s]');

subplot(3,1,2);
plot(t,x_Q(6,:));
xlabel('t [sec]');
ylabel('\omega_2 [rad/s]');

subplot(3,1,3);
plot(t,x_Q(7,:));
xlabel('t [sec]');
ylabel('\omega_3 [rad/s]');

% --> System input (Control torque):
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - INDI + Time-Scale - Control torques (inputs)');

subplot(3,1,1);
plot(t,u(1,:));
xlabel('t [sec]');
ylabel('T_c_1 [N m]');

subplot(3,1,2);
plot(t,u(2,:));
xlabel('t [sec]');
ylabel('T_c_2 [N m]');

subplot(3,1,3);
plot(t,u(3,:));
xlabel('t [sec]');
ylabel('T_c_3 [N m]');

% --> Quaternion error:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - INDI + Time-Scale - Quaternions error');

subplot(4,1,1);
plot(t,e_Q(1,:));
xlabel('t [sec]');
ylabel('e_1 [-]');

subplot(4,1,2);
plot(t,e_Q(2,:));
xlabel('t [sec]');
ylabel('e_2 [-]');

subplot(4,1,3);
plot(t,e_Q(3,:));
xlabel('t [sec]');
ylabel('e_3 [-]');

subplot(4,1,4);
plot(t,e_Q(4,:));
xlabel('t [sec]');
ylabel('e_4 [-]');

% --> Derivative of the quaternion error:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Quaternions - INDI + Time-Scale - Derivative of the quaternion error');

subplot(4,1,1);
plot(t,de_Q(1,:));
xlabel('t [sec]');
ylabel('de_1/dt [s^{-1}]');

subplot(4,1,2);
plot(t,de_Q(2,:));
xlabel('t [sec]');
ylabel('de_2/dt [s^{-1}]');

subplot(4,1,3);
plot(t,de_Q(3,:));
xlabel('t [sec]');
ylabel('de_3/dt [s^{-1}]');

subplot(4,1,4);
plot(t,de_Q(4,:));
xlabel('t [sec]');
ylabel('de_4/dt [s^{-1}]');

% ------ COMPARISON OF ALL RESPONSES:
% % ---- For comparison:
% -- q_1 q_2 and q_3:
figure(fig1);
subplot(3,1,1);
plot(t,x_Q(1,:));
axis([495 550 -1.5 1.5],'auto y');
xlabel('t [sec]');
ylabel('q_1 [-]');
hold on;
subplot(3,1,2);
plot(t,x_Q(2,:));
axis([495 550 -1.5 1.5],'auto y');
xlabel('t [sec]');
ylabel('q_2 [-]');
hold on;
subplot(3,1,3);
plot(t,x_Q(3,:));
axis([495 550 -1.5 1.5],'auto y');
xlabel('t [sec]');
ylabel('q_3 [-]');
legend('Ref','PD','NDI','NDI-Time','INDI-Time');
hold on;

% -- Quaternion error (first 3):
figure(fig2);
subplot(3,1,1);
plot(t,e_Q(1,:));
axis([495 550 -1.5 1.5],'auto y');
xlabel('t [sec]');
ylabel('e_1 [-]');
hold on;
subplot(3,1,2);
plot(t,e_Q(2,:));
axis([495 550 -1.5 1.5],'auto y');
xlabel('t [sec]');
ylabel('e_2 [-]');
hold on;
subplot(3,1,3);
plot(t,e_Q(3,:));
axis([495 550 -1.5 1.5],'auto y');
xlabel('t [sec]');
ylabel('e_3 [-]');
legend('PD','NDI','NDI-Time','INDI-Time');
hold on;

% -- Angular velocities:
figure(fig3);
subplot(3,1,1);
plot(t,x_Q(5,:));
axis([495 550 -1.5 1.5],'auto y');
xlabel('t [sec]');
ylabel('\omega_1 [rad/s]');
hold on;
subplot(3,1,2);
plot(t,x_Q(6,:));
axis([495 550 -1.5 1.5],'auto y');
xlabel('t [sec]');
ylabel('\omega_2 [rad/s]');
hold on;
subplot(3,1,3);
plot(t,x_Q(7,:));
axis([495 550 -1.5 1.5],'auto y');
xlabel('t [sec]');
ylabel('\omega_3 [rad/s]');
legend('PD','NDI','NDI-Time','INDI-Time');
hold on;

% -- Torques:
figure(fig4);
subplot(3,1,1);
plot(t,u(1,:));
axis([495 550 -1.5 1.5],'auto y');
xlabel('t [sec]');
ylabel('T_c_1 [N m]');
hold on;
subplot(3,1,2);
plot(t,u(2,:));
axis([495 550 -1.5 1.5],'auto y');
xlabel('t [sec]');
ylabel('T_c_2 [N m]');
hold on;
subplot(3,1,3);
plot(t,u(3,:));
axis([495 550 -1.5 1.5],'auto y');
xlabel('t [sec]');
ylabel('T_c_3 [N m]');
legend('PD','NDI','NDI-Time','INDI-Time');
hold on;