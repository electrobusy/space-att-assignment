% TuDelft - Faculty of Aerospace Engineering
% Spacecraft Attitude Dynamics & Control 
% Spacecraft Attitude Dynamics & Control: Exercise
% Rohan Camlesh Chotalal -> Student Number: 4746317
% ASSIGNMENT: Assignment 2

clear all;
close all;

% -- Figure numbers chosen to compare all methods
% -- Euler angles:
fig1 = 100;
% -- Euler angle error (first 3):
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
% q_ini = att2q(att_ini); % - initial quaternions

% Reference commands:
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
% q_ref_t = att2q(att_ref_t); % compute ref. q --> IN OTHER SCRIPT

% % Plot of the reference angles:
fig_num = 1;
% figure(fig_num);
% suptitle('Euler - Reference commands for attitude control');
% subplot(3,1,1);
% plot(t,att_ref_t(1,:));
% xlabel('t [sec]');
% ylabel('\theta_{ref_1} [rad]');
% 
% subplot(3,1,2);
% plot(t,att_ref_t(2,:));
% xlabel('t [sec]');
% ylabel('\theta_{ref_2} [rad]');
% 
% subplot(3,1,3);
% plot(t,att_ref_t(3,:));
% xlabel('t [sec]');
% ylabel('\theta_{ref_3} [rad]');


%% Ex_1: Model Implementation

% --> Euler Angles: 
% Implemented in function EulerKinematics.m (or EulerKinematicsMat.m) & 
% EulerDynamics.m (or EulerDynamicsMat.m)

%% Ex_2, Ex_3 and Ex_4: Design of linear PD controller including simulations

% Variable Initialization:
rot_vel_ini = pi/6*[1 1 1]'; % [Ask the initial angular velocity]

x_E = zeros(6,length(t));

x_E_dot = zeros(6,length(t));

T_c_E = zeros(3,length(t));

e_E = zeros(3,length(t));
de_E = zeros(3,length(t));

x_E(:,1) = [att_ini; rot_vel_ini];

% Proportional and Derivative gains for:
% -- Euler model
% Kp_E = [I(1,1); I(2,2); I(3,3)*0.5];
% Kd_E = [I(1,1)*6; I(2,2)*6; I(3,3)];

Kp_E = [200; 200; 0.5];
Kd_E = [700; 700; 1];

for i = 1:length(t)-1
    % --> Current Euler Angles and Angular Velocity 
    att = x_E(1:3,i); % current Euler Angles
    rot_vel = x_E(4:6,i); % current Angular Velocity
    
    % --> Kinematics Equation:
    x_E_dot(1:3,i) = EulerKinematicsMat(att,rot_vel,n);
    
    % --> Computation of the error:
    e_E(:,i) = att_ref_t(:,i) - att; % e = theta_ref - theta
    
    % --> Computation of the derivative of the error:
    d_att = x_E_dot(1:3,i);
    de_E(:,i) = - d_att; % d_e = - de_theta
    T_c_E(:,i) = Kp_E.*e_E(:,i) + Kd_E.*de_E(:,i); 
    
    % --> Dynamics Equation:
    x_E_dot(4:6,i) = EulerDynamicsMat(att,rot_vel,I,n,T_d,T_c_E(:,i));
    
    % - Numerical Euler Integration to obtain the next state:
    x_E(:,i+1) = x_E(:,i) + dt*x_E_dot(:,i);
end

% % ---- Plots:
% 1) Euler:
% --> Angles:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Euler - PD -> Attitude');
subplot(3,1,1);
plot(t,x_E(1,:),t,att_ref_t(1,:));
xlabel('t [sec]');
ylabel('\theta_1 [rad]');
legend('Attitude','Command');

subplot(3,1,2);
plot(t,x_E(2,:),t,att_ref_t(2,:));
xlabel('t [sec]');
ylabel('\theta_2 [rad]');

subplot(3,1,3);
plot(t,x_E(3,:),t,att_ref_t(3,:));
xlabel('t [sec]');
ylabel('\theta_3 [rad]');


% --> Angular Velocity:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Euler - PD -> Angular velocities');

subplot(3,1,1);
plot(t,x_E(4,:));
xlabel('t [sec]');
ylabel('\omega_1 [rad/s]');

subplot(3,1,2);
plot(t,x_E(5,:));
xlabel('t [sec]');
ylabel('\omega_2 [rad/s]');

subplot(3,1,3);
plot(t,x_E(6,:));
xlabel('t [sec]');
ylabel('\omega_3 [rad/s]');

% --> System input (Control torque):
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Euler - PD - Control torques (inputs)');

subplot(3,1,1);
plot(t,T_c_E(1,:));
xlabel('t [sec]');
ylabel('T_c_1 [N m]');

subplot(3,1,2);
plot(t,T_c_E(2,:));
xlabel('t [sec]');
ylabel('T_c_2 [N m]');

subplot(3,1,3);
plot(t,T_c_E(3,:));
xlabel('t [sec]');
ylabel('T_c_3 [N m]');

% --> Attitude error:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Euler - PD - Attitude error');

subplot(3,1,1);
plot(t,e_E(1,:));
xlabel('t [sec]');
ylabel('e_1 [rad]');

subplot(3,1,2);
plot(t,e_E(2,:));
xlabel('t [sec]');
ylabel('e_2 [rad]');

subplot(3,1,3);
plot(t,e_E(3,:));
xlabel('t [sec]');
ylabel('e_3 [rad]');

% --> Derivative of the attitude error:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Euler - PD - Derivative of the attitude error');

subplot(3,1,1);
plot(t,de_E(1,:));
xlabel('t [sec]');
ylabel('de_1/dt [rad/s]');

subplot(3,1,2);
plot(t,de_E(2,:));
xlabel('t [sec]');
ylabel('de_2/dt [rad/s]');

subplot(3,1,3);
plot(t,de_E(3,:));
xlabel('t [sec]');
ylabel('de_3/dt [rad/s]');

% % ---- For comparison:
% -- Euler angles:
figure(fig1);
suptitle('Euler -> Attitude');
subplot(3,1,1);
plot(t,att_ref_t(1,:),t,x_E(1,:));
xlabel('t [sec]');
ylabel('\theta_1 [-]');
hold on;
subplot(3,1,2);
plot(t,att_ref_t(2,:),t,x_E(2,:));
xlabel('t [sec]');
ylabel('\theta_2 [-]');
hold on;
subplot(3,1,3);
plot(t,att_ref_t(3,:),t,x_E(3,:));
xlabel('t [sec]');
ylabel('\theta_3 [-]');
hold on;

% -- Euler angles error:
figure(fig2);
suptitle('Euler - Attitude error');
subplot(3,1,1);
plot(t,e_E(1,:));
xlabel('t [sec]');
ylabel('e_1 [-]');
hold on;
subplot(3,1,2);
plot(t,e_E(2,:));
xlabel('t [sec]');
ylabel('e_2 [-]');
hold on;
subplot(3,1,3);
plot(t,e_E(3,:));
xlabel('t [sec]');
ylabel('e_3 [-]');
hold on;

% -- Angular velocities:
figure(fig3);
suptitle('Euler -> Angular velocities');
subplot(3,1,1);
plot(t,x_E(4,:));
xlabel('t [sec]');
ylabel('\omega_1 [rad/s]');
hold on;
subplot(3,1,2);
plot(t,x_E(5,:));
xlabel('t [sec]');
ylabel('\omega_2 [rad/s]');
hold on;
subplot(3,1,3);
plot(t,x_E(6,:));
xlabel('t [sec]');
ylabel('\omega_3 [rad/s]');
hold on;

% -- Torques:
figure(fig4);
suptitle('Euler - Control torques (inputs)');
subplot(3,1,1);
plot(t,T_c_E(1,:));
xlabel('t [sec]');
ylabel('T_c_1 [N m]');
hold on;
subplot(3,1,2);
plot(t,T_c_E(2,:));
xlabel('t [sec]');
ylabel('T_c_2 [N m]');
hold on;
subplot(3,1,3);
plot(t,T_c_E(3,:));
xlabel('t [sec]');
ylabel('T_c_3 [N m]');
hold on;


%% Ex_5 to Ex_8: Design and Simulation of attitude control system of the 3 NDI variants
%% 1) NDI Approach

% Reinitialize variables:
x_E = zeros(6,length(t));

x_E_dot = zeros(6,length(t));

u = zeros(3,length(t));

x_E(:,1) = [att_ini; rot_vel_ini];

% Proportional and derivative gains for:
% -- Euler model
Kp_E = [1 1 1]';
Kd_E = [2 2 3]';

for i = 1:length(t)-1
    % --> Current Euler Angles and Angular Velocity 
    att = x_E(1:3,i); % current Euler Angles
    rot_vel = x_E(4:6,i); % current Angular Velocity
    
    % --> Kinematics Equation:
    x_E_dot(1:3,i) = EulerKinematicsMat(att,rot_vel,n);
    
    % --> Computation of the error:
    e_E(:,i) = att_ref_t(:,i) - att; % e = theta_ref - theta
    
    % --> Computation of the derivative of the error:
    d_att = x_E_dot(1:3,i);
    de_E(:,i) = - d_att; % d_e = - de_theta
    
    % --> Virtual Control choice:
    v = Kp_E.*e_E(:,i) + Kd_E.*de_E(:,i); % NOTE: d^2(theta_ref)/dt^2 = 0
    
    % --> NDI control:
    u(:,i) = NDI_Euler(att,rot_vel,T_d,I,n,v);
    
    % --> Dynamics Equation:
    x_E_dot(4:6,i) = EulerDynamicsMat(att,rot_vel,I,n,T_d,u(:,i));
    
    % - Numerical Euler Integration to obtain the next state:
    x_E(:,i+1) = x_E(:,i) + dt*x_E_dot(:,i);     
end

% % ---- Plots:
% 1) Euler:
% --> Angles:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Euler - NDI -> Attitude');
subplot(3,1,1);
plot(t,x_E(1,:),t,att_ref_t(1,:));
xlabel('t [sec]');
ylabel('\theta_1 [rad]');
legend('Attitude','Command')

subplot(3,1,2);
plot(t,x_E(2,:),t,att_ref_t(2,:));
xlabel('t [sec]');
ylabel('\theta_2 [rad]');

subplot(3,1,3);
plot(t,x_E(3,:),t,att_ref_t(3,:));
xlabel('t [sec]');
ylabel('\theta_3 [rad]');

% --> Angular Velocity:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Euler - NDI -> Angular velocities');

subplot(3,1,1);
plot(t,x_E(4,:));
xlabel('t [sec]');
ylabel('\omega_1 [rad/s]');

subplot(3,1,2);
plot(t,x_E(5,:));
xlabel('t [sec]');
ylabel('\omega_2 [rad/s]');

subplot(3,1,3);
plot(t,x_E(6,:));
xlabel('t [sec]');
ylabel('\omega_3 [rad/s]');

% --> System input (Control torque):
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Euler - NDI - Control torques (inputs)');

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

% --> Attitude error:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Euler - NDI - Attitude error');

subplot(3,1,1);
plot(t,e_E(1,:));
xlabel('t [sec]');
ylabel('e_1 [rad]');

subplot(3,1,2);
plot(t,e_E(2,:));
xlabel('t [sec]');
ylabel('e_2 [rad]');

subplot(3,1,3);
plot(t,e_E(3,:));
xlabel('t [sec]');
ylabel('e_3 [rad]');

% --> Derivative of the attitude error:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Euler - NDI - Derivative of the attitude error');

subplot(3,1,1);
plot(t,de_E(1,:));
xlabel('t [sec]');
ylabel('de_1/dt [rad/s]');

subplot(3,1,2);
plot(t,de_E(2,:));
xlabel('t [sec]');
ylabel('de_2/dt [rad/s]');

subplot(3,1,3);
plot(t,de_E(3,:));
xlabel('t [sec]');
ylabel('de_3/dt [rad/s]');

% % ---- For comparison:
% -- Euler angles:
figure(fig1);
subplot(3,1,1);
plot(t,x_E(1,:));
xlabel('t [sec]');
ylabel('\theta_1 [-]');
hold on;
subplot(3,1,2);
plot(t,x_E(2,:));
xlabel('t [sec]');
ylabel('\theta_2 [-]');
hold on;
subplot(3,1,3);
plot(t,x_E(3,:));
xlabel('t [sec]');
ylabel('\theta_3 [-]');
hold on;

% -- Euler angle Error:
figure(fig2);
subplot(3,1,1);
plot(t,e_E(1,:));
xlabel('t [sec]');
ylabel('e_1 [-]');
hold on;
subplot(3,1,2);
plot(t,e_E(2,:));
xlabel('t [sec]');
ylabel('e_2 [-]');
hold on;
subplot(3,1,3);
plot(t,e_E(3,:));
xlabel('t [sec]');
ylabel('e_3 [-]');
hold on;

% -- Angular velocities:
figure(fig3);
subplot(3,1,1);
plot(t,x_E(4,:));
xlabel('t [sec]');
ylabel('\omega_1 [rad/s]');
hold on;
subplot(3,1,2);
plot(t,x_E(5,:));
xlabel('t [sec]');
ylabel('\omega_2 [rad/s]');
hold on;
subplot(3,1,3);
plot(t,x_E(6,:));
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
x_E = zeros(6,length(t));

x_E_dot = zeros(6,length(t));

u = zeros(3,length(t));

rot_vel_ref = zeros(3,length(t));

x_E(:,1) = [att_ini; rot_vel_ini];

% Proportional and Derivative gains for:
% -- Euler model
Kp_E_o = [1 1 0.5]'; % outer loop gain

Kp_E_i = [10 10 5]'; % inner loop gain

for i = 1:length(t)-1
    % --> Current Euler Angles and Angular Velocity 
    att = x_E(1:3,i); % current Euler Angles
    rot_vel = x_E(4:6,i); % current Angular Velocity
    
    % --> NDI Time-Scale separation:
    [x_E(:,i+1),x_E_dot(:,i),e_E(:,i),u(:,i),rot_vel_ref(:,i+1)] = NDI_TScale_Euler(att,rot_vel,att_ref_t(:,i),rot_vel_ref(:,i),T_d,I,n,dt,Kp_E_o,Kp_E_i);    
end

% % ---- Plots:
% 1) Euler:
% --> Angles:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Euler - NDI + Time-Scale -> Attitude');
subplot(3,1,1);
plot(t,x_E(1,:),t,att_ref_t(1,:));
xlabel('t [sec]');
ylabel('\theta_1 [rad]');
legend('Attitude','Command')

subplot(3,1,2);
plot(t,x_E(2,:),t,att_ref_t(2,:));
xlabel('t [sec]');
ylabel('\theta_2 [rad]');

subplot(3,1,3);
plot(t,x_E(3,:),t,att_ref_t(3,:));
xlabel('t [sec]');
ylabel('\theta_3 [rad]');

% --> Angular Velocity:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Euler - NDI + Time-Scale -> Angular velocities');

subplot(3,1,1);
plot(t,x_E(4,:));
xlabel('t [sec]');
ylabel('\omega_1 [rad/s]');

subplot(3,1,2);
plot(t,x_E(5,:));
xlabel('t [sec]');
ylabel('\omega_2 [rad/s]');

subplot(3,1,3);
plot(t,x_E(6,:));
xlabel('t [sec]');
ylabel('\omega_3 [rad/s]');

% --> System input (Control torque):
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Euler - NDI + Time-Scale - Control torques (inputs)');

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

% --> Attitude error:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Euler - NDI + Time-Scale - Attitude error');

subplot(3,1,1);
plot(t,e_E(1,:));
xlabel('t [sec]');
ylabel('e_1 [rad]');

subplot(3,1,2);
plot(t,e_E(2,:));
xlabel('t [sec]');
ylabel('e_2 [rad]');

subplot(3,1,3);
plot(t,e_E(3,:));
xlabel('t [sec]');
ylabel('e_3 [rad]');

% --> Derivative of the attitude error:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Euler - NDI + Time-Scale - Derivative of the attitude error');

subplot(3,1,1);
plot(t,de_E(1,:));
xlabel('t [sec]');
ylabel('de_1/dt [rad/s]');

subplot(3,1,2);
plot(t,de_E(2,:));
xlabel('t [sec]');
ylabel('de_2/dt [rad/s]');

subplot(3,1,3);
plot(t,de_E(3,:));
xlabel('t [sec]');
ylabel('de_3/dt [rad/s]');

% % ---- For comparison:
% -- Euler angles:
figure(fig1);
subplot(3,1,1);
plot(t,x_E(1,:));
xlabel('t [sec]');
ylabel('\theta_1 [-]');
hold on;
subplot(3,1,2);
plot(t,x_E(2,:));
xlabel('t [sec]');
ylabel('\theta_2 [-]');
hold on;
subplot(3,1,3);
plot(t,x_E(3,:));
xlabel('t [sec]');
ylabel('\theta_3 [-]');
hold on;

% -- Euler angle error:
figure(fig2);
subplot(3,1,1);
plot(t,e_E(1,:));
xlabel('t [sec]');
ylabel('e_1 [-]');
hold on;
subplot(3,1,2);
plot(t,e_E(2,:));
xlabel('t [sec]');
ylabel('e_2 [-]');
hold on;
subplot(3,1,3);
plot(t,e_E(3,:));
xlabel('t [sec]');
ylabel('e_3 [-]');
hold on;

% -- Angular velocities:
figure(fig3);
subplot(3,1,1);
plot(t,x_E(4,:));
xlabel('t [sec]');
ylabel('\omega_1 [rad/s]');
hold on;
subplot(3,1,2);
plot(t,x_E(5,:));
xlabel('t [sec]');
ylabel('\omega_2 [rad/s]');
hold on;
subplot(3,1,3);
plot(t,x_E(6,:));
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
x_E = zeros(6,length(t));

x_E_dot = zeros(6,length(t));

u = zeros(3,length(t));

rot_vel_ref = zeros(3,length(t));

x_E(:,1) = [att_ini; rot_vel_ini];

% Proportional and Derivative gains for:
% -- Euler model
Kp_E_o = [0.5 0.5 1]'; % outer loop gain

Kp_E_i = [1 1 2]'; % inner loop gain

for i = 1:length(t)-1
    % --> Current Euler Angles and Angular Velocity 
    att = x_E(1:3,i); % current Euler Angles
    rot_vel = x_E(4:6,i); % current Angular Velocity
    
    % --> INDI Time-Scale separation:
    [x_E(:,i+1),x_E_dot(:,i+1),e_E(:,i),u(:,i+1),rot_vel_ref(:,i+1)] = INDI_TScale_Euler(att,rot_vel,att_ref_t(:,i),rot_vel_ref(:,i),x_E_dot(:,i),u(:,i),T_d,I,n,dt,Kp_E_o,Kp_E_i);    
end

% % ---- Plots:
% 1) Euler:
% --> Angles:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Euler - INDI + Time-Scale -> Attitude');
subplot(3,1,1);
plot(t,x_E(1,:),t,att_ref_t(1,:));
xlabel('t [sec]');
ylabel('\theta_1 [rad]');
legend('Attitude','Command')

subplot(3,1,2);
plot(t,x_E(2,:),t,att_ref_t(2,:));
xlabel('t [sec]');
ylabel('\theta_2 [rad]');

subplot(3,1,3);
plot(t,x_E(3,:),t,att_ref_t(3,:));
xlabel('t [sec]');
ylabel('\theta_3 [rad]');

% --> Angular Velocity:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Euler - INDI + Time-Scale -> Angular velocities');

subplot(3,1,1);
plot(t,x_E(4,:));
xlabel('t [sec]');
ylabel('\omega_1 [rad/s]');

subplot(3,1,2);
plot(t,x_E(5,:));
xlabel('t [sec]');
ylabel('\omega_2 [rad/s]');

subplot(3,1,3);
plot(t,x_E(6,:));
xlabel('t [sec]');
ylabel('\omega_3 [rad/s]');

% --> System input (Control torque):
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Euler - INDI + Time-Scale - Control torques (inputs)');

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

% --> Attitude error:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Euler - INDI + Time-Scale - Attitude error');

subplot(3,1,1);
plot(t,e_E(1,:));
xlabel('t [sec]');
ylabel('e_1 [rad]');

subplot(3,1,2);
plot(t,e_E(2,:));
xlabel('t [sec]');
ylabel('e_2 [rad]');

subplot(3,1,3);
plot(t,e_E(3,:));
xlabel('t [sec]');
ylabel('e_3 [rad]');

% --> Derivative of the attitude error:
fig_num = fig_num + 1;
figure(fig_num);
suptitle('Euler - INDI + Time-Scale - Derivative of the attitude error');

subplot(3,1,1);
plot(t,de_E(1,:));
xlabel('t [sec]');
ylabel('de_1/dt [rad/s]');

subplot(3,1,2);
plot(t,de_E(2,:));
xlabel('t [sec]');
ylabel('de_2/dt [rad/s]');

subplot(3,1,3);
plot(t,de_E(3,:));
xlabel('t [sec]');
ylabel('de_3/dt [rad/s]');

% ------ COMPARISON OF ALL RESPONSES:
% % ---- For comparison:
% -- theta_1 theta_2 and theta_3:
figure(fig1);
subplot(3,1,1);
plot(t,x_E(1,:));
axis([495 550 -1.5 1.5],'auto y');
xlabel('t [sec]');
ylabel('\theta_1 [-]');
hold on;
subplot(3,1,2);
plot(t,x_E(2,:));
axis([495 550 -1.5 1.5],'auto y');
xlabel('t [sec]');
ylabel('\theta_2 [-]');
hold on;
subplot(3,1,3);
plot(t,x_E(3,:));
axis([495 550 -1.5 1.5],'auto y');
xlabel('t [sec]');
ylabel('\theta_3 [-]');
legend('Ref','PD','NDI','NDI-Time','INDI-Time');
hold on;

% -- Euler angle error:
figure(fig2);
subplot(3,1,1);
plot(t,e_E(1,:));
axis([495 550 -1.5 1.5],'auto y');
xlabel('t [sec]');
ylabel('e_1 [-]');
hold on;
subplot(3,1,2);
plot(t,e_E(2,:));
axis([495 550 -1.5 1.5],'auto y');
xlabel('t [sec]');
ylabel('e_2 [-]');
hold on;
subplot(3,1,3);
plot(t,e_E(3,:));
axis([495 550 -1.5 1.5],'auto y');
xlabel('t [sec]');
ylabel('e_3 [-]');
legend('PD','NDI','NDI-Time','INDI-Time');
hold on;

% -- Angular velocities:
figure(fig3);
subplot(3,1,1);
plot(t,x_E(4,:));
axis([495 550 -1.5 1.5],'auto y');
xlabel('t [sec]');
ylabel('\omega_1 [rad/s]');
hold on;
subplot(3,1,2);
plot(t,x_E(5,:));
axis([495 550 -1.5 1.5],'auto y');
xlabel('t [sec]');
ylabel('\omega_2 [rad/s]');
hold on;
subplot(3,1,3);
plot(t,x_E(6,:));
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