function u = NDI_Quaternions(q,rot_vel,T_d,J,n,v)

% - Kinematic equation
N = (1/2)*[
          0 rot_vel(3) -rot_vel(2) rot_vel(1);
          -rot_vel(3) 0 rot_vel(1) rot_vel(2);
          rot_vel(2) -rot_vel(1) 0 rot_vel(3);
          -rot_vel(1) -rot_vel(2) -rot_vel(3) 0;
          ]; % rotational velocity made for quaternions
f11 = N; % Kinematics equation
f12 = (1/2)*n*[0 0 1 0;
               0 0 0 1;
              -1 0 0 0;
               0 -1 0 0;
              ]; % Gravity gradient term
g1 = zeros(4,1);

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
f22 = J^-1*3*n^2*a3_mat*J*a3; % Gravity gradient term
g2 = J^-1*T_d;

% - Computation of M = A
d_N_1 = f11 + f12;
d_N_2 = (1/2)*[
               q(4) -q(3) q(2);
               q(3) q(4) -q(1);
               -q(2) q(1) q(4);
               -q(1) -q(2) -q(3);
              ];
d_N = [d_N_1 d_N_2]; % last equation is ommited due to NDI constraint

M = d_N*[
         zeros(4,3);
         J^-1;
        ]; 

% - Computation of l = b
l_aux = [ 
        f11+f12 zeros(4,3); 
        zeros(3,4) f21;
        ]*[q; rot_vel] + [g1; g2] + [zeros(4,1); f22];
l = d_N*l_aux; % last equation is ommited due to NDI constraint

% - NDI controller
u = M(1:3,:)^-1*(v - l(1:3));

% NOTE: We must use m outputs to control the m inputs. The last quaternion
% is a combination of the other 3.

end