function u = NDI_Euler(att,omega,T_d,J,n,v)

% Matrix from kinematic equation
N = [
    1 sin(att(1))*tan(att(2)) cos(att(1))*tan(att(2));
    0 cos(att(1)) -sin(att(1));
    0 sin(att(1))/cos(att(2)) cos(att(1))/cos(att(2));
    ];
f11 = N; % Kinematics equation
f12 = n/cos(att(2))*[
                    sin(att(3));
                    cos(att(2))*cos(att(3));      
                    sin(att(2))*sin(att(3));
                   ]; % Gravity gradient term
g1 = zeros(3,1);

% - Dynamics equation
omega_mat = [
            0 -omega(3) omega(2);
            omega(3) 0 -omega(1);
            -omega(2) omega(1) 0;
            ]; % rotational velocity made for the external product
f21 = -J^-1*omega_mat*J; % Dynamics equation
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

% - Computation of M = A
d_N_2 = N;
d_N_1 = [
        (cos(att(1))*omega(2) - sin(att(1))*omega(3))*tan(att(2)) (sin(att(1))*omega(2) + cos(att(1))*omega(3))/(cos(att(2)))^2 0; 
        -sin(att(1))*omega(2) - cos(att(1))*omega(3) 0 0;
        (cos(att(1))*omega(2) - sin(att(1))*omega(3))/(cos(att(2))) (sin(att(1))*omega(2) + cos(att(1)*omega(3))*tan(att(2)))/cos(att(2)) 0;
        ];
d_N = [d_N_1 d_N_2];

M = d_N*[
        zeros(3);
        J^-1;
        ];

% - Computation of l = b
l_aux = [ 
        zeros(3) f11; 
        zeros(3) f21;
        ]*[att; omega] + [g1; g2] + [f12; f22];
l = d_N*l_aux;

% - NDI controller
u = M^-1*(v - l);

end