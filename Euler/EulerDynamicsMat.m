% Complete Dynamics (matrix notation)
function omega_dot = EulerDynamicsMat(theta,omega,I,n,T_d,T_c)

% - Angular velocity term of the spacecraft:
omega_mat = [
    0 -omega(3) omega(2);
    omega(3) 0 -omega(1);
    -omega(2) omega(1) 0;
    ];
f21 = -I^-1*omega_mat*I*omega;

% - Gravity gradient term:
a3 = [
    -sin(theta(2));
    sin(theta(1))*cos(theta(2));
    cos(theta(1))*cos(theta(2));
    ];
a3_mat = [
    0 -a3(3) a3(2);
    a3(3) 0 -a3(1);
    -a3(2) a3(1) 0;
    ];
f22 = I^-1*3*n^2*a3_mat*I*a3;

% - Sum of previous terms gives total angular velocity of the Spacecraft:
omega_dot = f21 + f22 + I^-1*T_d + I^-1*T_c;

end