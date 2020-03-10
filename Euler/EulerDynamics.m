% Complete Dynamics (- no matrix notation)
function omega_dot = EulerDynamics(theta,omega,J,n,T_d,T_c)

omega_dot = zeros(3,1);

omega_dot(1) = -(J(3,3)-J(2,2))/J(1,1)*omega(3)*omega(2) + 1/J(1,1)*T_d(1) + 3*n^2/J(1,1)*(J(3,3)-J(2,2))*sin(theta(1))*cos(theta(1))*(cos(theta(2)))^2 + 1/J(1,1)*T_c(1);
omega_dot(2) = -(J(1,1)-J(3,3))/J(2,2)*omega(1)*omega(3) + 1/J(2,2)*T_d(2) + 3*n^2/J(2,2)*(J(3,3)-J(1,1))*sin(theta(2))*cos(theta(1))*cos(theta(2)) + 1/J(2,2)*T_c(2);
omega_dot(3) = -(J(2,2)-J(1,1))/J(3,3)*omega(1)*omega(2) + 1/J(3,3)*T_d(3) + 3*n^2/J(3,3)*(J(1,1)-J(2,2))*sin(theta(1))*sin(theta(2))*cos(theta(2)) + 1/J(3,3)*T_c(3);
end