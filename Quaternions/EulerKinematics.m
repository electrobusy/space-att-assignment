% Complete Kinematics (- no matrix notation)
function theta_dot = EulerKinematics(theta,omega,n)

theta_dot = zeros(3,1);

theta_dot(1) = omega(1) + sin(theta(1))*tan(theta(2))*omega(2) + cos(theta(1))*tan(theta(2))*omega(3) + n*sin(theta(3))/cos(theta(2));
theta_dot(2) = cos(theta(1))*omega(2) - sin(theta(1))*omega(3) + n*cos(theta(3));
theta_dot(3) = sin(theta(1))/cos(theta(2))*omega(2) + cos(theta(1))/cos(theta(2))*omega(3) + n*tan(theta(2))*sin(theta(3));

end