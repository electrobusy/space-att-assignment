% Complete Kinematics (matrix notation)
function theta_dot = EulerKinematicsMat(theta,omega,n)

% - Body angular velocity term:
N = [
    1 sin(theta(1))*tan(theta(2)) cos(theta(1))*tan(theta(2));
    0 cos(theta(1)) -sin(theta(1));
    0 sin(theta(1))/cos(theta(2)) cos(theta(1))/cos(theta(2));
    ];
f11 = N*omega;

% - Gravity gradient term:
f12 = n*[
    sin(theta(3))/cos(theta(2));
    cos(theta(3));
    tan(theta(2))*sin(theta(3));
    ];

% - Sum of previous terms gives total angular velocity of the Spacecraft:
theta_dot = f11 + f12;

end