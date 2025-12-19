function ds = maxRadiusOrbitTransferEqAndAdjointEq(t, s)
global T mdot
a = T/(1-mdot*t);

r = s(1);
u = s(2);
v = s(3);
lambda_r = s(4);
lambda_u = s(5);
lambda_v = s(6);

% Choose u such as Hu = 0
theta = atan2(-lambda_u, -lambda_v);

ds = zeros(6,1);

% Motivation in report.
ds(1) = u;
ds(2) = v^2/r - 1/r^2 + a*sin(theta);
ds(3) = -u*v/r + a*cos(theta);
ds(4) = lambda_u*(v^2/(r^2) - 2/(r^3)) - lambda_v*(u*v/(r^2));
ds(5) = -lambda_r + lambda_v*v/r;
ds(6) = -lambda_u*2*v/r + lambda_v*u/r;
