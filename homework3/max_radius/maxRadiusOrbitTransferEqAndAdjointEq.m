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
ds(1) = ???;
ds(2) = ???;
ds(3) = ???;
ds(4) = ???;
ds(5) = ???;
ds(6) = ???;
