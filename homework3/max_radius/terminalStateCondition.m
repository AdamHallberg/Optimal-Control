function mu = terminalStateCondition(lambda0, x0, tf)
% terminalStateCondition - Computes deviation from the final constrints for the joint system, (10.4) in the lecture notes
% lambda0 - Adjoint variable at time t=0
% x0      - State at time t=0
% tf      - Final time

% Integrate the system and adjoint equations forward in time
[t,s] = ode23(@maxRadiusOrbitTransferEqAndAdjointEq, [0 tf], [x0;lambda0]); 
r_tf = s(end,1);
u_tf = s(end,2);
v_tf = s(end,3);
lambda_r_tf = s(end, 4);
lambda_u_tf = s(end, 5);
lambda_v_tf = s(end, 6);
mu = zeros(3,1);
mu(1) = ???; % final state constraint 1
mu(2) = ???; % final state constraint 2 
mu(3) = ???; % final state constraint 3

