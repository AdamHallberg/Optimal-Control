% Maximum orbit radius transfer problem solved with a shooting method.
global T mdot
T = 0.1405;
mdot = 0.07489;
tf = 3.3155; % final time
xi = [1,0,1]'; % initial constraint
%%% Optimization options
optim = optimset('fsolve');
optim.Display = 'off';
optim.Algorithm = 'levenberg-marquardt';

%%% Start guess
lambda0 = [???; ???; ???];



%%% Solve the problem
lambda0 = fsolve(@terminalStateCondition, lambda0, optim, xi, tf); % finds a lambda0 such that terminalStateCondition(lambda0) = 0

%%% Check the terminal constraints are satisfied
maxConditionError=1e-3;
if norm( terminalStateCondition(lambda0, xi, tf) ) > maxConditionError
  disp('Warning! Terminal constraints not satisfied!'); 
end

[t,s] = ode23(@maxRadiusOrbitTransferEqAndAdjointEq, [0 tf], [xi;lambda0]); 

finalRadius = ???
fprintf('Final radius: %f \n', finalRadius)
plotResult(t, s);