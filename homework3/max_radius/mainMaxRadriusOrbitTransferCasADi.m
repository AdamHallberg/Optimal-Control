% Max radius orbit transfer in a given time
% ----------------------
% An optimal control problem (OCP),

import casadi.*

N = 100; % number of control intervals

%% Model Parameters

TT = 0.1405;
mdot = 0.07489;
tf=3.3155;

%% Optimization Problem

opti = casadi.Opti(); % Optimization problem

% ---- decision variables ---------
X = opti.variable(3,N+1); % state trajectory
r = X(1,:);
u = X(2,:);
v = X(3,:);
Theta = opti.variable(1,N);   % control trajectory (throttle)

% ---- objective          ---------
opti.minimize(???); % -r(tf)

% ---- dynamic constraints --------
f = @(t,x,u)[???]; % dx/dt = f(x,u)

dt = tf/N; % length of a control interval
for k=1:N % loop over control intervals
   % Runge-Kutta 4 integration
   t(k)=(k-1)*dt;
   k1 = f((k-1)*dt,X(:,k),         Theta(:,k));
   k2 = f((k-1)*dt+0.5*dt,X(:,k)+dt/2*k1, Theta(:,k));
   k3 = f((k-1)*dt+0.5*dt, X(:,k)+dt/2*k2, Theta(:,k));
   k4 = f((k-1)*dt+dt, X(:,k)+dt*k3,   Theta(:,k));
   x_next = X(:,k) + dt/6*(k1+2*k2+2*k3+k4); 
   opti.subject_to(X(:,k+1)==x_next); % close the gaps
end

% ---- boundary conditions -------

opti.subject_to(r(1)==???);   % Initial Condition
opti.subject_to(u(1)==???);   % Initial Condition
opti.subject_to(v(1)==???);   % Initial Condition
opti.subject_to(u(N+1)==???); % final Constraint
opti.subject_to(???); % final Constraint

% ---- initial values for solver ---
opti.set_initial(???, ???); % Hint: The initial value for the solver is zero. However, in the system's equation, you have 1/r. Therefore, please adjust the initial conditions of the solver for 'r' to a different value for all sampling times.

% ---- solve NLP              ------
opti.solver('ipopt'); % set numerical backend
sol = opti.solve();   % actual solve

% ---- post-processing        ------

s=[sol.value(r(1:N));sol.value(u(1:N));sol.value(v(1:N))]';
plotResult(t(1:N), s,sol.value(Theta(1:N))');
r_tf=sol.value(r(???))
