function plotResult(t, s, theta, newFigs)
% plotResult - Plot the orbit of the satilite
% t - time vector
% t - time vector

r = s(:,1); u = s(:,2); v = s(:,3); 
if ~exist('theta','var') || isempty(theta)
  theta = atan2(-s(:,5), -s(:,6));
end
if ~exist('newFigs','var') 
  newFigs = true;
end
phi = cumtrapz(t,v./r); % angle to the shuttle in a polar coord system
x = r.*cos(phi); y = r.*sin(phi); 
n = length(r);
psi = pi/2*ones(n,1) + phi - theta; % angle of the thrust relative a global ref system
%%% Plot states
if newFigs
  figure;
else
  figure(1);
end
subplot(2,1,1), hold on
plot(t,r,'b.-',t,u,'r.-',t,v,'g.-',"LineWidth",1.2); grid
legend('r','u','v');
%%% Plot 
subplot(2,1,2), hold on
plot(t,mod(theta,2*pi),'b.-',"LineWidth",1.2); grid; 
ylabel('\theta (deg)')
xlabel('time')
%%% Plot thrust directions
if newFigs
  figure; 
else
  figure(2);
end
hold on
tt = linspace(t(1),t(end), 20);
x2 = interp1(t,x,tt);
y2 = interp1(t,y,tt);
psi2 = interp1(t,mod(psi,2*pi),tt);
%xt = x + 0.2*cos(psi);
%yt = y + 0.2*sin(psi); 
%plot([x,xt]',[y,yt]','r-');
xt2 = x2 + 0.2*cos(psi2);
yt2 = y2 + 0.2*sin(psi2); 
plot([x2;xt2],[y2;yt2],'r-');
%%% Plot Earth/Mars orbits etc.
plot(x,y,x(n),y(n),'go',0,0,'go'); grid; hold on 
q = (0:90-1)*pi/90; 
plot(cos(q),sin(q),'g--',r(end)*cos(q),r(end)*sin(q),'g--',1,0,'go');
text(-.6,.65,'Earth orbit');
text(-.1,.1,'Sun') 
text(.9,1.3,'Max orbit')
axis([-1.6 1.6 -.5 2.2]); axis equal
ylabel('y/r_e'), xlabel('x/r_e');