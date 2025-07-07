function y=flowkick(x0,y0,tau,k) 
% Stommel's two-box model with kicks
% k is the salinity kick parameter
% tau is flow time

%undisturbed system
f= @(~,x) [(1/6)*(1-x(1))-5*abs(2*x(1)-x(2))*x(1); 1-x(2)-5*abs(2*x(1)-x(2))*x(2)];

%flow
[~,X]=ode45(f,[0,tau],[x0,y0]); %compute flow for time tau
    
%kick radially with u as scaling factor.
y=X(end,:)+[k, 0];
end