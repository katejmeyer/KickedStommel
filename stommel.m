function [dxdt] = stommel(~,x)
%STOMMEL vector field for Stommel's two-box model
%   Created to make a phase portrait
%   x(1)=x=salinity anomaly
%   x(2)=y=temperature anomaly
R=2;
delta=1/6;
lambda=0.2;
dxdt=zeros(2,1);
dxdt(1) = delta*(1-x(1))-(1/lambda)*abs(R*x(1)-x(2))*x(1) ;
dxdt(2)= 1-x(2)-(1/lambda)*abs(R*x(1)-x(2))*x(2);
end