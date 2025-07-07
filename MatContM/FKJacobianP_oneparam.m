function y=FKJacobianP_oneparam(x0,y0,tau,u)
%This is function file to calculate first derivatives of the of the
%flow-kick map with respect to the parameters using central differences. 
h=10^(-6);
%%%%%First flow-kick iterates at perturbations of the parameters.
tauph=flowkick(x0,y0,tau+h,u); 
taumh=flowkick(x0,y0,tau-h,u); 
uph=flowkick(x0,y0,tau,u+h);
umh=flowkick(x0,y0,tau,u-h);

%Derivatives of F (x-coord of flowkick) and G (y-coord of flowkick)
DFdtau=(tauph(1)-taumh(1))/(2*h); %F(tau+h)-F(tau-h)/2h
DGdtau=(tauph(2)-taumh(2))/(2*h); %G(tau+h)-G(tau-h))/2h
Da=(uph-umh)/(2*h);

y=[DFdtau, Da(1); DGdtau, Da(2)];