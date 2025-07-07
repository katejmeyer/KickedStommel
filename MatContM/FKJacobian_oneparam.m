function yprime=FKJacobian_oneparam(x0,y0,tau,u)
%This function is going to compute the jacobian of the flow-kick map. 
%we need to find the x and y partial derivatives of the time integral of
%the underlying DE.
%We'll compute these derivatives using central differences.
h=10^(-4);
%%%%%First calculate flow-kick map at points close by by calling the
%%%%%flowkick function file. output y is a 1x2 vector with 2d solution,
%%%%%inputs are (t,x-IC,y-IC,tau,u)
xph=flowkick(x0+h,y0,tau,u); %flow-kick with IC (x0+h,y0)
xmh=flowkick(x0-h,y0,tau,u); %flow-kick with IC (x0-h,y0)
yph=flowkick(x0,y0+h,tau,u); %flow-kick with IC (x0,y0+h)
ymh=flowkick(x0,y0-h,tau,u); %flow-kick with IC (x0,y0-h)

%F is x-coord of flow-kick map, G is y-coord of flow-kick map
%Jacobian is [DFdx DFdy; DGdx DGdy]
DFdx=(xph(1)-xmh(1))/(2*h); %DFdx is the change in the x-cood of flow-kick map with respect to x: (F(x+h,y)-F(x-h,y))/2h
DFdy=(yph(1)-ymh(1))/(2*h); %F(x,y+h)-F(x,y-h)/2h
DGdx=(xph(2)-xmh(2))/(2*h); %G(x+h,y)-G(x-h),y)/2h
DGdy=(yph(2)-ymh(2))/(2*h); %G(x,y+h)-G(x,y-h)/2h

yprime=[DFdx DFdy; DGdx DGdy];