function yprime=FKHessian_oneparam(x0,y0,tau,u)
%This function is going to compute the hessian of the flow-kick map. 
%we need to find the x and y second partial derivatives of the time integral of
%the underlying DE.
%We'll compute these derivatives using central differences.
%f''=(f(x+h)-2f(x)+f(x-h))/h^2
h=10^(-4);
%%%%%First calculate flow-kick map at points close by by calling the
%%%%%flowkick function file. output y is a 1x2 vector with 2d solution,
%%%%%inputs are (x-IC,y-IC,kick x-coord, kick, y-coord,tau)
xy=flowkick(x0,y0,tau,u); %flow-kick with IC (x0,y0)
xph=flowkick(x0+h,y0,tau,u); %flow-kick with IC (x0+h,y0)
xmh=flowkick(x0-h,y0,tau,u); %flow-kick with IC (x0-h,y0)
yph=flowkick(x0,y0+h,tau,u); %flow-kick with IC (x0,y0+h)
ymh=flowkick(x0,y0-h,tau,u); %flow-kick with IC (x0,y0-h)
xmhymh=flowkick(x0-h,y0-h,tau,u);
xphyph=flowkick(x0+h,y0+h,tau,u);
xphymh=flowkick(x0+h,y0-h,tau,u);
xmhyph=flowkick(x0-h,y0+h,tau,u);

%F is x-coord of flow-kick map, G is y-coord of flow-kick map
%Hessian for the 2d flow kick is going to be a 2x2x2 array [ [[F_xx
%F_yx],[G_xx G_yx]]; [[F_xy F_yy],[G_xy G_yy]] ]
F_xx=(xph(1)-2*xy(1)+xmh(1))/h^2;
F_yx=(xmhymh(1)+xphyph(1)-xphymh(1)-xmhyph(1))/(4*h^2);
G_xx=(xph(2)-2*xy(2)+xmh(2))/h^2;
G_yx=(xmhymh(2)+xphyph(2)-xphymh(2)-xmhyph(2))/(4*h^2);
F_xy=F_yx;
F_yy=(yph(1)-2*xy(1)+ymh(1))/h^2;
G_xy=G_yx;
G_yy=(yph(2)-2*xy(2)+ymh(2))/h^2;
yprime(:,:,1)=[[F_xx, F_yx];[G_xx, G_yx]];
yprime(:,:,2)= [[F_xy, F_yy];[G_xy G_yy]];
%yprime=[ [[F_xx, F_yx];[G_xx, G_yx]]; [[F_xy, F_yy];[G_xy G_yy]] ];