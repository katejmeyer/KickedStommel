function y=FKDer3_oneparam(x0,y0,tau,u)
%This function is going to compute the third order derivatives of the
%flow-kick map.
%We'll compute these derivatives using central differences.
%
h=10^(-6);
%%%%%First calculate flow-kick map at points close by by calling the
%%%%%flowkick function file. output y is a 1x2 vector with 2d solution,
%%%%%inputs are (x-IC,y-IC,kick in x-coord, kick in y-coord,tau)
xy=flowkick(x0,y0,tau,u); %flow-kick with IC (x0,y0)
xph=flowkick(x0+h,y0,tau,u); %flow-kick with IC (x0+h,y0)
xmh=flowkick(x0-h,y0,tau,u); %flow-kick with IC (x0-h,y0)
yph=flowkick(x0,y0+h,tau,u); %flow-kick with IC (x0,y0+h)
ymh=flowkick(x0,y0-h,tau,u); %flow-kick with IC (x0,y0-h)
xp2h=flowkick(x0+2*h,y0,tau,u);
xm2h=flowkick(x0-2*h,y0,tau,u);
yp2h=flowkick(x0,y0+2*h,tau,u);
ym2h=flowkick(x0,y0-2*h,tau,u);

xmhymh=flowkick(x0-h,y0-h,tau,u);
xphyph=flowkick(x0+h,y0+h,tau,u);
xphymh=flowkick(x0+h,y0-h,tau,u);
xmhyph=flowkick(x0-h,y0+h,tau,u);

%F is x-coord of flow-kick map, G is y-coord of flow-kick map
%Hessian for the 2d flow kick is going to be a 2x2x2 array [ [[F_xx
%F_yx],[G_xx G_yx]]; [[F_xy F_yy],[G_xy G_yy]] ]
F_xxx=(-xm2h(1)+2*xmh(1)-2*xph(1)+xp2h(1))/(2*h^3);
F_xxy=(xmhyph(1)-2*yph(1)+xphyph(1)-xmhymh(1)+2*ymh(1)-xphymh(1))/(2*h^3);
F_xyy=(-xmhyph(1)+xphyph(1)+2*xmh(1)-2*xph(1)-xmhymh(1)+xphymh(1))/(2*h^3);
F_yyy=(-ym2h(1)+2*ymh(1)-2*yph(1)+yp2h(1))/(2*h^3);

G_xxx=(-xm2h(2)+2*xmh(2)-2*xph(2)+xp2h(2))/(2*h^3);
G_xxy=(xmhyph(2)-2*yph(2)+xphyph(2)-xmhymh(2)+2*ymh(2)-xphymh(2))/(2*h^3);
G_xyy=(-xmhyph(2)+xphyph(2)+2*xmh(2)-2*xph(2)-xmhymh(2)+xphymh(2))/(2*h^3);
G_yyy=(-ym2h(2)+2*ymh(2)-2*yph(2)+yp2h(2))/(2*h^3);

y(:,:,1,1)=[[F_xxx,F_xxy];[G_xxx,G_xxy]];
y(:,:,1,2)=[[F_xxy,F_xyy];[G_xxy,G_xyy]];
y(:,:,2,1)=[[F_xxy,F_xyy];[G_xxy,G_xyy]];
y(:,:,2,2)=[[F_xyy,F_yyy];[G_xyy,G_yyy]];