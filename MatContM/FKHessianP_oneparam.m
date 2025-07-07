function y=FKHessianP_oneparam(x0,y0,tau,u)
%MatCont wants the HessianP function to give a 3d array of parital derivatives that are parameter derivatives of the Jacobian. 
%so H(i,j,k) is the kth parameter derivatives of the ij-th entry in the
%Jacobian. 

%These are mixed 2nd order partials using the mixed 2nd order central
%difference approximation.
h=10^(-6);
%%%%%First do flow-kick iterates at perturbations of the parameters.
%x or y and tau
xmhtaumh=flowkick(x0-h,y0,tau-h,u);
xphtauph=flowkick(x0+h,y0,tau+h,u);
xphtaumh=flowkick(x0+h,y0,tau-h,u);
xmhtauph=flowkick(x0-h,y0,tau+h,u);

ymhtaumh=flowkick(x0,y0-h,tau-h,u);
yphtauph=flowkick(x0,y0+h,tau+h,u);
yphtaumh=flowkick(x0,y0+h,tau-h,u);
ymhtauph=flowkick(x0,y0-h,tau+h,u);

%x or y and u
xmhamh=flowkick(x0-h,y0,tau,u-h);
xphaph=flowkick(x0+h,y0,tau,u+h);
xphamh=flowkick(x0+h,y0,tau,u-h);
xmhaph=flowkick(x0-h,y0,tau,u+h);

ymhamh=flowkick(x0,y0-h,tau,u-h);
yphaph=flowkick(x0,y0+h,tau,u+h);
yphamh=flowkick(x0,y0+h,tau,u-h);
ymhaph=flowkick(x0,y0-h,tau,u+h);

%old notes, I think:
%Jacobian [F_x F_y; G_x G_y];
%HessianP is 2x2x3 array: [D[J]/dk1; D[J]/dk2; D[J]/dtau ]

F_xtau=(xmhtaumh(1)+xphtauph(1)-xphtaumh(1)-xmhtauph(1))/(4*h^2);
F_ytau=(ymhtaumh(1)+yphtauph(1)-yphtaumh(1)-ymhtauph(1))/(4*h^2);
G_xtau=(xmhtaumh(2)+xphtauph(2)-xphtaumh(2)-xmhtauph(2))/(4*h^2);
G_ytau=(ymhtaumh(2)+yphtauph(2)-yphtaumh(2)-ymhtauph(2))/(4*h^2);

D_xa=(xmhamh+xphaph-xphamh-xmhaph)/(4*h^2);
D_ya=(ymhamh+yphaph-yphamh-ymhaph)/(4*h^2);


y(:,:,1)=[[F_xtau,F_ytau];[G_xtau,G_ytau]];
y(:,:,2)=[[D_xa(1),D_ya(1)];[D_xa(2),D_ya(2)]];
