function out = Stommel
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(Stommel);
y0=[0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,tau,u)
    dydt=transpose(flowkick(kmrgd(1),kmrgd(2),tau,u));
% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,tau,u)
    jac=FKJacobian_oneparam(kmrgd(1),kmrgd(2),tau,u);
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,tau,u)
    jacp=FKJacobianP_oneparam(kmrgd(1),kmrgd(2),tau,u);
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,tau,u)
    hess=FKHessian_oneparam(kmrgd(1),kmrgd(2),tau,u);
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,tau,u)
    hessp=FKHessianP_oneparam(kmrgd(1),kmrgd(2),tau,u);
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,tau,u)
    tens3=FKDer3_oneparam(kmrgd(1),kmrgd(2),tau,u);
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,tau,u)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,tau,u)
