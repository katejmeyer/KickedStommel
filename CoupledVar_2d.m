function [ dzdt ] = CoupledVar_2d(~,z)
% gives original ODE, along with varational equation to be solved
% simultaneously

dzdt = zeros(6,1); %set up structure to hold vector field functional output

R=2;
delta=1/6;
lambda=0.2;
vfield = [delta*(1-z(1))-(1/lambda)*abs(R*z(1)-z(2))*z(1),1-z(2)-(1/lambda)*abs(R*z(1)-z(2))*z(2)]; %defining the Stommel vector field

dzdt(1)= vfield(1); 
dzdt(2)= vfield(2); %together these solve Stommel ODE

%Jacobian of Stommel's vector field
Jac=[-delta-(1/lambda)*abs(R*z(1)-z(2))-(R*z(1)/lambda)*sign(R*z(1)-z(2)),(z(1)/lambda)*sign(R*z(1)-z(2));-(R*z(2)/lambda)*sign(R*z(1)-z(2)),-1-(1/lambda)*abs(R*z(1)-z(2))+(z(2)/lambda)*sign(R*z(1)-z(2))];

lhs= Jac*[z(3),z(5);z(4),z(6)];
dzdt(3)=lhs(1,1);
dzdt(4)=lhs(2,1);
dzdt(5)=lhs(1,2);
dzdt(6)=lhs(2,2);

end