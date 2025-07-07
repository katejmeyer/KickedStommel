function [ Newton_pts, flag ] = Newton_cts(x0,r) %put in 2d column for x0
%NEWTON Implementation of Newton's Method to find zeros of Stommel+r vector
%field
%   These are ODE fixed points / equilibria

%Parameters
R=2;
delta=1/6;
lambda=0.2;

%Set up an initial condition
x = x0;

%Set up a place to hold Newton iterates x_n 
Newton_pts=[];

%Set up a flag
flag=1; % will be set to zero if Newton converges, otherwise stays 1.

%Set up counter
cnt=0;

while cnt<100 &&  flag==1
    
    cnt=cnt+1;
        
    %Compute Df
    DF=[-delta-(1/lambda)*abs(R*x(1)-x(2))-(R*x(1)/lambda)*sign(R*x(1)-x(2)),...
        (x(1)/lambda)*sign(R*x(1)-x(2));...
        -(R*x(2)/lambda)*sign(R*x(1)-x(2)),...
        -1-(1/lambda)*abs(R*x(1)-x(2))+(x(2)/lambda)*sign(R*x(1)-x(2))];
    
    %Compute F(x)
    F_of_x=[delta*(1-x(1))-(1/lambda)*abs(R*x(1)-x(2))*x(1)+r;...
        1-x(2)-(1/lambda)*abs(R*x(1)-x(2))*x(2)];
    
    %Compute next Newton iterate
    x =  x - (inv(DF)*F_of_x);
    %disp(inv(DF)*F_of_x) %for debugging
    
    Newton_pts(1:2,cnt) = x;
    
    if norm(F_of_x) < 10^(-10) %check whether x is nearly at equilibrium
        flag=0;
    end

end
endpoint=Newton_pts(end,1);

end