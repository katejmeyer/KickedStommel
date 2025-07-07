function [ Newton_pts, flag ] = Newton_2d(tau,k,x0) %put in 2d column for x0 and k
%NEWTON Implementation of Newton's Method to find zeros of F(x)=phi^tau(x)+k-x
%   These are flow-kick fixed points / equilibria
%   You will need the function CoupledVar_2d.m to run Newton_2d.m

% 5/20/25 Working 
% (have to input kick and initial position as COLUMN VECTORS!!!!)


%F_of_x=1;

%Set up an initial condition
x = x0;
ODEinit = [x(1); x(2); 1; 0; 0; 1];

%Set up a place to hold Newton iterates x_n 
Newton_pts=[];

%Set up a flag
flag=1; % will be set to zero if Newton converges, otherwise stays 1.

%Set up counter
cnt=0;

while cnt<1 &&  flag==1

cnt=cnt+1;
    
%Compute solution to variational equation along trajectory starting at
%x and lasting for time tau
[~,X] = ode45(@CoupledVar_2d,[0,tau],ODEinit);

post_flow=X(end,1:2)'; %this is phi^tau of x

Vsoln = [X(end,3), X(end,5); X(end,4), X(end,6)]; %Variational equation solution in matrix form

%Compute DF
DF=Vsoln - eye(2);

%Compute F(x)
F_of_x= post_flow + k - x;

%Compute next Newton iterate
x =  x - (inv(DF)*F_of_x);
%disp(inv(DF)*F_of_x); %for debugging

ODEinit = [x(1); x(2); 1; 0; 0; 1];


Newton_pts(1:2,cnt) = x;

if norm(F_of_x) < 10^(-10) %check whether x is nearly at equilibrium
    flag=0;
end

end

flag;

end