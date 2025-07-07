%This script creates Figure 5 in the manuscript 
% "Recovery-induced tipping in Stommel’s kicked ocean box model".
%It calls stommel.m and Newton_cts.m

%% Branch of equilibria from A

%parameter
r=0;
r_incr=0.001;

%initial guess for Newton
x=[0.135;0.4836];

%tolerance for movement of points
tol=0.001;

%count
count=1;

%hold products
productsA=[];

while r_incr>0.0000001 && count<200
    newton=Newton_cts(x,r);
    x2=newton(:,end);
    if abs(x2-x)<tol
        productsA(1:2,count)=x;
        productsA(3,count)=r;
        x=x2;
        r=r+r_incr;
        count=count+1;
    else
        r_incr=r_incr/2;
        r=r-r_incr;
    end
end

%% Branch of equilibria from B

%parameter
r=0;
r_incr=0.001;

%initial guess for Newton
x=[0.3518;0.7651];

%tolerance for movement of points
tol=0.001;

%count
count=1;

%hold products
productsB=[];

while r_incr>0.0000001 && count<250
    newton=Newton_cts(x,r);
    x2=newton(:,end);
    if abs(x2-x)<tol
        productsB(1:2,count)=x;
        productsB(3,count)=r;
        x=x2;
        r=r+r_incr;
        count=count+1;
    else
        r_incr=r_incr/2;
        r=r-r_incr;
    end
end

%% Branch of equilibria from C

%parameter
r=0;
r_incr=0.001;

%initial guess for Newton
x=[0.4321;0.8203];

%tolerance for movement of points
tol=0.001;

%count
count=1;

%hold products
productsC=[];

while r_incr>0.0000001 && r<1
    newton=Newton_cts(x,r);
    x2=newton(:,end);
    if abs(x2-x)<tol
        productsC(1:2,count)=x;
        productsC(3,count)=r;
        x=x2;
        r=r+r_incr;
        count=count+1;
    else
        r_incr=r_incr/2;
        r=r-r_incr;
    end
end

%% Fig 5a (y vs. r)
figure
hold on
plot(productsA(3,:),productsA(2,:),'k-','LineWidth',2)
plot(productsB(3,:),productsB(2,:),'k--','LineWidth',2)
plot(productsC(3,:),productsC(2,:),'k-','LineWidth',2)
axis([0 0.05 0 1])

%% Fig 5b (y vs. x)
figure
hold on
plot(productsA(1,:),productsA(2,:),'k-','LineWidth',2)
plot(productsB(1,:),productsB(2,:),'k:','LineWidth',3)
plot(productsC(1,:),productsC(2,:),'k-','LineWidth',2)
%separatrix
[~,X]=ode45(@stommel,[0 5],[0.422 1]);
plot(X(:,1),X(:,2),'k--',LineWidth=1.5)
[~,X]=ode45(@stommel,[0 3],[1 0.794]);
plot(X(:,1),X(:,2),'k--',LineWidth=1.5)
%circulation=0 line
plot([0 0.5],[0 1],'Color',[0.5 0.5 0.5])
%special values of r
r=0.0277;
newton=Newton_cts([0.22781;0.592928],r);
x=newton(:,end);
plot(x(1),x(2),'c.')
%
r=0.0277;
newton=Newton_cts([0.414;0.768],r);
x=newton(:,end);
plot(x(1),x(2),'c.')
%
r=0.421;
newton=Newton_cts([0.3503;0.398],r);
x=newton(:,end);
plot(x(1),x(2),'c.')

axis([0 1 0 1])