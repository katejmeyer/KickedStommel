%This script can be used to make panels b-e of Figure 3 in 
% "Recovery-induced tipping in Stommel’s kicked ocean box model"
% by modifying tau and/or kick below.

tau=0.1; % 0.1 for panel b, 1 for panel c, 9 for panel d, 5 for panel e
kick=0.1; % 0.1 for panels b,c,d; 0.4 for panel e

figure
hold on

%separatrix
[~,X]=ode45(@stommel,[0,5.5],[0.422 1]);
plot(X(:,1),X(:,2),'k--','LineWidth',1)
[~,X]=ode45(@stommel,[0,6],[1 0.794]);
plot(X(:,1),X(:,2),'k--','LineWidth',1)

%flow=0 line
plot([0 0.5],[0 1],'color',[0,0,0]+0.5,'LineWidth',1)

x=[0.135 0.4835]; %equilibrium A
plot(x(1),x(2),'k.')
for i=1:50 % flow-kick trajectory
    [~,xf]=ode45(@stommel,[0,tau],x);
    plot(xf(:,1),xf(:,2),'k-','LineWidth',1.5)
    x=xf(end,:)+[kick 0];
    plot([x(1), xf(end,1)],[x(2), xf(end,2)],'k:','LineWidth',1.5)
end

%final flow-kick of trajectory in red
[~,xf]=ode45(@stommel,[0,tau],x);
plot(xf(:,1),xf(:,2),'r-','LineWidth',2)
x=xf(end,:)+[kick 0];
plot([x(1), xf(end,1)],[x(2), xf(end,2)],'r:','LineWidth',2)

axis([0 1 0 1])
pbaspect([1 1 1])