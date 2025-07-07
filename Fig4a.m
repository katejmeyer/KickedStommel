%This script creates Figure 4a in the manuscript 
% "Recovery-induced tipping in Stommel’s kicked ocean box model"
% It calls the functions stommel.m and Newton_2d.m, which calls
% CoupledVar_2d.m.

k=[0.1;0];
tau=9;

figure
hold on
axis([0 1 0 1])

%% Separatrix
[~,X]=ode45(@stommel,[0 2],[0.422 1]);
plot(X(:,1),X(:,2),'Color',[0.75 0.75 0.75],LineWidth=1,LineStyle='--')
[~,X]=ode45(@stommel,[0 2],[1 0.793]);
plot(X(:,1),X(:,2),'Color',[0.75 0.75 0.75],LineWidth=1,LineStyle='--')

%% Equilibria for Stommel's undisturbed model
plot([0.1350],[0.4836],'o',MarkerSize=5,MarkerFaceColor=[0.75 0.75 0.75],MarkerEdgeColor=[0.75 0.75 0.75]) %equilibrium A
plot([0.3518],[0.7651],'o',MarkerSize=5,MarkerFaceColor=[1 1 1],MarkerEdgeColor=[0.75 0.75 0.75]) %equilibrium B
plot([0.4320],[0.8203],'o',MarkerSize=5,MarkerFaceColor=[0.75 0.75 0.75],MarkerEdgeColor=[0.75 0.75 0.75]) %equilibrium C

%% A-like flow-kick fixed point
x0=[0.23;0.58]; %nearby guess

%find fixed point numerically
newtonpts=Newton_2d(tau,k,x0);
x0=newtonpts(:,end);

%plot fixed point
plot(x0(1),x0(2),'o',MarkerSize=5,MarkerFaceColor=[0 0 1],MarkerEdgeColor=[0 0 1])

%plot flow
[~,X]=ode45(@CoupledVar_2d,[0 tau],[x0' 1 0 0 1]);
plot(X(:,1),X(:,2),'Color',[0 0 1],LineWidth=1.5)

%plot kick
plot([X(end,1),X(1,1)],[X(end,2),X(1,2)],'Color',[0 0 1],LineWidth=1.5,LineStyle=':')

%% B-like flow-kick fixed point
x0=[0.29;0.55]; %nearby guess

%find fixed point numerically
newtonpts=Newton_2d(tau,k,x0);
x0=newtonpts(:,end);

%plot fixed point
plot(x0(1),x0(2),'o',MarkerSize=5,MarkerFaceColor=[1 1 1],MarkerEdgeColor=[0.33 0 0.83])

%plot flow
[~,X]=ode45(@CoupledVar_2d,[0 tau],[x0' 1 0 0 1]);
plot(X(:,1),X(:,2),'Color',[0.33 0 0.83],LineWidth=1.5)

%plot kick
plot([X(end,1),X(1,1)],[X(end,2),X(1,2)],'Color',[0.33 0 0.83],LineWidth=1.5,LineStyle=':')

%% C-like fixed point
x0=[0.53;0.92]; %nearby guess

%find fixed point numerically
newtonpts=Newton_2d(tau,k,x0);
x0=newtonpts(:,end);

%plot fixed point
plot(x0(1),x0(2),'o',MarkerSize=5,MarkerFaceColor=[1 0.4 0],MarkerEdgeColor=[1 0.4 0])

%plot flow
[~,X]=ode45(@CoupledVar_2d,[0 tau],[x0' 1 0 0 1]);
plot(X(:,1),X(:,2),'Color',[1 0.4 0],LineWidth=1.5)

%plot kick
plot([X(end,1),X(1,1)],[X(end,2),X(1,2)],'Color',[1 0.4 0],LineWidth=1.5,LineStyle=':')

