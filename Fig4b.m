%This script creates Figure 4b in the manuscript 
% "Recovery-induced tipping in Stommel’s kicked ocean box model".
% It calls minusstommel.m and Newton_2d.m, which calls CoupledVar_2d.m
k=[0.1;0];
figure
hold on

%% Separatrix
[~,X]=ode45(@minusstommel,[0 4],[0.353,0.767]);
plot(X(:,1),X(:,2),'--','Color',[0.7 0.7 0.7],LineWidth=1.5)
[~,X]=ode45(@minusstommel,[0 4],[0.352,0.765]);
plot(X(:,1),X(:,2),'--','Color',[0.7 0.7 0.7],LineWidth=1.5)

%% A-branch of flow-kick fixed points
x0=[0.23;0.58];
tau=9;
tau_incr=0.1;

Abranch=[];

%get a really exact initial guess so it passes the if/else
newtonpts=Newton_2d(tau,k,x0);
x0=newtonpts(:,end);

for i=1:500
    newtonpts=Newton_2d(tau,k,x0);
    if abs(newtonpts(:,end)-x0)<0.01
        x0=newtonpts(:,end);
        tau=tau-tau_incr;
        Abranch=[Abranch,[x0;tau]];
    else
        tau_incr=tau_incr/2;
    end

end

plot(Abranch(1,:),Abranch(2,:),'b-',LineWidth=1.5)

%% B-branch of flow-kick fixed points
mypurple=[88/255 6/255 213/255];

x0=[0.29;0.55];
tau=9;
tau_incr=0.1;
Bbranch=[];

%get a really exact initial guess so it passes the if/else
newtonpts=Newton_2d(tau,k,x0);
x0=newtonpts(:,end);

for i=1:500
    newtonpts=Newton_2d(tau,k,x0);
    if abs(newtonpts(:,end)-x0)<0.009
        x0=newtonpts(:,end);
        tau=tau-tau_incr;
        Bbranch=[Bbranch,[x0;tau]];
    else
        tau_incr=tau_incr/2;
    end

end

plot(Bbranch(1,:),Bbranch(2,:),'Color',mypurple,LineWidth=1.5)

%% C-branch of flow-kick fixed points
myorange=[1 102/255 0];

x0=[0.53;0.92];
tau=9;
tau_incr=0.01;
Cbranch=[];

%get a really exact initial guess so it passes the if/else
newtonpts=Newton_2d(tau,k,x0);
x0=newtonpts(:,end);

for i=1:891
    newtonpts=Newton_2d(tau,k,x0);
    x0=newtonpts(:,end);
    Cbranch=[Cbranch,[x0;tau]];
    tau=tau-tau_incr;
end

plot(Cbranch(1,:),Cbranch(2,:),'Color',myorange,'LineWidth',1.5)

%% Additional plot features

axis([0 1 0 1])

mylightorange=[1 102/255 0 0.5];

%plot initial flow-kick cycle on C branch
[~,X]=ode45(@CoupledVar_2d,[0 9],[0.532062898865720 0.820344611402048 1 0 0 1]);
plot(X(:,1),X(:,2),'Color',mylightorange,'LineWidth',1.5)

%plot flow-kick cycle at tau=1 on C branch
[~,X]=ode45(@CoupledVar_2d,[0 1],[0.469957166272185 0.714812784995343 1 0 0 1]);
plot(X(:,1),X(:,2),'Color',mylightorange,'LineWidth',1.5)

%plot flow-kick cycle at separatrix on C branch
[~,X]=ode45(@CoupledVar_2d,[0 0.233],[0.4103 0.4058 1 0 0 1]);
plot(X(:,1),X(:,2),'Color',mylightorange,'LineWidth',1.5)

%plot flow-kick cyle at tau=0.1 on C branch
[~,X]=ode45(@CoupledVar_2d,[0 0.1],[0.460557117680800 0.272874204896324 1 0 0 1]);
plot(X(:,1),X(:,2),'Color',mylightorange,'LineWidth',1.5)
