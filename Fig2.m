%This script creates a elements of a phase portrait for Stommel's 
%nondimensionalized two-box model of ocean circulation.
%It calls the function stommel.m.
figure
hold on

initCond=[0 0.15; 0.05 1; 0.23 1; 0.35 1; 0.422 1; 0.45 1; 0.7 1;  1 0.794; 1 0.35;  0.25 0; 0.05 0];
for i=1:size(initCond,1)
    [~,X]=ode45(@stommel,[0 100],initCond(i,:));
    plot(X(:,1),X(:,2),'k-',LineWidth=0.8)
    %X(end,:)
end
axis([0 1 0 1])

%equilibrium A
eqA=plot(0.1350,0.4835,'ko');
eqA.MarkerFaceColor="black";
%equilibrium B
eqB=plot(0.352,0.765,'ko');
eqB.MarkerFaceColor="white";
%equilibrium C
eqC=plot(0.4321,0.8203,'ko');
eqC.MarkerFaceColor="black";

%flow=0 line
flowline=plot([0 0.5],[0 1]);
flowline.Color=[0.5 0.5 0.5];
flowline.LineWidth=1.25;

