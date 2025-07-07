%This script creates Figure 3a in the manuscript 
% "Recovery-induced tipping in Stommel’s kicked ocean box model"

%create grid of tau and kappa values
krange=0:0.005:1;
trange=0.05:0.05:10;
%create arrays to hold basin and circulation outcomes for each (tau, kick)
basinArray=zeros(length(krange),length(trange));
circArray=zeros(length(krange),length(trange));

kindex=1;
for kick=krange
tindex=1;    
for tau=trange    
    x=[0.135 0.4835]; %start at equilibrium A
    % flow-kick trajectory from A
    for i=1:100 
        [~,xf]=ode45(@stommel,[0,tau],x);
        x=xf(end,:)+[kick 0];
    end
    assignments = zeros(length(xf),1); %classify circulation direction along flow
    for i = 1:length(xf)
        if xf(i,2)<2*xf(i,1)
            assignments(i) = 0;
        elseif 1 == 1
            assignments(i) = 1;
        end
    end
    %check circulation
    if x(end,2)<2*x(end,1) && xf(end,2)<2*xf(end,1)
        circArray(kindex,tindex)=3; %reversed circulation
    elseif x(end,2)<2*x(end,1) && xf(end,2)>2*xf(end,1)
        circArray(kindex,tindex)=5; %mixed circulation
    elseif x(end,2)>2*x(end,1) && xf(end,2)<2*xf(end,1)
        circArray(kindex,tindex)=7; %mixed circulation
    elseif x(end,2)>2*x(end,1) && xf(end,2)>2*xf(end,1)
        circArray(kindex,tindex)=11; %preserved circulation
    end

    % Assign mixed circulation even when start and end circulations are the same, but the
    % flow travels through a region with different circulation 
    if sum(assignments) ~= length(assignments) && sum(assignments) ~= 0
        circArray(kindex, tindex) = 5;
    end
    
    %check basin
    [~,Y]=ode45(@stommel,[0, 200],x(end,:)); %run flow forward to check long-term behavior
    if abs(Y(end,:)-[0.135,0.4835])<0.01 %check if trajectory converges to near A
        basinArray(kindex,tindex)=1;
    elseif abs(Y(end,:)-[0.4321,0.8202])<0.01 %check if trajectory converges to near C
        basinArray(kindex,tindex)=2;
    end    
    tindex=tindex+1;
end
    kindex=kindex+1;
end

bothArray=basinArray.*circArray; %multiply prime elements from each array to create unique identifying number
bothArray=flipud(bothArray);

%% Process unique identifiers into matrices of 1s and 0s to plot
bothArray3=bothArray==3;
bothArray5=bothArray==5;
bothArray7=bothArray==7;
bothArray11=bothArray==11;
bothArray6=bothArray==6;
bothArray10=bothArray==10;
bothArray14=bothArray==14;
bothArray22=bothArray==22;

[r3,c3]=find(bothArray3);
[r5,c5]=find(bothArray5);
[r7,c7]=find(bothArray7);
[r11,c11]=find(bothArray11);
[r6,c6]=find(bothArray6);
[r10,c10]=find(bothArray10);
[r14,c14]=find(bothArray14);
[r22,c22]=find(bothArray22);

%%
figure
hold on
scatter(c3, r3, 40, 'sb', 'filled')
scatter(c5, r5, 40, 'sg', 'filled')
scatter(c7, r7, 40, 'sk', 'filled')
scatter(c11, r11, 40, 'sc', 'filled')
scatter(c6, r6, 40, 'sw', 'filled')
scatter(c10, r10, 40, 'sy', 'filled')
scatter(c14, r14, 40, 'sr', 'filled')
scatter(c22, r22, 40, 'sr', 'filled')
set(gca, 'Color','k', 'YDir','reverse')
axis([0 size(basinArray,2) 0 size(basinArray,1)])

%% lines to compare to continuous disturbance
plot([0, 200],[200, 200-(0.028*200*10)],'k-')
plot([0, 200],[200, 200-(0.421*200*10)],'k-')

%% specific parameter combinations for figure panels
plot(2,180,'ko') %b
plot(20,180,'ko') %c
plot(180,180,'ko') %d
plot(40,100,'ko') %e

%% rearrange layers for clarity of figure
h=get(gca,'Children');
uistack(h(11),'up',1);
uistack(h(13),'up',3);
uistack(h(9),'down',2);
uistack(h(14),'up',4);

