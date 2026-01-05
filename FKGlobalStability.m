%This script is intended to test the presence of a globally-attracting
%flow-kick fixed point
% in the manuscript 
% "Recovery-induced tipping in Stommel’s kicked ocean box model"

%It returns an array of zeros and/or ones, with rows corresponding to tau
%values and columns corresponding to kick values. An entry of 0 means that
%flow-kick trajectories from a grid of different initial conditions all
%converged to within errorTol of each other after n iterations. An entry of
%1 means that they did not.

n=100; %number of iterations of flow-kick
errorTol=0.01; %error tolerance

%grid of initial conditions
xrange=0:0.2:1;
yrange=0:0.2:1;

%tau and kick values to test
tauRange=0.1:0.1:10;
kRange=0:0.1:1;
%array to store outcomes (0: trajectories converge to nearly same point; 1:
%they do not)
homogEndpoints=ones(length(tauRange),length(kRange));

tauCount=1;

for tau=tauRange
    kCount=1;
    for k=kRange
        FinalPositions=zeros(length(xrange),length(yrange),2); %store final positions from different initial conditions        
        xindex=1;
        for x=xrange
            xhold=x; %save for reset in y loop
            yindex=1;
            for y=yrange
                trajectory=zeros(2,n+1);
                trajectory(1,1)=x;
                trajectory(2,1)=y;
                % flow-kick trajectory               
                for i=1:n % flow-kick trajectory
                    [~,xf]=ode45(@stommel,[0,tau],[x y]);
                    x=xf(end,1)+k;
                    y=xf(end,2);
                end
                %store flow-kick trajectory endpoint
                FinalPositions(xindex,yindex,1)=x;
                FinalPositions(xindex,yindex,2)=y;
                yindex=yindex+1;
                x=xhold; %reset x
            end
            xindex=xindex+1;
        end

        if max(FinalPositions(:,:,1),[],"all")-min(FinalPositions(:,:,1),[],"all")<errorTol...
                && max(FinalPositions(:,:,2),[],"all")-min(FinalPositions(:,:,2),[],"all")<errorTol
            homogEndpoints(tauCount,kCount)=0;
        end
        kCount=kCount+1;
    end
    tauCount=tauCount+1;
end

homogEndpoints
