%% 0. continuation on Stommel flow-kick model
    %adapt code from Alanna Hoyer-Leitzel to Stommel

%% 1. initialize stuff

%matcont directory
cd ~/Documents/MATLAB/matcontm5p4copy/

%use the right flowkick system
status = copyfile('Systems/flowkick_Stommel.m', 'Systems/flowkick.m');

%initialize matcont
init
global opt cds fpmds

%set systems file to be used by matcont
sval=Stommel; %it contains the map and its derivatives as functions

%run will be struct - clear it so it'll work if it is previously defined.
clear('run')



%% 2. set up continuation parameters
% params are in order [tau (flow time); k (salinity kick)] 
contparam = 1; % 1 for first parameter, 2 for second parameter
params=[0.41;0.01]; %initial parameter values

opt = contset;

%% 3. want to find a stable fixed point. For tau = 5.5 and k=0.1, there should be one at (0.2452,0.4969) 
x1 = [0.2 0.45]; %close to expected stable flow-kick fixed point

fun_eval = sval{2}; %the second function in the sval structure is the flow-kick map
% Iterate the map and hope to converge to a fixed point
for j = 1:100
    x1=feval(fun_eval,0,x1,params(1),params(2)); % the second input is t
end
%can leave x1= ... without ; to check progress or suppress output with ;


%% 4. continue stable fixed point in tau

n=1; %sets one iteration of flow-kick map rather than two-iterations if n=2, etc.

[x0,v0]=init_FPm_FPm(@Stommel,x1, params, contparam, n); %init_FPm_FPm is a built-in function in the FixedPointMap folder

%indices keep track of which direction the continuation runs
tauindex=1 %print to keep track of output
opt=contset(opt,'MaxNumPoints',100);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Multipliers',1);
opt=contset(opt,'MaxStepSize', 0.1);
opt=contset(opt,'Backward',0); %run forward
[run(tauindex).x,run(tauindex).v,run(tauindex).s,run(tauindex).h,run(tauindex).f]=cont(@fixedpointmap,x0,v0,opt); %cont is built-in function in Continuer folder

%run backward
tauindex=2
opt=contset(opt,'Backward',1); 
opt=contset(opt,'MaxNumPoints',100);
[run(tauindex).x,run(tauindex).v,run(tauindex).s,run(tauindex).h,run(tauindex).f]=cont(@fixedpointmap,x0,v0,opt);

% (NEXT, CONTINUE STABLE FIXED POINT A'; REPEAT FOR B' AND C')

%% 5. Plot branches of continued fixed points 
% have to use View > Rotate3D tool in Figure window to see in 3D

figure
hold on
plot3(params(1),x1(1),x1(2),'ro') %starting point
plot3(run(1).x(3,:),run(1).x(1,:),run(1).x(2,:),'k') %forward continuation
plot3(run(2).x(3,:),run(2).x(1,:),run(2).x(2,:),'b') %backward continuation

xlabel('tau') %put name of parameter that is being varied (from section 2)
ylabel('x')
zlabel('y')

%% 6. Find LP (limit point = saddle-node bifurcation) coordinates from stored data
 
run_index = 1;
LPpoint=find(strcmp({run(run_index).s.label},'LP  ')) %where LP occurs in run(-).s list
l=length(LPpoint); % l>1 if multiple LP points
LPindex=zeros(1,l); %where along continuation curve NS occur(s)
xx=zeros(2,l);
pp=zeros(2,l);
%save coordinates for input in initializing NS continuation.
for i=[1:l] %(pull out info for each LP point)
    LPindex(i)=run(run_index).s(LPpoint(i)).index;
    %xx will be coordinates for LP
    xx(:,i)=run(run_index).x(1:2,LPindex(i));
    %pp will be parameters for LP
    pp(:,i)=params; %initial parameter values
    pp(contparam,i)=run(run_index).x(3,LPindex(i)); %update
end
LPindex
xx
pp


%% 7. LP (saddle-node) continuation

contparam2 = [1,2]; %continuing in tau and u
i=1; %choose which of NS points (if multiple) to continue
[x0,v0] = init_LPm_LPm(@Stommel,xx(:,i),pp(:,i),contparam2,1);

%indices again track direction of continuation
%forward
jindex=1 
opt = contset;
opt = contset(opt, 'MaxStepSize', 0.5);
opt=contset(opt,'MaxNumPoints',100);
opt=contset(opt,'Singularities',1);
%opt = contset(opt,'IgnoreSingularity',5);
%opt = contset(opt,'Multipliers',1);
opt=contset(opt,'Backward',0);
[lpstruct(jindex).x,lpstruct(jindex).v,lpstruct(jindex).s,lpstruct(jindex).h,lpstruct(jindex).f]=cont(@limitpointmap,x0,v0,opt);

%backward
jindex=2
opt=contset(opt,'MaxNumPoints',50);
opt = contset(opt, 'MaxStepSize', 0.1);
opt=contset(opt,'Backward',1);
[lpstruct(jindex).x,lpstruct(jindex).v,lpstruct(jindex).s,lpstruct(jindex).h,lpstruct(jindex).f]=cont(@limitpointmap,x0,v0,opt);


%% 8. Plot LP curve
figure
hold on
for j = [1,2]
    plot(lpstruct(j).x(3,:),lpstruct(j).x(4,:),'LineWidth',2)
end

legend({'runA','runB'})
xlabel('\tau')
ylabel('k')

%% Afterwards, convert coordinates on LP curve to interface with Fig3 array 
taus1=lpstruct(1).x(3,:);
taus1=20*taus1;
kappas1=lpstruct(1).x(4,:);
kappas1=200-200*kappas1;
taus2=lpstruct(2).x(3,:);
taus2=20*taus2;
kappas2=lpstruct(2).x(4,:);
kappas2=200-200*kappas2;

figure
hold on
plot(taus1,kappas1)
plot(taus2,kappas2)
set(gca, 'YDir','reverse')