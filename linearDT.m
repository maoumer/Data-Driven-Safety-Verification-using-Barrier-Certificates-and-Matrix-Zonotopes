% linearDT - computes the data driven reachable set of discrete time systems
% x(k+1) = Ax(k) + Bu(k) + w(k)
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: Amr Alanwar, Anne Koch, Frank Allgöwer, Karl Johansson "Data Driven Reachability Analysis Using Matrix Zonotopes"
%
%
%
% Author:       Amr Alanwar
% Written:      28-October-2020
% Last update:  
% Last revision:---

%------------- BEGIN CODE --------------

rand('seed',1);

clear all; close all; clc
%% system dynamics
dim_x = 5;
A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B_ss = ones(dim_x,1);
C = [1,0,0,0,0];
D = 0;
%%%%%%%%%%% To take smaller dimensions if needed, added by Mohammed
dim_x = 5; % change to 2 for running harder system verification method. change line 157/8 too.
A = A(1:dim_x,1:dim_x); B_ss = B_ss(1:dim_x); C = C(1:dim_x);
%%%%%%%%%%%
% define continuous time system
sys_c = ss(A,B_ss,C,D);
% convert to discrete system
samplingtime = 0.05;
sys_d = c2d(sys_c,samplingtime);

%Number of trajectories
initpoints =1;
%Number of time steps
steps = 120;
totalsamples = initpoints*steps;
%% initial set and input
X0 = zonotope(ones(dim_x,1),0.1*diag(ones(dim_x,1)));
U = zonotope(10,0.25);

%noise zontope W
W = zonotope(zeros(dim_x,1),0.005*ones(dim_x,1));

%Construct matrix zonotpe \mathcal{M}_w
index=1;
for i=1:size(W.generators,2)
    vec=W.Z(:,i+1);
    GW{index}= [vec,zeros(dim_x,totalsamples-1)];
    for j=1:totalsamples-1
        GW{j+index}= [GW{index+j-1}(:,2:end) GW{index+j-1}(:,1)];
    end
    index = j+index+1;
end
Wmatzono= matZonotope(zeros(dim_x,totalsamples),GW);


% randomly choose constant inputs for each step / sampling time
for i=1:totalsamples
    u(i) = randPoint(U);
end



%simulate the system to get the data
x0 = X0.center;
x(:,1) = x0;
index=1;
for j=1:dim_x:initpoints*dim_x
    x(j:j+dim_x-1,1) = randPoint(X0);
    for i=1:steps
        utraj(j,i) = u(index);
        x(j:j+dim_x-1,i+1) = sys_d.A*x(j:j+dim_x-1,i) + sys_d.B*u(index) + randPoint(W);      
        index=index+1;
    end
end


% concatenate the data trajectories 
index_0 =1;
index_1 =1;
for j=1:dim_x:initpoints*dim_x
    for i=2:steps+1
        x_meas_vec_1(:,index_1) = x(j:j+dim_x-1,i);
        index_1 = index_1 +1;
    end
    for i=1:steps
        u_mean_vec_0(:,index_0) = utraj(j,i);
        x_meas_vec_0(:,index_0) = x(j:j+dim_x-1,i);
        index_0 = index_0 +1;
    end
end

% X_+ is X_1T
% X_- is X_0T
U_full = u_mean_vec_0(:,1:totalsamples); %same as u 
X_0T = x_meas_vec_0(:,1:totalsamples);
X_1T = x_meas_vec_1(:,1:totalsamples);


% plot simulated trajectory
figure('Renderer', 'painters', 'Position', [10 50 1400 750]);
t = 0:totalsamples;
for i = 1:dim_x
    subplot(1,dim_x,i); hold on; box on; grid on;
    plot(t,x(i,:),'b'); 
    xlabel('t'); 
    ylabel('x_'+string(i));
end
% saveas(gcf,"states_plot"+string(dim_x)+".png") % uncomment if needed
close;


X1W_cen =  X_1T - Wmatzono.center;
X1W = matZonotope(X1W_cen,Wmatzono.generator);

% set of A and B
AB = X1W * pinv([X_0T;U_full]);

% validate that A and B are within AB
intAB11 = intervalMatrix(AB);
intAB1 = intAB11.int;
valid = all(intAB1.sup >= [sys_d.A,sys_d.B],'all') & all(intAB1.inf <= [sys_d.A,sys_d.B],'all')


%% compute next step sets from model / data

% set number of steps in analysis
totalsteps = 5;
X_model = cell(totalsteps+1,1);
X_data = cell(totalsteps+1,1);
% init sets for loop
X_model{1} = X0; X_data{1} = X0;

for i=1:totalsteps
    
    % 1) model-based computation
    X_model{i,1}=reduce(X_model{i,1},'girard',400);
    X_model{i+1,1} = sys_d.A * X_model{i} + sys_d.B * U+W;
    % 2) Data Driven approach
    X_data{i,1}=reduce(X_data{i,1},'girard',400);
    X_data{i+1,1} = AB * (cartProd(X_data{i},U)) +W;
    
    
end


%% visualization

projectedDims = {[1 2],[3 4],[4 5]}; % for 5d system
% projectedDims = {[1 2],[2 1]}; % for 2d system
axx{1} = [0.75,1.5,0.5,4]; axx{2} = [0.75,3,0.8,2.2];axx{3} = [0.75,2.3,0.75,2.8];
index=1;
numberofplots = 5;%length(X_model)
for plotRun=1:length(projectedDims)
    
    figure('Renderer', 'painters', 'Position', [10 50 600 750])
    
    % set axis
      
    index=index+1;
    % plot initial set
    handleX0 = plot(X0,projectedDims{plotRun},'k-','LineWidth',2);
    hold on;
    
    
    % plot reachable sets starting from index 2, since index 1 = X0
    
    % plot reachable sets from model
    for iSet=2:numberofplots
        handleModel=  plot(X_model{iSet},projectedDims{plotRun},'b','Filled',true,'FaceColor',[.8 .8 .8],'EdgeColor','b');
    end
    
    % plot reachable sets from data
    for iSet=2:numberofplots
        handleData= plot(X_data{iSet},projectedDims{plotRun},'r');
    end
    
    % label plot
    xlabel(['x_{',num2str(projectedDims{plotRun}(1)),'}']);
    ylabel(['x_{',num2str(projectedDims{plotRun}(2)),'}']);
    %axis(axx{plotRun});
    % skip warning for extra legend entries
    warOrig = warning; warning('off','all');
    legend([handleX0,handleModel,handleData],...
        'Initial Set','Set from Model','Set from Data','Location','northwest');
    warning(warOrig);
    ax = gca;
    ax.FontSize = 22;
    %set(gcf, 'Position',  [50, 50, 800, 400])
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

end

%% Save relevant variables for Julia - added by Mohammed
close all;
xmin = [-4; 0; 0; 0; 0]; xmax = [2;5;5;3;6]; % state set
xumin = [1.5;0;0;2.5;0]; xumax = xumin + 0.5; % unsafe states
%%%%%%%%%%% To take subset of system - added by Mohammed
xmin = xmin(1:dim_x); xmax = xmax(1:dim_x);
Xint = interval(xmin, xmax);
xumin = xumin(1:dim_x); xumax = xumax(1:dim_x);
%%%%%%%%%%%
X0int = interval(X0); % zonotope to interval of initial set
x0min = X0int.inf; x0max = X0int.sup; % get min and max of interval
Uint = interval(U); % zonotope to interval of input set
umin = Uint.inf; umax = Uint.sup; % get min and max of interval
wcenter = W.center;
Wint = interval(W);
wmin = Wint.inf; wmax = Wint.sup;
ABcenter = AB.center;
Acenter = ABcenter(:,1:end-1); Bcenter = ABcenter(:,end);
ABmin = intAB1.inf; ABmax = intAB1.sup; % see line 129
Amin = ABmin(:,1:end-1); Bmin = ABmin(:,end);
Amax = ABmax(:,1:end-1); Bmax = ABmax(:,end);
Aint = intAB1(:,1:dim_x); Bint = intAB1(:,dim_x+1); % last column B
Ax = Aint*Xint; Bu = Bint*Uint;
AxBuW = Ax + Bu + Wint;
AxBuWmin = AxBuW.inf; AxBuWmax = AxBuW.sup; 
zero = zeros(size(AB.center)); % create zero matrix same size as [A B]
AgBg = matZonotope(zero,AB.generator); % zero out the center of [A B]
intAgBg = intervalMatrix(AgBg); % convert to interval matrix from mat zonotope
intAgBg = intAgBg.int; % change to matrix of intervals
Ag = intAgBg(:,1:dim_x); Bg = intAgBg(:,dim_x+1); % extract generators for A and B from [A B] generators
AgX = Ag*Xint; BgU = Bg*Uint; % dynamics based on generators
uncertainty = AgX + BgU + Wint; %total uncertainty of generators + noise zonotope
uncertaintyMin = uncertainty.inf; uncertaintyMax = uncertainty.sup;

% save variables to file
save("vars_"+string(dim_x)+".mat","Amin","Bmin","Amax","Bmax","x0min","x0max","xumin","xumax","xmin","xmax","umin","umax","wcenter","wmin","wmax","Acenter","Bcenter","uncertaintyMin","uncertaintyMax")


%------------- END OF CODE --------------
