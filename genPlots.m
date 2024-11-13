% This script compares the POD and OpInf ROM for the specified examples. 
% The subfigures of Figures 5.2 and 5.3 are produced by this script. 
% To switch betwen the examples set FOM.eqtype to "Heat", "2dHeat" or
% "ConvectionReaction"
% 
% To switch betwen the used snapshot-matrices set snapshotType to "moment"
% or "state"

clear;
rng(0);
addpath(genpath("./"))

%% define FOM example

FOM.eqtype = "Heat";  % set the example
switch 1
  case strcmp(FOM.eqtype,"Heat")
    N = 100;
  case strcmp(FOM.eqtype,"2dHeat")
    N = 20;
  case strcmp(FOM.eqtype,"ConvectionReaction")
    N = 1; % will be overwritten
end

% get matricies of specified example
[FOM.E,FOM.A,FOM.B,FOM.Bil,FOM.M,FOM.K,FOM.ind] = getMatrices(N,1/(N+1),FOM.eqtype);

% time-step size
FOM.h=1e-4;

% add the stepping function of the time-discretised ODE
FOM = AddStepFuncToFOM(FOM);

% make sure we use the correct N. this is only relevant for the 2d Heat ex.
FOM.N = size(FOM.A,1);


%% setup parameters and allocate storage

[n,m] = size(FOM.B);

% number of sample trajectories to generate
L = 1e2;

% number of time steps for each trajectory
s=100;

% the time steps t_0,...,t_{s-1} at which the trajectories are observed
FOM.t = (0:(s-1))*FOM.h;


%% compute subspace

% decide which snapshot matrix to use: "moment" or "state"
snapshotType = "state"; 

% polynomial with random coefficients
u = ppval(spline(linspace(0,s*FOM.h,11),randn(11,1)),FOM.t);

switch 1
  case strcmp(snapshotType,"moment")
    % moment snapshots
    [EV,CV] = computeModel(FOM_reduced,zeros(n,1),eye(n),FOM.t,u,s,L);
    [V,S,~] = svd([EV, reshape(CV,n,[])],"econ");

  case strcmp(snapshotType,"state")
    % state snapshots
    X = queryBB(FOM.step,zeros(n,L),u,L);
    [V,S,~] = svd(reshape(X,n,[]),"econ");
  otherwise
    error("please specify snapshotType as 'moment' or 'state'.")
end

% the ranks for which we want to compute the ROMs. 
rmax = 20;
ranks = [1:rmax];

% the subspace of is spanned by the columns of Vr
Vr = V(:,ranks);

%% generate training data

% create storage objects
FOM.EObs = cell(1,m+rmax);
FOM.CObs = cell(1,m+rmax);
FOM.uObs = cell(1,m+rmax);

% FOM structure without EObs,CObs and uObs fields
% using this instead of FOM makes the computeModel not slower as these
% fields are filled
FOM_reduced = rmfield(FOM, {'EObs', 'CObs', 'uObs'});

% iterate over linearly independent initial condition - input combinations
% to ensure that the data-matrix has full column-rank
X0 = Vr;
for ii=1:(m+rmax)
    disp("ii=" + ii + " of " + (m+rmax))
    u = zeros(m,s);
    if ii<=m
        u(ii,:) = ones(1,s);
        x0 = zeros(n,1);
    else
        x0 = X0(:,ii-m);
    end
    [FOM.EObs{ii},FOM.CObs{ii}] = computeModel(FOM_reduced,x0,eye(n),FOM.t,u,s,L);
    FOM.uObs{ii} = u;
end




%% construct ROMs



% construct the ROMs
[ROMs] = buildROMs(FOM,V(:,ranks));


%% test ROMs

% testing parameters: L= samples size, s = number of time-steps
LTest = 1e4;
sTest = s;
tTest = [0:(sTest-1)]*FOM.h;

% input and initial conditions
uTest = rand*ones(m,sTest);
x0Test = zeros(n,1);

% compute FOM reference
[ExpFOM,CovFOM,fFOM] = computeModel(FOM_reduced,x0Test,eye(FOM.N),tTest,uTest,sTest,LTest);

% compute errors of the ROMs
[errE,errC,errf] = testROMs(ROMs,...
  V,ranks,ExpFOM,CovFOM,fFOM,x0Test,tTest,uTest,sTest,LTest);

% store errors
errors.errE = errE;
errors.errC = errC;
errors.errf = errf;
errors.FOMeqtype = FOM.eqtype;
errors.L = LTest;
errors.s = sTest;

if ~exist("./data","dir")
  mkdir("./data")
end
save("./data/err"+errors.FOMeqtype,'-struct','errors','-v7.3');


%% plot

% singular values of snapshot matrix
f1 = figure(1);
grid on
semilogy(diag(S))
xlabel('$\sigma_r$','Interpreter','latex')
ylabel('$\sigma_i$','Interpreter','latex')
title("singular values of snapshot matrix "+ FOM.eqtype + ' equation','Interpreter','latex')
set(f1,'Position',[100 100 500 500])

% error in expectation
f2 = figure(2);
errEmat = cell2mat(errE);
errEmat(errEmat >=1) = NaN;
semilogy(ranks,errEmat(:,1),'b-s','LineWidth',2)
hold on
plot(ranks,errEmat(:,2),'r-o','LineWidth',2)
hold off
grid on
xlabel('ROM dimension r','Interpreter','latex')
ylabel('relative error','Interpreter','latex')
title("relative errors of expectation, " + FOM.eqtype + ' equation','Interpreter','latex')
legend(["POD", "OpInf"])
axis([1 max(ranks) 1e-5 1e0])
set(f2,'Position',[100 100 500 500])

% error in covariance
f3 = figure(3);
errCmat = cell2mat(errC);
errCmat(errCmat >=1) = NaN;
semilogy(ranks,errCmat(:,1),'b-s','LineWidth',2)
hold on
plot(ranks,errCmat(:,2),'r-o','LineWidth',2)
hold off
grid on
xlabel('ROM dimension r','Interpreter','latex')
ylabel('relative error','Interpreter','latex')
title("relative errors of covariance, " + FOM.eqtype + ' equation','Interpreter','Latex')
legend(["POD", "OpInf"])
axis([1 max(ranks) 1e-5 1e0])
set(f3,'Position',[100 100 500 500])


% extract weak errors of \Phi_1 at last end-time T
for ii=1:numel(ranks)
  errfPlot(ii,:) = errf{ii,1}(:,end)';
end

% weak error for \Phi_1
f4 = figure(4);
semilogy(ranks,errfPlot(:,1),'b-s','LineWidth',2)
hold on
plot(ranks,errfPlot(:,2),'r-o','LineWidth',2)
hold off
grid on
xlabel('ROM dimension r','Interpreter','latex')
ylabel('relative error','Interpreter','latex')
title("relative weak error $e_{\Phi_1}$, " + FOM.eqtype + ' equation','Interpreter','latex')
legend(["POD", "OpInf"])
axis([1 max(ranks) 1e-5 1e0])
set(f4,'Position',[100 100 500 500])


% extract weak errors of \Phi_2 at last end-time T
for ii=1:numel(ranks)
  errfPlot(ii,:) = errf{ii,2}(:,end)';
end

% weak error for \Phi_2
f5 = figure(5);
semilogy(ranks,errfPlot(:,1),'b-s','LineWidth',2)
hold on
plot(ranks,errfPlot(:,2),'r-o','LineWidth',2)
hold off
grid on
xlabel('ROM dimension r','Interpreter','latex')
ylabel('relative error','Interpreter','latex')
title("relative weak error $e_{\Phi_2}$, " + FOM.eqtype + ' equation','Interpreter','latex')
legend(["POD", "OpInf"])
axis([1 max(ranks) 1e-5 1e0])
set(f5,'Position',[100 100 500 500])
