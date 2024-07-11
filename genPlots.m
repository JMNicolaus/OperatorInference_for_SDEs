clear;
rng(0);
addpath(genpath("./"))
addpath(genpath("~/matlab2tikz/src"))

%% define FOM example

FOM.eqtype = "Heat";  % set the example
snapshotType = "moment"; % decide which snapshot matrix to use: "moment" or "state"
switch 1
    case strcmp(FOM.eqtype,"Heat")
        N = 100;
    case strcmp(FOM.eqtype,"2dHeat")
        N = 20;
    case strcmp(FOM.eqtype,"AdvDiff")
        N = 100;
    otherwise
        N = 1; % will be overwritten
end

[FOM.E,FOM.A,FOM.B,FOM.Bil,FOM.M,FOM.K,FOM.ind] = getMatrices(N,1/(N+1),FOM.eqtype);

% time-step size
FOM.h=1e-4;

% add the stepping function of the time-discretised ODE
FOM = AddStepFuncToFOM(FOM);

% make sure we use the correct N. this is only relevant for the 2d Heat ex.
FOM.N = size(FOM.A,1);


%% generate data

[n,m] = size(FOM.B);

% number of sample trajectories to generate
L = 1e2;

% number of time steps for each trajectory
s=10;

% the time steps t_0,...,t_{s-1} at which the trajectories are observed
FOM.t = (0:(s-1))*FOM.h;

% create storage objects
FOM.EObs = cell(1,m+n);
FOM.CObs = cell(1,m+n);
FOM.uObs = cell(1,m+n);

% iterate over linearly independent initial condition - input combinations
% to ensure that the data-matrix has full column-rank
X0 = eye(n);
X0 = X0(:,randperm(n));
for ii=1:(m+n)
    disp("ii=" + ii + " of " + (m+n))
    u = zeros(m,s);
    if ii<=m
        u(ii,:) = ones(1,s);
        x0 = zeros(n,1);
    else
        x0 = X0(:,ii-m);
    end
    [FOM.EObs{ii},FOM.CObs{ii},~] = computeModel(FOM,x0,eye(n),FOM.t,u,s,L);
    FOM.uObs{ii} = u;
end

% polynomial with random coefficients
u = ppval(spline(linspace(0,s*FOM.h,11),randn(11,1)),FOM.t);

switch 1
  case strcmp(snapshotType,"moment")
    % moment snapshots
    [EV,CV] = computeModel(FOM,zeros(n,1),eye(n),FOM.t,u,s,L);
    [V,S,~] = svd([EV, reshape(CV,n,[])],"econ");

  case strcmp(snapshotType,"state")
    % state snapshots
    X = queryBB(FOM.step,zeros(n,L),u,L);
    [V,S,~] = svd(reshape(X,n,[]),"econ");
  otherwise
    error("please specify snapshotType as 'moment' or 'state'.")
end


%% construct ROMs

ranks = [1:20];
[ROMs] = buildROMs(FOM,V(:,ranks));


%% test ROMs

%LTest = 1e6;
LTest = 1e4;
sTest = s;
tTest = [0:(sTest-1)]*FOM.h;
uTest = rand*ones(m,sTest);
x0Test = zeros(n,1);
[ExpFOM,CovFOM,fFOM] = computeModel(FOM,x0Test,eye(FOM.N),tTest,uTest,sTest,LTest);
[errE,errC,errf] = testROMs(ROMs,...
  V,ranks,ExpFOM,CovFOM,fFOM,x0Test,tTest,uTest,sTest,LTest);

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

f1 = figure(1);
grid on
semilogy(diag(S))
xlabel('$\sigma_r$','Interpreter','latex')
ylabel('$\sigma_i$','Interpreter','latex')
title("singular values of snapshot matrix "+ FOM.eqtype + ' equation','Interpreter','latex')
set(f1,'Position',[100 100 500 500])


f2 = figure(2);
errEmat = cell2mat(errE);
errEmat(errEmat >=1) = NaN;
semilogy(ranks,errEmat(:,1),'k-o','LineWidth',2)
hold on
plot(ranks,errEmat(:,2),'r--x','LineWidth',2)
hold off
grid on
xlabel('ROM dimension r','Interpreter','latex')
ylabel('relative error','Interpreter','latex')
title("relative errors of expectation, " + FOM.eqtype + ' equation','Interpreter','latex')
legend(["POD", "OpInf"])
axis([1 max(ranks) 1e-5 1e0])
set(f2,'Position',[100 100 500 500])


f3 = figure(3);
errCmat = cell2mat(errC);
errCmat(errCmat >=1) = NaN;
semilogy(ranks,errCmat(:,1),'k-o','LineWidth',2)
hold on
plot(ranks,errCmat(:,2),'r--x','LineWidth',2)
hold off
grid on
xlabel('ROM dimension r','Interpreter','latex')
ylabel('relative error','Interpreter','latex')
title("relative errors of covariance, " + FOM.eqtype + ' equation','Interpreter','latex')
legend(["POD", "OpInf"])
axis([1 max(ranks) 1e-5 1e0])
set(f3,'Position',[100 100 500 500])
