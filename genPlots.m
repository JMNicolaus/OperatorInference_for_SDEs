% This script compares the POD and OpInf ROM for the specified examples.
% The subfigures of Figures 5.2 and 5.3 are produced by this script.
% To switch betwen the examples set FOM.eqtype to "Heat", "2dHeat" or
% "ConvectionReaction"
%
% To switch betwen the used snapshot-matrices set snapshotType to "moment"
% or "state"

clear all;
close all;
rng(0);
addpath(genpath("./"))

%% setup parameters

% select experiment from "Heat", "2dHeat", "ConvectionReaction"
FOM.eqtype = "Heat";

% selct which snapshot matrix to use: "moment" or "state"
snapshotType = "state";

% set maximum rank
r_max = 20;

% number of samples
L_subspace  = 1e2;
L_train     = 1e2;
L_test      = 1e4;

% number of time-steps s, time-step size h and the times themselves t
s = 100;
h = 1e-4;
t = h*[0:s-1];

%% define FOM example

switch 1
  case strcmp(FOM.eqtype,"Heat")
    Ntemp = 1000;
  case strcmp(FOM.eqtype,"2dHeat")
    Ntemp = 20; % number of grid point for one spatial dimension
  case strcmp(FOM.eqtype,"ConvectionReaction")
    Ntemp = 1; % will be overwritten
end

% get matricies of specified example and set isBil flag
[FOM.E,FOM.A,FOM.B,FOM.Bil,FOM.M,FOM.K,FOM.ind] = getMatrices(Ntemp,1/(Ntemp+1),FOM.eqtype);
FOM.isBil = nnz(FOM.Bil)~=0;

% store time-step size in FOM
FOM.h = h;

% add the stepping function of the time-discretised ODE
FOM = AddStepFuncToFOM(FOM);

% the time steps t_0,...,t_{s-1} at which the trajectories are observed
FOM.t = t;

[N,m] = size(FOM.B);

% set number of extra training pairs depending on if bilinearity is present
if FOM.isBil
    k_extra = 10;
else
    k_extra = 0;
end



%% compute subspace

% polynomial with random coefficients
u = ppval(spline(linspace(0,s*FOM.h,11),randn(11,1)),FOM.t);

switch 1
  case strcmp(snapshotType,"moment")
    % moment snapshots
    [EV,CV] = computeModel(FOM,zeros(N,1),eye(N),FOM.t,u,s,L_subspace);
    [V,S,~] = svd([EV, reshape(CV,N,[])],"econ");

  case strcmp(snapshotType,"state")
    % state snapshots
    X = queryBB(FOM.step,zeros(N,L_subspace),u,L_subspace);
    [V,S,~] = svd(reshape(X,N,[]),"econ");
  otherwise
    error("please specify snapshotType as 'moment' or 'state'.")
end

% the ranks for which we want to compute the ROMs.
ranks = [1:r_max];

% the subspace of is spanned by the columns of Vr
Vr = V(:,ranks);

%% generate training data

% create storage objects
FOM.EObs = cell(1,m+r_max);
FOM.CObs = cell(1,m+r_max);
FOM.uObs = cell(1,m+r_max);

% FOM structure without EObs,CObs and uObs fields
% using this instead of FOM makes the computeModel not slower as these
% fields are filled
FOM_reduced = rmfield(FOM, {'EObs', 'CObs', 'uObs'});

% iterate over linearly independent initial condition - input combinations
% to ensure that the data-matrix has full column-rank

% define controls to train on
%uTrain = cell(1, m + 1 + k_extra);
for ii=1:(m+1)
  utemp = zeros(m,s);
  if ii~=1
    utemp(ii-1,:) = rand*ones(1,s);
  end
  uTrain{ii} = utemp;
end
clear utemp

% define initial conditions to train on
%x0Train = cell(1,r_max+1 +k_extra);
for ii=1:(r_max+1)
  if ii==1
    x0temp = zeros(N,1);
  else
    x0temp = Vr(:,ii-1);
  end
  x0Train{ii} = x0temp;
end
clear x0temp

% train on linearly independent pairs
idx = 1;
for ii=1:m+1
  u = uTrain{ii};
  for jj=1:r_max+1
    x0 = x0Train{jj};
    if (~FOM.isBil && ((ii-1)*(jj-1)~=0)) || (ii==1 && jj==1)
      % We dont need to sample pairs of non-zero IC and non-zero control
      % if the the FOM is not bilinear.
      continue
    end
    disp((ii-1) + " " + (jj-1))
    [EObs_temp,CObs_temp] = computeModel(FOM_reduced,x0,eye(N),t,u,s,L_train);
    % store only the projected moments
    FOM.EObs{idx} = Vr'*EObs_temp;
    FOM.CObs{idx} = pagemtimes(Vr',pagemtimes(CObs_temp,Vr));
    FOM.uObs{idx} = u;
    clear EObs_temp CObs_temp
    idx = idx +1;
  end
end

% generate additional training sets
for kk=1:k_extra

  % additional control and IC to train on
  u = repmat(rand(m,1),1,s);
  x0= Vr(:,randi(r_max));%abs(randn(N,1));

  % save to training set
  uTrain{m+1+kk} = u;
  x0Train{r_max+1+kk} = x0;

  disp(kk)
  [EObs_temp,CObs_temp] = computeModel(FOM_reduced,x0,eye(N),t,u,s,L_train);
  % store only the projected moments
  FOM.EObs{idx+kk-1} = Vr'*EObs_temp;
  FOM.CObs{idx+kk-1} = pagemtimes(Vr',pagemtimes(CObs_temp,Vr));
  FOM.uObs{idx+kk-1} = u;
  clear EObs_temp CObs_temp
end


%% construct ROMs

[ROMs] = buildROMs(FOM,Vr);


%% test ROMs

% input and initial conditions
u_test = rand*ones(m,s);
x0_test = zeros(N,1);

% compute FOM reference
[ExpFOM,CovFOM,fFOM] = computeModel(FOM_reduced,x0_test,eye(N),t,u_test,s,L_test);

% compute errors of the ROMs
[errE,errC,errf] = testROMs(ROMs,...
  V,ranks,ExpFOM,CovFOM,fFOM,x0_test,t,u_test,s,L_test);

%% store errors
errors.eqtype = FOM.eqtype;
errors.snapshotType = snapshotType;
errors.errE = errE;
errors.errC = errC;
errors.errf = errf;
errors.L_subspace = L_subspace;
errors.L_train = L_train;
errors.L_test = L_test;
errors.s = s;
errors.h = h;
errors.r_max = r_max;
errors.k_extra = k_extra;
errors.x0_test = x0_test;
errors.u_test = u_test;

% Create base data directory if necessary 
if ~exist("./data","dir")
    mkdir("./data")
end

% Generate a subfolder name
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
foldername = fullfile('./data', ...
    [num2str(errors.eqtype) '_' timestamp]);

% Create subfolder
mkdir(foldername);

% Save the errors struct 
save(fullfile(foldername, 'errors.mat'), '-struct', 'errors', '-v7.3');


%% plot

% singular values of snapshot matrix
f1 = figure(1);
grid on
semilogy(diag(S))
xlabel('$\sigma_r$','Interpreter','latex')
ylabel('$\sigma_i$','Interpreter','latex')
title("singular values of snapshot matrix "+ errors.eqtype + ' equation','Interpreter','latex')
set(f1,'Position',[100 100 500 500])

% error in expectation
f2 = figure(2);
errEmat = cell2mat(errors.errE);
errEmat(errEmat >=1) = NaN;
semilogy(ranks,errEmat(:,1),'b-s','LineWidth',2)
hold on
plot(ranks,errEmat(:,2),'r-o','LineWidth',2)
hold off
grid on
xlabel('ROM dimension r','Interpreter','latex')
ylabel('relative error','Interpreter','latex')
title("relative errors of expectation, " + errors.eqtype + ' equation','Interpreter','latex')
legend(["POD", "OpInf"])
axis([1 max(ranks) 1e-5 1e0])
set(f2,'Position',[100 100 500 500])

% error in covariance
f3 = figure(3);
errCmat = cell2mat(errors.errC);
errCmat(errCmat >=1) = NaN;
semilogy(ranks,errCmat(:,1),'b-s','LineWidth',2)
hold on
plot(ranks,errCmat(:,2),'r-o','LineWidth',2)
hold off
grid on
xlabel('ROM dimension r','Interpreter','latex')
ylabel('relative error','Interpreter','latex')
title("relative errors of covariance, " + errors.eqtype + ' equation','Interpreter','Latex')
legend(["POD", "OpInf"])
axis([1 max(ranks) 1e-5 1e0])
set(f3,'Position',[100 100 500 500])


% extract weak errors of \Phi_1 at last end-time T
for ii=1:numel(ranks)
  errfPlot(ii,1) = errors.errf{ii,1}(1,end)';
  errfPlot(ii,2) = errors.errf{ii,2}(1,end)';
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
title("relative weak error $e_{\Phi_1}$, " + errors.eqtype + ' equation','Interpreter','latex')
legend(["POD", "OpInf"])
axis([1 max(ranks) 1e-5 1e0])
set(f4,'Position',[100 100 500 500])


% extract weak errors of \Phi_2 at last end-time T
for ii=1:numel(ranks)
  errfPlot(ii,1) = errors.errf{ii,1}(2,end)';
  errfPlot(ii,2) = errors.errf{ii,2}(2,end)';
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
title("relative weak error $e_{\Phi_2}$, " + errors.eqtype + ' equation','Interpreter','latex')
legend(["POD", "OpInf"])
axis([1 max(ranks) 1e-5 1e0])
set(f5,'Position',[100 100 500 500])
