clear;close all;
rng(0);
addpath(genpath("./"))

%% define FOM
FOM.eqtype = "Heat";
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
FOM.h=1e-3;
FOM = AddStepFuncToFOM(FOM);
FOM.N = size(FOM.A,1);
[n,m] = size(FOM.B);

% number of samples
L = 1e4;

% number of time steps
s=1e3;
FOM.t = (0:(s-1))*FOM.h;

% input
u = ones(1,s);
x0 = zeros(n,L);

% obtain samples
X = queryBB(FOM.step,x0,u,L);

% svd of state-snapshots
[U1,S1,~] = svd(reshape(X,n,[]),"econ");

% svd of moment-snapshots
E = squeeze(mean(X,2));
C = page_cov(X,true);
C = reshape(C,n,[]);
[U2,S2,~] = svd([E C],"econ");

% turn the vectors in the same direction
for ii=1:n
  U1(:,ii) = sign(U1(1,ii))*U1(:,ii);
  U2(:,ii) = sign(U2(1,ii))*U2(:,ii);
end

figure(1)
for ii=1:12
  subplot(4,3,ii)
  hold on
  plot(U1(:,ii))
  plot(U2(:,ii))
  hold off
  legend
end

figure(2)
for ii=1:12
  subplot(4,3,ii)
  semilogy(abs(U1(:,ii)-U2(:,ii)))
  legend
end

err = zeros(10);
for ii=1:10
  for jj=1:10
    err(ii,jj) = norm(U1(:,ii)-U2(:,jj));
  end
end
round(log10(err),2)
