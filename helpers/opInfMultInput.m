function [Ahat,Bhat,Nhat,Diffhat,Khat] = opInfMultInput(Ecell,Ccell,ucell,Vr,tObs,isBil,sExp,sCov,dmax)
% Ecell and Ccell are cell array of size 1 x pmax
% each cell in Ecell and Ccell contains observation data of size r x s
% the inferred operators have the sizes:
% A: r x r
% B: r x m
% N: r x r x m
% Diffhat: r x k
% Khat: k x k
% tObs is a vector of times, at which the observations were taken.
% ifBil is a boolean used to enable/disable the inference of a bilinearity
% sExp is the number of observations to use for the expectation
% sCOv is the number of observations to use for the covariance
% dmax is the maximum dimension of the Wiener process

if ~exist("isBil","var")
  isBil = false;
end
if ~exist("sExp","var")
  sExp = numel(tObs)-1; % must be at least r+m+r*m
end
if ~exist("sCov","var")
  sCov = numel(tObs)-5;
end
if ~exist("dmax","var")
  dmax = size(Vr,2);
end


% number of initial condition - input pairs that were observed
pmax = numel(ucell);
m = size(ucell{1},1);

% rank of ROM
r = size(Vr,2);

% time difference between observations
dt = diff(tObs);

D = [];
rhs = [];

% construct least-squares matricies
for ii = 1:pmax
  if isempty(Ecell{ii})
    continue
  end
  % get state expectation
  Er = Ecell{ii}(1:r,:);
  % approximate time derivative
  [x,dxdt,ind] = central_finite_differences(Er(:,1:sExp),dt(1),2);
  % get supplied control
  u = ucell{ii}(:,ind);
  % construct data-matrix
  for jj=1:numel(ind)
    if isBil
      % system has bilinear term
      Dnew = [D,[x(:,jj); u(:,jj); kron(u(:,jj),x(:,jj))]];
    else
      % system does not have a bilinear term
      Dnew = [D,[x(:,jj); u(:,jj)]];
    end
    D = Dnew;
    rhs = [rhs,dxdt(:,jj)];
  end
end

% regularize by computing truncated SVD
disp("rank = "+ r + " condition number:" + cond(D))
tol = 1e-4;
[UD,SD,VD] = svd(D);
rD = find(diag(SD)./SD(1,1)>= tol,1,'last');


% solve for operators
operators = rhs*VD(:,1:rD)*diag(1./diag(SD(1:rD,1:rD)))*UD(:,1:rD)';

% extract operators
Ahat = operators(:,1:r);
Bhat = operators(:,r+1:r+m);
if isBil
  Nhat = operators(:,r+m+1:end);
else
  Nhat = zeros(r,r*m);
end

%% get covariances

Cres=zeros(r*pmax*sCov,r);

% compute residuals
for ii=1:pmax
   if isempty(Ccell{ii})
    continue
  end
  Cr = Ccell{ii}(1:r,1:r,:);
  Crdot = diff(Cr,1,3)./reshape(dt,1,1,[]);
  for jj = 1:sCov
    Clyap = (Ahat+Nhat*kron(u(:,jj),eye(r)))*Cr(:,:,jj) + Cr(:,:,jj)*(Ahat+Nhat*kron(u(:,jj),eye(r)))';
    Cres((ii-1)*sCov*r + (jj-1)*r+1  : (ii-1)*sCov*r + jj*r,:) = Crdot(:,:,jj)-Clyap;
  end
end

% solve least squares problem
H = repmat(eye(r),pmax*sCov,1)\Cres;

% make sure H is symmetric
Hsym = H/2+H'/2;

% get Eigendecomposition of H and sort
[U,S] = eig(Hsym);
[Ssort,I] = sort(diag(S),'descend');
S = diag(Ssort);
U = U(:,I);

% get the index of the last singular value that is greater or equal
% than tol*\sigma_max: (corresponds to 2-norm)
tol = max(S,[],'all')/1000;
d = find(diag(S)>= tol,1,'last');
d = min(d,dmax);
if isempty(d)
  d=0;
end

% truncate H and decompose
Diffhat = U(:,1:d) * sqrt(S(1:d,1:d));
Khat = eye(d);


end
