function [Ahat,Bhat,Nhat,Diffhat,Khat] = opInfMultInput(Ecell,Ccell,Vr,ucell,tObs,sExp,sCov,dmax)
% Ecell and Ccell are cell array of size 1 x pmax
% each cell in Ecell and Ccell contains observation data of size r x s
% the inferred operators have the sizes:
  % A: r x r
  % B: r x m 
  % N: r x r x m                
  % Diffhat: r x k
  % Khat: k x k
% tObs is a vector of times, at which the observations were taken.
% sExp is the number of observations to use for the expectation
% sCOv is the number of observations to use for the covariance
% dmax is the maximum dimension of the Wiener process

    if ~exist("sExp","var")
        sExp = numel(tObs)-1; % must be at least r+m+r*m
    end
    if ~exist("sCov","var")
        sCov = numel(tObs)-1;
    end
    if ~exist("dmax","var")
      dmax = size(Vr,2);
    end

% number of initial condition - input pairs that were observed
pmax = numel(ucell);

% rank of ROM
r = size(Vr,2);

% time difference between observations
dt = diff(tObs);

D = [];
rhs = [];
p = [];

% construct least-squares matricies
for ii = 1:pmax 
  Er = Ecell{ii}(1:r,:);
  u = ucell{ii};
  [x,dxdt,ind] = central_finite_differences(Er(:,1:sExp),dt(1),2);

  % in our examples the bilinearity is alway 0
  Dnew = [D,[x;u(ind)]];
  p = [p,ii];
  D = Dnew;
  rhs = [rhs,dxdt];
end
% check condition number
%disp("r="+num2str(size(Vr,2))+", condition: " + num2str(cond(D)))

% solve for operators
operators = D'\rhs';

% extract operators
Ahat = operators(1:r,:)';
Bhat = operators(r+1:r+1,:)';

% in our examples the bilinearity is alway 0
Nhat = zeros(size(Ahat));

%% get covariances

Cres=zeros(r*pmax*sCov,r);
H = zeros(r,r);

% compute residuals
for ii=1:pmax
  Cr = Ccell{p(ii)}(1:r,1:r,:);
  Crdot = diff(Cr,1,3)./reshape(dt,1,1,[]);
  for jj = 1:sCov
    Clyap = (Ahat+Nhat*u(jj))*Cr(:,:,jj) + Cr(:,:,jj)*(Ahat+Nhat*u(jj))';
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

