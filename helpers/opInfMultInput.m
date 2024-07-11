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

    if ~exist("sExp","var")
        sExp = numel(tObs)-1; % must be at least r+m+r*m
    end
    if ~exist("sCov","var")
        sCov = numel(tObs)-1;
    end
    if ~exist("dmax","var")
      dmax = size(Vr,2);
    end


pmax = numel(ucell); % number of inputs
r = size(Vr,2);
dt = diff(tObs);

D = [];
rhs = [];
p = [];

for ii = 1:pmax
  
  Er = Vr'*Ecell{ii};
  u = ucell{ii};
  %uKronEr = u.*Er;
  [x,dxdt,ind] = central_finite_differences(Er(:,1:sExp),dt(1),2);
  Dnew = [D,[x;u(ind)]];
%  Dnew = [D,[Er(:,1:sExp); u(1:sExp)]]; %uKronEr(:,1:sExp)]];
   %if cond(Dnew)<cond(D) ||ii<=2 % only select inputs that are beneficial
     p = [p,ii];
     D = Dnew;
 %    Erdot = diff(Er,1,2)./dt;
    rhs = [rhs,dxdt];
%     rhs =[rhs,Erdot(:,1:sExp)];   
   %end
end
disp("r="+num2str(size(Vr,2))+", condition: " + num2str(cond(D)) + ", used inputs: " + num2str(p))
operators = D'\rhs';

Ahat = operators(1:r,:)';
Bhat = operators(r+1:r+1,:)';
Nhat = zeros(size(Ahat));%operators(r+2:end,:)';

%% get covariances

numP = numel(p);
Cres=zeros(r*numP*sCov,r);
H = zeros(r,r);
for ii=1:numP
  Cr = pagemtimes(pagemtimes(Vr',Ccell{p(ii)}),Vr);
  Crdot = diff(Cr,1,3)./reshape(dt,1,1,[]);
%   [y,dydt,indy] = central_finite_differences(Cr(:,:,1:sCov),dt(1),8,3);
%   u = ucell{p(ii)};
%   psi = Ahat+Nhat.*reshape(u(indy),1,1,[]);
%   psiC = psi.*y;
%   H = H + 1/numP*(mean(dydt-psiC-pagetranspose(psiC),3));
  for jj = 1:sCov
    Clyap = (Ahat+Nhat*u(jj))*Cr(:,:,jj) + Cr(:,:,jj)*(Ahat+Nhat*u(jj))';
    Cres((ii-1)*sCov*r + (jj-1)*r+1  : (ii-1)*sCov*r + jj*r,:) = Crdot(:,:,jj)-Clyap;
  end


end
H = repmat(eye(r),numP*sCov,1)\Cres;
Hsym = H/2+H'/2;
[U,S] = eig(Hsym);
[Ssort,I] = sort(diag(S),'descend');
S = diag(Ssort);
U = U(:,I);

% get the index of the last singular value that is greater or equal
% - tol*\sigma_max: (corresponds to 2-norm)
tol = max(S,[],'all')/1000;
d = find(diag(S)>= tol,1,'last'); 
d = min(d,dmax);
if isempty(d) 
  d=0;
end

% - sum of remaining eigenvalues is lower than tol: (Frobenius norm)
% tol = 1e-3;
% d=find(sum(diag(S))-cumsum(diag(S))<=tol,1,'first');

% truncate H and decompose
Diffhat = U(:,1:d) * sqrt(S(1:d,1:d));
Khat = eye(d);





%[Diffhat,d] = chol(sparse(H),'lower');
%Khat = eye(d-1);

% modified cholesky for indefinite matrices
% [Diffhat,K] = mchol(H);
% Diffhat = Diffhat*sqrt(K);
% Khat = eye(size(K));

% Pierre's proposal: just use one of the residuals
% Cr(:,:,1) = Vr'*Ccell{1}(:,:,1)*Vr;
% Cr(:,:,2) = Vr'*Ccell{1}(:,:,2)*Vr;
% Crdot = (Cr(:,:,2)-Cr(:,:,1))/dt;
% Clyap = (Ahat+Nhat*u(ii))*Cr(:,:,1) + Cr(:,:,1)*(Ahat+Nhat*u(1))';
%[Diffhat,d] = chol(sparse(Crdot-Clyap),'lower');


end

