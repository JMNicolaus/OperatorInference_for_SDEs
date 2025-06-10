function [Mhat,Khat] = infer_diffusion(C_train,u_train,h,Ahat,Nhat)

m = size(u_train{1},1);
p = numel(u_train);

% the solution of the least squares problem is the mean over the residuals 
Hhat = zeros(size(Ahat));
for ii=1:p
  [Cr,Cr_dot,ind] = central_finite_differences(C_train{ii},h,2,3);
  u = u_train{ii};
  for jj=1:numel(ind)
    Clyap = (Ahat+Nhat*kron(u(:,jj),eye(m)))*Cr(:,:,jj);
    Hhat = Hhat + 1/(p*numel(ind))*(Cr_dot(:,:,jj) - Clyap - Clyap'); 
  end
end
% make sure that the H is symmetric
Hhat = Hhat/2+Hhat'/2;

% get Eigendecomposition of H and sort
[HU,HS] = eig(Hhat);
[Ssort,I] = sort(diag(HS),'descend');
HS = diag(Ssort);
HU = HU(:,I);

% get the index of the last singular value that is greater or equal
% than tol*\sigma_max: (corresponds to 2-norm)
tol = max(HS,[],'all')/1000;
d = find(diag(HS)>= tol,1,'last');
if isempty(d)
  d=0;
end

% truncate H and decompose
Mhat = HU(:,1:d) * sqrt(HS(1:d,1:d));
Khat = eye(d);
end

