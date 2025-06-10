function [Ehat,Ahat,Bhat,Nhat] = infer_drift(E_train,u_train,h,isbilinear,s)

r = size(E_train{1}(:,1),1);
D = [];
rhs = [];
m = size(u_train{1},1);
p = numel(u_train);
for ii=1:p
  % central finite difference quotient of accuracy 2 along dimension 2
  [Er,Er_dot,ind] = central_finite_differences(E_train{ii},h,2,2);
  u = u_train{ii};
  for jj=1:min(numel(ind),s)
   if isbilinear
      % system has bilinear term
      D_new = [D,[Er(:,jj); u(:,jj); kron(u(:,jj),Er(:,jj))]];
    else
      % system does not have a bilinear term
      D_new = [D,[Er(:,jj); u(:,jj)]];
   end
   if ii<2%r+m
     % use full trajectories
     D = D_new;
   else
     % only use this time point if it reduces the condition number
     if cond(D_new)<cond(D)
       D = D_new;
     else
       continue
     end
   end
  rhs = [rhs,Er_dot(:,jj)];
  end
end
disp("cond(D) = " + cond(D))

% solve least-squares problem

ops = rhs/D;
Ahat = ops(1:r,1:r);
Bhat = ops(:,r+1:r+m);
if isbilinear
  Nhat = ops(:,r+m+1:end);
else 
  Nhat = zeros(size(Ahat));
end
% in our examples E is always the identity
Ehat = speye(r);
end

