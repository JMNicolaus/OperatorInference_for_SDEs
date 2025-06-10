function [E_train,C_train,u_train] = train3(f,Vr,m,s,L,h,isbilinear)

L_max = max(L);

[n,rmax] = size(Vr);

% generate training data
E_train = cell(numel(L),1);
C_train = cell(numel(L),1);
u_train = cell(numel(L),1);

k = 21;
us = linspace(-2,2,k);

idx = 0;
for ii=1:k
  idx = idx+1;
  disp("x0 = 0, u = e_"+ii);
  x0 = zeros(n,1);
  u = us(ii)*ones(m,s);
  X_train = pagemtimes(Vr',stepSDE(f,x0,u,L_max));
  for jj = 1:numel(L)
    LL = L(jj);
    E_train{jj}{idx} = squeeze(mean(X_train(:,1:LL,:),2));
    C_train{jj}{idx} = page_cov(X_train(:,1:LL,:),true) ;
    u_train{jj}{idx} = u;
  end
end

for ii=1:rmax
  idx = idx +1;
  disp("x0 = Vr(:,"+ii+"), u = 0");
  x0 = Vr(:,ii);
  u = zeros(m,s);
  X_train = pagemtimes(Vr',stepSDE(f,x0,u,L_max));
  for jj = 1:numel(L)
    LL = L(jj);
    E_train{jj}{idx} = squeeze(mean(X_train(:,1:LL,:),2));
    C_train{jj}{idx} = page_cov(X_train(:,1:LL,:),true) ;
    u_train{jj}{idx} = u;
  end
end


end


