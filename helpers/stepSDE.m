function X = stepSDE(f,x0,u,L)
  n = size(x0,1);
  s = size(u,2);
  X = zeros(n,L,s);
  X(:,:,1) = f(x0,u(:,1),L);
  for ii=2:size(u,2)
    X(:,:,ii) = f(X(:,:,ii-1),u(ii),L);
  end
end