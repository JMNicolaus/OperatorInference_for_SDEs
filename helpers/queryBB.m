function X = queryBB(f,x,u,L)
% evolves the black-box dynamical system x_{k+1} = f(x_k,u_k)
  X = zeros(size(x,1),L,size(u,2));
  X(:,:,1) = x;
  for ii=1:size(u,2)-1
    X(:,:,ii+1) = f(X(:,:,ii),u(:,ii),L);
  end
end