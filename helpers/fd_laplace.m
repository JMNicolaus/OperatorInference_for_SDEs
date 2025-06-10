function A = fd_laplace(N,dx)
  % finite-difference discretisation of the 1d Laplace operator
  A = (-2*diag(ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1))/dx^2;
  A = sparse(A);
end