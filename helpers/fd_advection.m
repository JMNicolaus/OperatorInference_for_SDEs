function A = fd_advection(N, dx)
  A=(diag(ones(N-1,1),1) + diag(-1*ones(N-1,1),-1))/2/dx; %central FD 
  A(1,1:2) =[-1 1]/dx; % Forward FD
  A(end,end-1:end) = [-1 1]/dx; % Backward FD
  A = sparse(A);
end