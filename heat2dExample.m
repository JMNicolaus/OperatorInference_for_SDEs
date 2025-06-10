function [E,A,B,N,M,f,n,m,h,s] = heat2dExample()
  n = 34;
  dx = 1/(n+1);
  h = 1e-2;
  mu1 = 1e-2;
	mu2 = 0;% 1e-1;
  s = 100;
  % 1d Laplace
    A1d = mu1*sparse(fd_laplace(n,dx));
    %A1d = A1d + mu2*fd_advection(n,dx);
    % 2d Laplace
    Afull = kron(speye(n),A1d)+kron(A1d,speye(n));

    % define coords of the rectangle that is removed
    cwidth = 10;
    clength = 14;
    cstart = n/2;
    cind = ([cstart:cstart+cwidth-1]) + ([cstart:cstart+clength-1]'*n);

    % solution is only computed on these indices. the complement is 0.
    ind = setdiff(1:n^2,cind);

    % dynamics in rectangle removed
    A = Afull(ind,ind);
    E = speye(size(A));

    % define input rectangle
    istart = 1;
    k = 12;
    Bfull = zeros(n^2,1);
    for ii=n/2-k:n/2+k
      Bfull((istart:(istart+k-1)) + n*(ii-1)) = 1;
    end
    B = sparse(Bfull(ind,:));
    [n,m] = size(B);

    % no bilinearity
    N = sparse(zeros(size(A)));

    % the 1d noise acts on the same rectangle as the input
    %M = 0.01*B;
    M = B/norm(B);
		%M = zeros(size(M));
    K = 1;
  
    % drift implicit Euler-Maruyama
    f = @(x0,u,L) (E-h*A-h*N*u)\(x0+h*B*u+sqrt(h)*M*randn(size(M,2),L));
end
