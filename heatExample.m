function [E,A,B,N,M,f,n,m,h,s] = heatExample()
  n = 100;
  dx = 1/(n+1);
  nx = [1:n]'*dx;
  h = 1e-3;
  mu = 0.1;
	s = 1000;
  % identity on the left side
    E = speye(n);

    % discretise Laplace
    A = mu*fd_laplace(n,dx);

    % B models BC
    m = 1;
    B = mu*[1; zeros(n-2,1); 1]/dx^2;
    B = sparse(B);

   % bilinear term
    
    N =1*fd_advection(n,dx);
    %Ndiag = -100*sin(10*pi*nx);
    %N = spdiags(Ndiag,0,n,n);

    % Wiener process is of dimension d=2
    d = 2;

    % spatial discretisation of the noise
    
    M = zeros(n,d);
    M(:,1) = exp(-10*(nx-1/2).^2);
    M(:,2) = sin(nx*2*pi);
    M = 0.1*M;
    %M = 20*M/norm(M);
		%M = zeros(size(M));

    % we use uncorrelated noise
    % K = eye(d);

    % drift implicit Euler-Maruyama
    f = @(x0,u,L) (E-h*A-h*N*u)\(x0+h*B*u+sqrt(h)*M*randn(size(M,2),L));
end