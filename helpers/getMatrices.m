function [E,A,B,Bil,Diff,K,ind] = getMatrices(N,dx,type)
% This function constructs the examples specified in the publication
% The FOM takes the form
% E*x' = [A*x+B*u+Bil(:,:,1)*x*u_1 + ... + Bil(:,:,m)*x*u_m]*dt + M*dW(t)
% %
% % E is NxN and is in these examples always the identity.
% % A is NxN and is the discretisation of the linear differential operator
% % B is Nxm and encodes the boundary conditions
% % Bil is NxNxm and is the collection of bilinear termns and is in these
% % examples always 0.
% % d is the dimension of the Wiener process
% % Diff is Nxd and corresponds to the spatial discretisation of the noise
% % K is dxd and specifies the correlation of the noise components
% % ind is a list of grid-indicies on which the solution is computed. this
% only not the full grid in the 2d Heat example.

switch 1
  case strcmp(type,'Heat')
    % identity on the left side
    E = speye(N);

    % discretise Laplace
    A = fd_laplace(N,dx);

    % B models BC
    B = [1/dx^2; zeros(N-2,1); 1/dx^2];

    % no bilinear term
    Bil = zeros(N);

    % Wiener process is of dimension d=2
    d = 2;

    % spatial discretisation of the noise
    Diff = zeros(N,d);
    Nx = [1:N]*dx;
    Diff(:,1) = exp(-10*(Nx-1/2).^2);
    Diff(:,2) = sin(Nx*2*pi);
    Diff = 0.1*Diff;

    % we use uncorrelated noise
    K = eye(d);

    % the solution is computed on all grid-indices
    ind = 1:N;


  case strcmp(type,'2dHeat')
    
    % 1d Laplace
    A1d = sparse(fd_laplace(N,dx));
    % 2d Laplace
    Afull = kron(speye(N),A1d)+kron(A1d,speye(N));

    % define coords of the rectangle that is removed
    cwidth = 5;
    clength = 7;
    cstart = N/2;
    cind = ([cstart:cstart+cwidth-1]) + ([cstart:cstart+clength-1]'*N);

    % solution is only computed on these indices. the complement is 0.
    ind = setdiff(1:N^2,cind);

    % dynamics in rectangle removed
    A = Afull(ind,ind);
    E = speye(size(A));

    % define input rectangle
    istart = 1;
    k = 6;
    Bfull = zeros(N^2,1);
    for ii=N/2-k:N/2+k
      Bfull((istart:(istart+k-1)) + N*(ii-1)) = 1;
    end
    B = sparse(Bfull(ind,:));

    % no bilinearity
    Bil = zeros(size(A));

    % the 1d noise acts on the same rectangle as the input
    Diff = 0.1*B;
    K = 1;


  case strcmp(type,"ConvectionReaction")
    if ~exist("pde.mat","file")
      error("please download the pde.mat file from here: 'https://morwiki.mpi-magdeburg.mpg.de/morwiki/index.php/Convection_Reaction'" + ...
        " and place it in the 'helpers' directory. ")
    end
    % get predefined operators
    pde = load("pde.mat");
    A = pde.A;
    B = pde.B;
    E = speye(size(A));

    % no bilinearity
    Bil = zeros(size(A));

    % the same noise as in the 1d Heat eq.
    d = 2;
    Diff = zeros(size(A,1),d);
    Nx = linspace(0,1,size(A,1))*dx;
    Diff(:,1) = exp(-10*(Nx-1/2).^2);
    Diff(:,2) = sin(Nx*2*pi);
    Diff = 0.1*Diff;
    K = eye(d);

    % the solution is computed on all indices
    ind = 1:size(A,1);
  otherwise
    error("No such equation type found! Heat,AdvDiff and DampedWave are available.")
end
end

function A = fd_laplace(N,dx)
% finite-difference discretisation of the 1d Laplace operator
A = (-2*diag(ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1))/dx^2;
end
