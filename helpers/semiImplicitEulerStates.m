function [X,u] = semiImplicitEulerStates(x0,E,A,B,N,M,U,t,L,numObs)
% simulates the SDE dX = AX + Bu + sum(N_i*u_i*X,3) dt + M dW
% and returns the states of X(t,w)
% t(i) is time t_i
% U(:,i) is input at time t_i
% L is number of samples
% numObs is the number of observations to make. length(t) >= numObs. 

assert(length(t)>=numObs,"length(t) >= numObs not satisfied");

n = size(A,1);
T = length(t);
k = size(M,2);

% Default number of obs is 1e4
if ~exist("numObs","Var")
  numObs = 1e4; 
end
numSteps = length(U);
q = numSteps/numObs;


X = zeros(n,numObs,L);
u = zeros(size(B,2),numObs);
dt = diff(t);

switch size(x0,2) 
  case 1
    state = repmat(x0,1,L);
  case L
    state = x0;
  otherwise
    error("x0 has to be of size n-by-1 or n-by-L")
end
X(:,1,:) = state;
u(:,1) = U(:,1);
if range(dt)<=1e-15 && range(U,2)==zeros(size(U,1)) % dt and U have constant entries

  [Q,R] = qr(E/dt(1)-A-N*U(:,1));
  Bu = Q'*B*U(:,1);
  M = Q'*1/sqrt(dt(1))*M;
  E = Q'*E;

  for ii=2:T
    b = E*state/dt(1)+Bu+ M*randn(k,L);
    state = R\(b);
    if mod(ii,q) == 0
      jj = ii*numObs/numSteps;
      X(:,jj,:) = state;
      u(:,jj) = U(:,ii);
    end
  end
else
  for ii=2:T
    b = E*state/dt(ii-1)+B*U(:,ii)+ 1/sqrt(dt(ii-1))*M*randn(k,L);
    G = (E/dt(ii-1)-A-N*U(:,ii));
    state = G\b;
    if mod(ii,q) == 0
      jj = ii*numObs/numSteps;
      X(:,jj,:) = state;
      u(:,jj) = U(:,ii);
    end
  end
end
end