function [FOM] = AddStepFuncToFOM(FOM)
% This function creates a "stepping function" FOM.step
% that can be called to advance the time-discretised dynamical system
% one time step.
%
% The FOM structure contains the matricies for computing
% E*x' = (A*x+B*u)*dt + M*dW(t)
% The time-discretisation is performed by a drift-implicit Euler-Maruyama
% method.

% compute time-step size h
  if ~isfield(FOM,'h')
    h = FOM.T/FOM.tsteps;
    FOM.h = h;
  else 
      h = FOM.h;
  end

  A = FOM.A;
  E = eye(size(A))/h-A;
  [Q,R] = qr(E);

  A = Q'*1/h;
  B = Q'*FOM.B;
  M = Q'*FOM.M*1/sqrt(h);
  FOM.step = @(x0,u,L) R\ (A*x0+B*u+M*randn(size(M,2),L));

end