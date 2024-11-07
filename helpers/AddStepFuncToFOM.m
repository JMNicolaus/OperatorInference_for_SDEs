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

  E = FOM.E;
  A = FOM.A;
  EAh = E-h*A;
  B = FOM.B;
  M = FOM.M;

  FOM.step = @(x0,u,L) EAh\(x0+h*B*u+sqrt(h)*M*randn(size(M,2),L));

end