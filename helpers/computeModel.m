function [E,C,f] = computeModel(Mdl,x0,Vr,t,u,s,L)
% Mdl is a structure providing either 
%   - a stepping function that can be used with queryBB
% or 
%   - the SDE coefficients E,A,B,N,M
% x0 is the initial condition and of size n-by-1
% Vr is the subspace for lifting/projection
% t are the time steps
% u are the input values at the time steps
% numObs are the number of observations to take
% numSamples is the number of Samples for the MC approximations

MaxLBatchSize = 1e3;
LBatches = L./MaxLBatchSize;
LBatches(LBatches <1)=1;
LBatchSize = min(L,MaxLBatchSize);


switch size(x0,2)
  case 1
    x0 = repmat(x0,1,LBatchSize);
  otherwise
    error("x0 has to be of size N-by-1. The possibility of providing ..." + ...
      "samples of an nondeterministic initial value is not implemented yet.")
end

if isfield(Mdl,'step')
    % this means that Mdl has a stepping function, for example the FOM
    MdlFlag = true;
else
    % Mdl doesnt have a stepping function, for example POD.
    % Provided operators are used and we can use SIEM to get samples
    MdlFlag = false;
    Er = Mdl.E;
    Ar = Mdl.A;
    Br = Mdl.B;
    Bilr = Mdl.Bil;
    Mr = Mdl.M;
end

n = size(x0,1);

E = zeros(n,s);
C = zeros(n,n,s);
f = zeros(3,s);
f1 = zeros(1,s);
f2 = zeros(1,s);
f3 = zeros(1,s);


if MdlFlag
  getStates = @(x0,u) queryBB(Mdl.step,Vr'*x0,u,LBatchSize);
else
  getStates = @(x0,u) permute(semiImplicitEulerStates(Vr'*x0,Er,Ar,Br,Bilr,Mr,u,t,LBatchSize,s),[1 3 2]);
end


parfor kk=1:LBatches
  %disp("Sampling Batch " + kk + " of " + LBatches + " of Model")
  states = getStates(x0,u)
  m = mean(states,2);
  c = 1/(LBatchSize-1)*pagemtimes(states-m,pagetranspose(states-m));
  E = E + 1/LBatches*Vr*permute(m,[1 3 2]);
  C = C + 1/LBatches*pagemtimes(Vr,pagemtimes(c,Vr'));
  states = permute(states,[1 3 2]);
  VrX = pagemtimes(Vr,states);
  f1 = f1 + 1/LBatches*mean(vecnorm(VrX,2).^2,3);
  f2 = f2 + 1/LBatches*mean(max(max(VrX,0),[],1),3);
  f3 = f3 + 1/LBatches*mean(mean(exp(VrX).*(VrX.^3),1),3);
end
f(1,:) = f1;
f(2,:) = f2;
f(3,:) = f3;
end
