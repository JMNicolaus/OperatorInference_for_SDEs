function [E,C,f1,f2] = compute_model(f,V,x0,u,L)

n = size(x0,1);
[m,s] = size(u);

% test parameters
batchSize = 1e4;  % Number of samples per batch
numBatches = ceil(L / batchSize);  % Total number of batches
if L <= batchSize
  numBatches = 1;
  batchSize = L;
end

E = zeros(n,s);
C = zeros(n,n,s);
f1 = 0;
f2 = 0;
for batch = 1:numBatches
  %batch
  %tic
  [E_temp,C_temp,f1_temp,f2_temp] = estimate(f,V,V'*x0,u,batchSize);
  E = E + 1/numBatches*E_temp;
  C = C + 1/numBatches*C_temp;
  f1 = f1 + 1/numBatches*f1_temp;
  f2 = f2 + 1/numBatches*f2_temp;
%toc;
end

end