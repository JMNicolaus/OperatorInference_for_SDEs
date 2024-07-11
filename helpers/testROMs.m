function [relErrE,relErrC,relErrf] = testROMs(ROMs,V,ranks,ExpFOM,CovFOM,fFOM,x0,t,u,numObs,numSamples)

% ROMs is a n-by-m cell. n 'types' of ROMs for 
% the m ranks specified in the array 'ranks'

[m,n] = size(ROMs);

Exp = cell(m,n);
Cov = cell(m,n);
f = cell(m,n);

relErrE = cell(m,n);
relErrC = cell(m,n);
relErrf = cell(m,n);


for ii=1:n
% iterate over ROMs
  for rr=1:m
    %iterate over ranks
    disp("Testing ROM " + ii + " of rank r="+ranks(rr));
    Vr = V(:,1:ranks(rr));
    Mdl = ROMs{rr,ii};
    [Exp{rr,ii},Cov{rr,ii},f{rr,ii}] = computeModel(Mdl,x0,Vr,t,u,numObs,numSamples);  
    relErrE{rr,ii} = norm(ExpFOM-Exp{rr,ii},"fro")/norm(ExpFOM,"fro");
    relErrC{rr,ii} = norm(CovFOM-Cov{rr,ii},"fro")/norm(CovFOM,"fro");
    relErrf{rr,ii} = (abs(fFOM-f{rr,ii}))./abs(fFOM);
    disp("relErrE = "+relErrE{rr,ii} + ", relErrC="+relErrC{rr,ii} + ", relErrf = "+relErrf{rr,ii}(:,end))
  end
end

end