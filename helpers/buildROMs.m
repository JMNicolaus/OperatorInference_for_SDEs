function [ROMs] = buildROMs(FOM,V)
% This function constructs the POD and OpInf ROMs up to r=rank(V)=size(V,2)
% ROM1 = POD ROM
% ROM2 = OpInf ROM
% The object ROMs is an rx2 array of cells 
% The columns ROMs(i,1) and ROMs(i,2) contain the POD and OpInf ROM of rank
% r = i, respectively. 
% The entries of ROMs(i,j) contain only the respective system-operators and
% not a stepping function. This way we are not fixed on one time-stepping
% size for testing.

rmax = size(V,2);
ROMPOD = cell(1,rmax);
ROMOI = cell(1,rmax);

for rr=1:rmax
  Vr = V(:,1:rr);

  % construct POD ROM
  ROMPOD{rr}.E = Vr'*FOM.E*Vr;
  ROMPOD{rr}.A = Vr'*FOM.A*Vr;
  ROMPOD{rr}.B = Vr'*FOM.B;
  ROMPOD{rr}.Bil = Vr'*FOM.Bil*Vr;
  ROMPOD{rr}.M = Vr'*FOM.M;

  % construct OpInf ROM
  ROMOI{rr}.E = eye(rr);
  [ROMOI{rr}.A,ROMOI{rr}.B,ROMOI{rr}.Bil,ROMOI{rr}.M,~] = opInfMultInput(FOM.EObs,FOM.CObs,Vr,FOM.uObs,FOM.t);

end

ROMs = cell(rmax,2);
ROMs(:,1) = ROMPOD;
ROMs(:,2) = ROMOI;
