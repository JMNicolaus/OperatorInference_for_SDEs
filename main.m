clear all;
addpath(genpath("./"))
rng(2)

%[E,A,B,N,M,f,n,m,h,s] = heatExample();
[E,A,B,N,M,f,n,m,h,s] = heat2dExample();
isbilinear = full(any(N,'all'));

% maximal ROM dimension
rmax = 10;


% number of samples
L_subspace = 1e4;
L_train = [1e1 1e2 1e4];
L_test = 1e6;


%% get subspace
disp("Computing subspace Vr")
u_subspace = cos(linspace(0,2*pi,s));
x0_subspace = zeros(n,1);
if ~exist("V","var")
  X_subspace = stepSDE(f,x0_subspace,u_subspace,L_subspace);
  [V,S,~] = svd(reshape(X_subspace,n,[]),"econ");
end
Vrmax = V(:,1:rmax);


%% training
[E_train,C_train,u_train] = train(f,Vrmax,m,s,L_train,h,isbilinear);


%% test set
u_test = cos(linspace(0,5*pi,s));
x0_test = x0_subspace;
seed_test = 42;

% FOM reference
disp("computing FOM reference")
rng(seed_test)
[EFOM,CFOM,f1FOM,f2FOM] = compute_model(f,eye(n),x0_test,u_test,L_test);

%% OpInf
E_error_OpInf = zeros(rmax,numel(u_train));
C_error_OpInf = zeros(rmax,numel(u_train));
f1_error_OpInf = zeros(rmax,numel(u_train));
f2_error_OpInf = zeros(rmax,numel(u_train));

for kk = 1:numel(E_train)
    % infer using different number of samples
    [Ehat,Ahat,Bhat,Nhat] = infer_drift(E_train{kk},u_train{kk},h,isbilinear,s);
    [Mhat,Khat] = infer_diffusion(C_train{kk},u_train{kk},h,Ahat,Nhat);

  for ii=1:rmax
    Vr = Vrmax(:,1:ii);
    Ehatr = Ehat(1:ii,1:ii);
    Ahatr = Ahat(1:ii,1:ii);
    Bhatr = Bhat(1:ii,:);
    Nhatr = Nhat(1:ii,1:ii);
    Mhatr = Mhat(1:ii,:);
    fhatr = @(x0,u,L) (Ehatr-h*Ahatr-h*Nhatr*u)\(x0+h*Bhatr*u+sqrt(h)*Mhatr*randn(size(Mhatr,2),L));
      disp(kk+" " + ii)
      rng(seed_test)
      [EROM,CROM,f1ROM,f2ROM]=compute_model(fhatr,Vr,x0_test,u_test,L_test);
      E_error_OpInf(ii,kk)=norm(EROM-EFOM,"fro")/norm(EFOM,"fro");
      C_error_OpInf(ii,kk)=page_norm(CROM-CFOM)/page_norm(CFOM);
      f1_error_OpInf(ii,kk)=abs(f1ROM-f1FOM)/abs(f1FOM);
      f2_error_OpInf(ii,kk)=abs(f2ROM-f2FOM)/abs(f2FOM);
  end
end


%% POD
E_error_POD = zeros(rmax,1);
C_error_POD = zeros(rmax,1);
f1_error_POD = zeros(rmax,1);
f2_error_POD = zeros(rmax,1);

for ii=1:rmax
  Vr = Vrmax(:,1:ii);
  Er = Vr'*E*Vr;
  Ar = Vr'*A*Vr;
  Br = Vr'*B;
  Nr = Vr'*N*Vr;
  Mr = Vr'*M;
  fr = @(x0,u,L) (Er-h*Ar-h*Nr*u)\(x0+h*Br*u+sqrt(h)*Mr*randn(size(Mr,2),L));

    disp("ROM dimension " + ii )
    rng(seed_test)
    [EROM,CROM,f1ROM,f2ROM]=compute_model(fr,Vr,x0_test,u_test,L_test);
    E_error_POD(ii)=norm(EROM-EFOM,"fro")/norm(EFOM,"fro");
    C_error_POD(ii)=page_norm(CROM-CFOM)/page_norm(CFOM);
    f1_error_POD(ii)=abs(f1ROM-f1FOM)/abs(f1FOM);
    f2_error_POD(ii)=abs(f2ROM-f2FOM)/abs(f2FOM);
end

%% save errors


error.E_error_POD = E_error_POD;
error.C_error_POD = C_error_POD;
error.f1_error_POD = f1_error_POD;
error.f2_error_POD = f2_error_POD;

error.E_error_OpInf = E_error_OpInf;
error.C_error_OpInf = C_error_OpInf;
error.f1_error_OpInf = f1_error_OpInf;
error.f2_error_OpInf = f2_error_OpInf;


%save( "./errors/error"+datestr(now, 'yyyy-mm-dd_HH-MM-SS')+".mat", '-struct', 'error', '-v7.3');

%% plot

ymin = 1e-6;

clf;
subplot(2,2,1)
semilogy(error.E_error_OpInf,'LineWidth',2)
legend
hold on
plot(error.E_error_POD,'LineWidth',2)
title("Expectation")
axis([1 rmax ymin 1])

subplot(2,2,2)
semilogy(error.C_error_OpInf,'LineWidth',2)
hold on
plot(error.C_error_POD,'LineWidth',2)
title("Covariance")
axis([1 rmax ymin 1])

subplot(2,2,3)
semilogy(error.f1_error_OpInf,'LineWidth',2)
hold on
plot(error.f1_error_POD,'LineWidth',2)
title("f1")
axis([1 rmax ymin 1])

subplot(2,2,4)
semilogy(error.f2_error_OpInf,'LineWidth',2)
hold on
plot(error.f2_error_POD,'LineWidth',2)
title("f2")
axis([1 rmax ymin 1])

%norm(EFOM(:,end))
%norm(CFOM(:,:,end))

%%
norm(EFOM,"fro")
page_norm(CFOM)

% for ii=1:2:2*rmax
% subplot(rmax,2,ii)
% plot(V(:,ii))
% hold on
% plot(V2(:,ii))
% subplot(rmax,2,ii+1)
% semilogy(abs(V(:,ii)-sign(V(1,ii)/V2(1,ii))*V2(:,ii)))
% end

function [y] = page_norm(A)
	y = sum(abs(A).^2,'all');
	y = sqrt(y);
end
