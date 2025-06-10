function [E,C,f1,f2] = estimate(f,Vr,x0,u,L)
  
  % sample model
  %X = pagemtimes(Vr,stepSDE(f,x0,u,L));
	X = stepSDE(f,x0,u,L);

  % estimate quantities
  E = Vr*reshape(mean(X,2),size(X,[1 3]));
  C = pagemtimes(pagemtimes(Vr,page_cov(X,true)),Vr');
	VrX = Vr*X(:,:,end);
	f1 = mean(vecnorm(VrX.^2));
	f2 = mean(VrX.^3./exp(VrX),'all');

end

