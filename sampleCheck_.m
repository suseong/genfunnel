function m=sampleCheck_(x,V,Vdot,rho,rhodot)
  if (deg(V,x)>2) error('only checks quadratics'); end
  
  n=length(x);
  K=100;
  X = randn(n,K);
  X = X./repmat(sqrt(sum(X.^2,1)),n,1);

  H = double(0.5*diff(diff(V,x)',x));
  b = -0.5*(H\double(subs(diff(V,x),x,0*x)'));

  try 
  X = repmat(b,1,K) + (H/(double(rho-subs(V,x,b))+eps))^(-1/2)*X;
  catch
    keyboard;
  end
  m=max(double(msubs(Vdot,x,X))) - rhodot;
  if (m>0)
    warning('found a positive Vdot');
  end
end