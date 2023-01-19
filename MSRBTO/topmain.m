
%
itte = 0;
while itte < maxite
  iter = iter+1;
  itte = itte+1;

  [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
  mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2, ...
  f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a,c,d);

  xold2 = xold1;
  xold1 = xval;
  xval = xmma;

  [f0val,df0dx,df0dx2,fval,dfdx,dfdx2] = three(xval);
  outvector = [iter f0val fval' xval']';
end
