function xxran=transforminv(u,type,mu,sigma,x)
[nn,nr]=size(u);
for i=1:1:nr
    t=type(i);
    s=sigma(i);
    if t==1
        xxran(i)=sigma(i)*u(i)+mu(i);
    elseif t==2
        sln=sqrt(log(1+(sigma(i)/mu(i))^2));
        mln=log(mu(i))-sln^2/2;
        cdfx=logncdf(x(i),mln,sln);
        pdfx=lognpdf(x(i),mln,sln);
        nc=norminv(cdfx);
        sigmaxx=normpdf(nc)/pdfx;
        muxx=x(i)-nc*sigmaxx;
        xxran(i)=u(i)*sigmaxx+muxx
    elseif t==3
        aev=sqrt(6)*sigma(i)/pi;
        uev=-psi(1)*aev-mu(i);
        cdfx=1-evcdf(-x(i),uev,aev);
        pdfx=evpdf(-x(i),uev,aev);
        nc=norminv(cdfx);
        sigmaxx=normpdf(nc)/pdfx;
        muxx=x(i)-nc*sigmaxx;
        xxran(i)=u(i)*sigmaxx+muxx
    end
end
