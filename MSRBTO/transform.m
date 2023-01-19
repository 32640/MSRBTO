function [u,gu]=transform(xxran,gx,type,mu,sigma)
[nn,nr]=size(xxran);
u=zeros(1,nr);
for i=1:1:nr
    t=type(i); %实际上t一直等于1；
    if t==1
        u(i)=(xxran(i)-mu(i))/sigma(i);
        gu(i)=gx(i)*sigma(i);  %sigma的作用？
    elseif t==2
        sln=sqrt(log(1+(sigma(i)/mu(i))^2));
        mln=log(mu(i))-sln^2/2;
        cdfx=logncdf(xxran(i),mln,sln);
        pdfx=lognpdf(xxran(i),mln,sln);
        nc=norminv(cdfx);
        sigmaxx=normpdf(nc)/pdfx;
        muxx=xxran(i)-nc*sigmaxx;
        u(i)=(xxran(i)-muxx)/sigmaxx;
        gu(i)=gx(i)*sigmaxx;
    elseif t==3
        aev=sqrt(6)*sigma(i)/pi;
        uev=-psi(1)*aev-mu(i);
        cdfx=1-evcdf(-xxran(i),uev,aev);
        pdfx=evpdf(-xxran(i),uev,aev);
        nc=norminv(cdfx);
        sigmaxx=normpdf(nc)/pdfx;
        muxx=x-nc*sigmaxx;
        u(i)=(xxran(i)-muxx)/sigmaxx;
        gu(i)=gx(i)*sigmaxx;
    end
end
