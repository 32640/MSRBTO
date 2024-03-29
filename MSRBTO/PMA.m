function [fvalx,Tmax01,Tmax02,Temfield,f0val,df0dx,df0dx2,fval,dfdx,dfdx2,xran,gnum,dgnum] = PMA(elecount_macro,nodecount_macro,elecount,nodecount,lx,ly,x,xPhys,xPhys_macro,xran,nelx_macro,nely_macro,penal_macro,nely,nelx,penal,rmin_macro,iK,jK,edofMat,iK_T,jK_T,edofMat_T,nu,hou,xishu,vxishu,betat,betatfz,n,ns,nf,nrd,gnum,dgnum,nfan,ft,xTilde,xTilde_macro,eta,beta,H,Hs,H_macro,Hs_macro,volfrac,volfrac_macro,m,ngd)
% xval=x_macro %宏观设计变量
% xran=mu %随机变量和模糊变量
% nelx %划分网格的数x方向
% nely %划分网格的数y方向
% penal %罚函数
% rmin %灵敏度过滤半径
% nu %泊松比
% hou %厚度
% xishu %约束前的系数，防止优化违反约束
% vxishu %约束前的系数，防止优化违反约束
[mu,sigma,type]=distribution;
[n]=length(xPhys_macro);  %即： n=nelx_macro*nely_macro
if ns>0         % ns=0; %设计变量中是随机变量参数的个数,即设计变量相对密度是与随机变量无关的
nsi=n-ns+1;
 for i=nsi:1:n
     mu(i-nsi+1)=xPhys_macro(i);
 end
else
end

[ng,nr]=size(xran);  %提取xran的维度大小，即 ng=4,nr=1 
gnum=gnum+1;
dgnum=dgnum+n;

for i=ngd+1:ng
 xxran=xran(i,:);
 itte=0;
 maxitte=1;
 precisionbeta=1;
 while (itte<maxitte&&precisionbeta>1e-6)
  itte=itte+1;
  xold=xxran;
  [~,~,~,~,~,~,~,g,gd,gx]=top(elecount_macro,nodecount_macro,elecount,nodecount,lx,ly,x,xPhys,xPhys_macro,xxran,nelx_macro,nely_macro,penal_macro,nely,nelx,penal,rmin_macro,nu,hou,iK,jK,edofMat,iK_T,jK_T,edofMat_T,ft,xTilde,xTilde_macro,eta,beta,H,Hs,H_macro,Hs_macro,volfrac,volfrac_macro,m);
  gnum=gnum+1;
  dgnum=dgnum+nr;
  gxx=gx(i,:);
  gdd=gd(i,:);  %位移约束灵敏度
  [u,gu]=transform(xxran,gxx,type,mu,sigma); %转变成标准正态分布
  %beta=(g(i)-gu*u')/norm(gu)
%   numfz=1:nf; %%%模糊变量数目 nf=0
  numpro=(nf+1):length(mu); %%%概率变量数目
%   p=nfan;
% u(numfz)=-betatfz*sign(gu(numfz)).*(abs(gu(numfz))/norm(gu(numfz),p/(p-1))).^(1/(p-1));
% beta1=(g(i)-gu(numpro)*u(numpro)')/norm(gu(numpro));
  u(numpro)=-betat*gu(numpro)/norm(gu(numpro));
  xxran=transforminv(u,type,mu,sigma,xxran); %从正态空间转化到原本的随机分布空间
  fvalx(i,1)=-g(i);     %混合可靠度指标的约束条件值
  betadx(i,:)=-gdd;
  precisionbeta=norm(xxran-xold)/norm(xxran);  %判别收敛的条件
 end
 dgnum=dgnum+(n-ns);
 dfdx(i,:)=betadx(i,:);
 xran(i,:)=xxran;
end

 for i=1:ngd
  xxran=mu;
  %% [f0val,df0dx,df0dx2,fval,dfdx,dgdx] = top(xval,xran,nelx,nely,penal,rmin,nu,hou,iK,jK,edofMat,W)
  [c,Tmax01,Tmax02,Temfield,f0val,df0dx,df0dx2,g,gd,gx]=top(elecount_macro,nodecount_macro,elecount,nodecount,lx,ly,x,xPhys,xPhys_macro,xxran,nelx_macro,nely_macro,penal_macro,nely,nelx,penal,rmin_macro,nu,hou,iK,jK,edofMat,iK_T,jK_T,edofMat_T,ft,xTilde,xTilde_macro,eta,beta,H,Hs,H_macro,Hs_macro,volfrac,volfrac_macro,m);
  gnum=gnum+1;
  dgnum=dgnum+nr;
  gdd=gd(i,:);        %位移约束灵敏度
  fvalx(i,1)=-g(i);   %混合可靠度指标的约束条件值
  dfdx(i,:)=-gdd;
  dgnum=dgnum+(n-ns);
 end

dfdx2=0*dfdx;  %原本此处为dfdx不是dfdx'，是我自己改的
fval=fvalx*vxishu;
dfdx=dfdx*vxishu;%原本此处为dfdx不是dfdx'，是我自己改的

end
