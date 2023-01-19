function [fvalx,Tmax01,Tmax02,Temfield,f0val,df0dx,df0dx2,fval,dfdx,dfdx2,xran,gnum,dgnum] = PMA(elecount_macro,nodecount_macro,elecount,nodecount,lx,ly,x,xPhys,xPhys_macro,xran,nelx_macro,nely_macro,penal_macro,nely,nelx,penal,rmin_macro,iK,jK,edofMat,iK_T,jK_T,edofMat_T,nu,hou,xishu,vxishu,betat,betatfz,n,ns,nf,nrd,gnum,dgnum,nfan,ft,xTilde,xTilde_macro,eta,beta,H,Hs,H_macro,Hs_macro,volfrac,volfrac_macro,m,ngd)
% xval=x_macro %�����Ʊ���
% xran=mu %���������ģ������
% nelx %�����������x����
% nely %�����������y����
% penal %������
% rmin %�����ȹ��˰뾶
% nu %���ɱ�
% hou %���
% xishu %Լ��ǰ��ϵ������ֹ�Ż�Υ��Լ��
% vxishu %Լ��ǰ��ϵ������ֹ�Ż�Υ��Լ��
[mu,sigma,type]=distribution;
[n]=length(xPhys_macro);  %���� n=nelx_macro*nely_macro
if ns>0         % ns=0; %��Ʊ�������������������ĸ���,����Ʊ�������ܶ�������������޹ص�
nsi=n-ns+1;
 for i=nsi:1:n
     mu(i-nsi+1)=xPhys_macro(i);
 end
else
end

[ng,nr]=size(xran);  %��ȡxran��ά�ȴ�С���� ng=4,nr=1 
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
  gdd=gd(i,:);  %λ��Լ��������
  [u,gu]=transform(xxran,gxx,type,mu,sigma); %ת��ɱ�׼��̬�ֲ�
  %beta=(g(i)-gu*u')/norm(gu)
%   numfz=1:nf; %%%ģ��������Ŀ nf=0
  numpro=(nf+1):length(mu); %%%���ʱ�����Ŀ
%   p=nfan;
% u(numfz)=-betatfz*sign(gu(numfz)).*(abs(gu(numfz))/norm(gu(numfz),p/(p-1))).^(1/(p-1));
% beta1=(g(i)-gu(numpro)*u(numpro)')/norm(gu(numpro));
  u(numpro)=-betat*gu(numpro)/norm(gu(numpro));
  xxran=transforminv(u,type,mu,sigma,xxran); %����̬�ռ�ת����ԭ��������ֲ��ռ�
  fvalx(i,1)=-g(i);     %��Ͽɿ���ָ���Լ������ֵ
  betadx(i,:)=-gdd;
  precisionbeta=norm(xxran-xold)/norm(xxran);  %�б�����������
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
  gdd=gd(i,:);        %λ��Լ��������
  fvalx(i,1)=-g(i);   %��Ͽɿ���ָ���Լ������ֵ
  dfdx(i,:)=-gdd;
  dgnum=dgnum+(n-ns);
 end

dfdx2=0*dfdx;  %ԭ���˴�Ϊdfdx����dfdx'�������Լ��ĵ�
fval=fvalx*vxishu;
dfdx=dfdx*vxishu;%ԭ���˴�Ϊdfdx����dfdx'�������Լ��ĵ�

end