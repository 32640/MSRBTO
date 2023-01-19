function [KHe,FHe,KHe_T,FHe_T,H]=plani4e(ex,ey,ep,D,DDY_T)
%%  �ı��εȲ�Ԫ
hou=0.004; %������ȣ���������Ƴߴ�ʱ���ú�Ȳ��ɺ���
%% ���룺�ڵ�����ex ey���������ep�����Ծ���D
 t=ep(2);  ir=ep(3);  ngp=ir*ir;
%% ��ʼ��
KHe=zeros(8,8);
FHe=zeros(8,3);
KHe_T=zeros(4,4);
FHe_T=zeros(4,2);
H=zeros(8,4);
bb=12*10^(-6);
A=[bb bb 0];
%% ��˹��
[ gp, w] = GaussIntegration( ir );
wp=w(:,1).*w(:,2);
xsi=gp(:,1);  eta=gp(:,2);
r2=ngp*2;
%% �κ���
N(:,1)=(1-xsi).*(1-eta)/4;  N(:,2)=(1+xsi).*(1-eta)/4;
N(:,3)=(1+xsi).*(1+eta)/4;  N(:,4)=(1-xsi).*(1+eta)/4;

dNr(1:2:r2,1)=-(1-eta)/4;      dNr(1:2:r2,2)= (1-eta)/4;
dNr(1:2:r2,3)= (1+eta)/4;     dNr(1:2:r2,4)=-(1+eta)/4;
dNr(2:2:r2,1)=-(1-xsi)/4;      dNr(2:2:r2,2)=-(1+xsi)/4;
dNr(2:2:r2,3)= (1+xsi)/4;     dNr(2:2:r2,4)= (1-xsi)/4;
%% �ſ˱Ⱦ���
JT=dNr*[ex;ey]';
%% ���㣺��˹��ѭ��
for i=1:ngp
    indx=[ 2*i-1; 2*i ];
    detJ=det(JT(indx,:));
    if detJ<10*eps
        disp('Jacobideterminant equal or less than zero!')
    end
    JTinv=inv(JT(indx,:));
    dNx=JTinv*dNr(indx,:);
    
    B(1,1:2:8)  =  dNx(1,:);
    B(2,2:2:8)  =  dNx(2,:);
    B(3,1:2:8)  =  dNx(2,:);
    B(3,2:2:8)  =  dNx(1,:);
    B_T(1,1:4)  =  dNx(1,:);
    B_T(2,1:4)  =  dNx(2,:);
    KHe=KHe+B'*D*B*detJ*wp(i)*t;
    FHe=FHe-B'*D*detJ*wp(i)*t;
    KHe_T=KHe_T+B_T'*DDY_T*B_T*detJ*wp(i)*t;
    FHe_T=FHe_T-B_T'*DDY_T*detJ*wp(i)*t;
    H=H+B'*D*A'*N(i,:)*detJ*wp(i)*t;
end
H=H*hou;
KHe= KHe*hou;
KHe_T= KHe_T*hou;
FHe=FHe*hou;
FHe_T=FHe_T*hou;
%--------------------------end--------------------------------
