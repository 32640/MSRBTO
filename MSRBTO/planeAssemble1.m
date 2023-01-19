function [C0e,C0e_T]=planeAssemble1(ex,ey,ep,D,He,DDY_T,He_T)
%%  四边形等参元
%% 输入：节点坐标ex ey，计算参数ep，弹性矩阵D
 t=ep(2);  ir=ep(3);  ngp=ir*ir;
%% 初始化
C0e=zeros(3,3);
C0e1=zeros(3,3);
C0e_T=zeros(2,2);
%% 高斯点
[ gp, w] = GaussIntegration( ir );
wp=w(:,1).*w(:,2);
xsi=gp(:,1);  eta=gp(:,2);  r2=ngp*2;
%% 形函数
N(:,1)=(1-xsi).*(1-eta)/4;  N(:,2)=(1+xsi).*(1-eta)/4;
N(:,3)=(1+xsi).*(1+eta)/4;  N(:,4)=(1-xsi).*(1+eta)/4;

dNr(1:2:r2,1)=-(1-eta)/4;      dNr(1:2:r2,2)= (1-eta)/4;
dNr(1:2:r2,3)= (1+eta)/4;      dNr(1:2:r2,4)=-(1+eta)/4;
dNr(2:2:r2,1)=-(1-xsi)/4;      dNr(2:2:r2,2)=-(1+xsi)/4;
dNr(2:2:r2,3)= (1+xsi)/4;      dNr(2:2:r2,4)= (1-xsi)/4;
%% 雅克比矩阵
JT=dNr*[ex;ey]';
I=eye(3);
I_T=eye(2);
%% 计算：高斯点循环
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
      C0e1=C0e1+(-D*B*He+D)*detJ*wp(i)*t;
    C0e=C0e+(I+B*He)'*D*(I+B*He)*detJ*wp(i)*t;
    C0e_T=C0e_T+(I_T+B_T*He_T)'*DDY_T*(I_T+B_T*He_T)*detJ*wp(i)*t;
end
%--------------------------end--------------------------------