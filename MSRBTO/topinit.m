clc;clear;
format long
tic
global nelx nely nelx_macro nely_macro
nelx_macro=100;
nely_macro=130;
[edofMat,iK,jK]=Pre_FEA(nelx_macro,nely_macro);
[iK_T,jK_T,edofMat_T]=pre_heat(nelx_macro,nely_macro);
%% 微观模型参数
lx=1;
ly=1;
nelx=60;
nely=60;
elecount=nelx*nely;
nodecount=(nely+1)*(nelx+1);
dofcount=2*nodecount;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elecount_macro=nelx_macro*nely_macro;
nodecount_macro=(nely_macro+1)*(nelx_macro+1);
%% 材料参数
ptype=1;
thickness=1;
integrationRule=2;
ep=[ptype thickness integrationRule];
materialNum=1;
%% 优化模型
rmin=3.5;
rmin_macro=3.5;
volfrac=0.3;
volfrac_macro=0.4;
penal=3;
penal_macro=3;
ft=3;
%% INITIALIZE ITERATION 初始化迭代
x = repmat(volfrac,nely,nelx);
x_macro=repmat(volfrac_macro,nely_macro,nelx_macro);
%  for i = 1:nelx
%      for j = 1:nely
%          if sqrt((i-nelx/2-0.5)^2+(j-nely/2-0.5)^2) < min(nelx,nely)/3   %生成单胞，中心是圆，可以对此时的x值画图看出来；
%              x(j,i) = volfrac/4;
%          end
%      end
%  end
for i = 1:nelx
    for j = 1:nely
        if sqrt((i-nelx/4)^2+(j-nely/4)^2) < min(nelx,nely)/8
            x(j,i) = volfrac/4;
        end
        if sqrt((i-3*nelx/4)^2+(j-nely/4)^2) < min(nelx,nely)/8
            x(j,i) = volfrac/4;      
        end
        if sqrt((i-3*nelx/4)^2+(j-3*nely/4)^2) < min(nelx,nely)/8
            x(j,i) = volfrac/4;      
        end
        if sqrt((i-nelx/4)^2+(j-3*nely/4)^2) < min(nelx,nely)/8
            x(j,i) = volfrac/4;      
        end
    end
end
xPhys = x;
xPhys_macro=x_macro;
change = 1;
beta = 1;
loopbeta = 0;
eta=0.5;
[H,Hs,H_macro,Hs_macro]=filter0(nelx,nely,nelx_macro,nely_macro,rmin,rmin_macro);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  优化参数 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%
m = 4;
ngd=1;%确定性约束
n=elecount_macro+elecount; %设计变量个数，在此处表示第一个微观体积约束
nfan=10000;  
ns=0;   %设计变量中是随机变量参数的个数,即设计变量相对密度是与随机变量无关的
[mu,sigma,type]=distribution;
nf=0;  %%%模糊变量数目
nrd=length(mu)-nf;   %%%概率变量数目
betat=2; %%%混合可靠性指标下限
betatfz=1; %%%模糊目标可靠性，一般是1
xishu=1e0;
vxishu=1e0;
nu=0.3; %泊松比u
hou=0.004; %厚度h
xran=[mu;mu;mu;mu];  %长度，宽度，弹性模量，载荷  的均值。
c = 100000*ones(m,1);
d = zeros(m,1);
a0 = 1;
a = zeros(m,1);
iter = 0;
gnum=0;
dgnum=0;
move=0.3;
maxite=300;
itte = 0;
while itte < maxite && change > 0.01
  iter = iter+1;
  itte = itte+1;
   xval=[x(:); x_macro(:)];
  xmax = min(1,xval+move);
  xmin = max(0.001,xval-move);
  
    if itte==1
        low=0.001*ones(n,1);
        upp=1*ones(n,1);
        xold1=xval;
        xold2=xval;
    end
   %% 均匀化分析  （光滑处理）
    if ft == 1
        xPhys = x;
        xPhys_macro = x_macro;
    elseif ft == 2
        xPhys(:) = (H*x(:))./Hs;
        xPhys_macro(:) = (H_macro*x_macro(:))./Hs_macro;
    elseif ft == 3
        xTilde = x;
        xTilde(:) = (H*x(:))./Hs;
        xPhys(:) = Heaviside(xTilde,beta,eta);
        xTilde_macro = x_macro;
        xTilde_macro(:) = (H_macro*x_macro(:))./Hs_macro;
        xPhys_macro(:) = Heaviside(xTilde_macro,beta,eta);
    end
  [fvalx,Tmax01,Tmax02,Temfield,f0val,df0dx,df0dx2,fval,dfdx,dfdx2,xran,gnum,dgnum] = PMA(elecount_macro,nodecount_macro,elecount,nodecount,lx,ly,x,xPhys,xPhys_macro,xran,nelx_macro,nely_macro,penal_macro,nely,nelx,penal,rmin_macro,iK,jK,edofMat,iK_T,jK_T,edofMat_T,nu,hou,xishu,vxishu,betat,betatfz,n,ns,nf,nrd,gnum,dgnum,nfan,ft,xTilde,xTilde_macro,eta,beta,H,Hs,H_macro,Hs_macro,volfrac,volfrac_macro,m,ngd);
  [xmma,ymma,zmma,lam,xsi,eta1,mu,zet,s,low,upp] = mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2, ...
  f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a,c,d);

    xold2=xold1;
    xold1=xval;       
 change = max(abs(xmma-xval));
  xval=xmma;
  x(:)=xval(1:elecount);
  x_macro(:)=xval(elecount+1:end);
  ix=51:nely_macro;
  iy1=1:nelx_macro/4 ;
  iy2=3*nelx_macro/4+1:nelx_macro ;
 for i=ix
    x_macro(i,iy1)=1e-9;
    x_macro(i,iy2)=1e-9;
    Temfield(i+1,iy1)=1e-9;
    Temfield(i,iy2)=1e-9;
 end
 %% UPDATE HEAVISIDE REGULARIZATION PARAMETER
    if ft == 3 && beta < 512 && (loopbeta >= 50 || change <= 0.01)
        beta = 2*(beta);
        loopbeta = 0;
        change = 1;
        fprintf('Parameter beta increased to %g.\n',beta);
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  fprintf(' It.:%3i   v_macro.:%11.4f   v.:%5.3f   fval.:%5.3f\n',itte,mean(xPhys_macro(:)), mean(xPhys(:)),fvalx/100);
 fprintf(' It.:%3i   obj.:%11.4f   v.:%5.3f   C_macro.:%5.3f  Tmax1.:%7.3f   Tmax2.:%7.3f\n',itte,f0val/elecount_macro, fval(1), fval(2),fval(3),fval(4));
    xxx=[xPhys xPhys;xPhys xPhys];
    title('beta=2')
    figure(1)
    subplot(2,2,1); colormap(gray); imagesc(1-xxx); caxis([0 1]); axis equal; axis off; drawnow;
    subplot(2,2,[3 4]); colormap(gray); imagesc(1-xPhys_macro); caxis([0 1]); axis equal; axis off; drawnow;
    subplot(2,2,2); colormap(jet); imagesc(Temfield); axis equal; axis tight; axis off;
    iteration(itte)=f0val/elecount_macro;
end
save('beat2.mat','xval')
toc
