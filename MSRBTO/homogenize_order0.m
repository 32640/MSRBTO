function [ DH0,DH0DX,H,DH0_T,DH0DX_T,H_T] = homogenize_order0(lx,ly,x,DD,ep,penal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lx        = Unit cell length in x-direction. x方向上单胞的长度
% ly        = Unit cell length in y-direction. y方向上单胞的长度
% x         = Material indicator matrix. Size used to determine nelx/nely
% DD      =matrial stiffness for both materials
% ep       =analysis type；thickness; integration rule;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%热传导的弹性矩阵
k=56;  %热传导系数
DD_T=[k 0;0 k];
%% INITIALIZE
% Deduce discretization 推导出离散化
[nely, nelx] = size(x);
dx = lx/nelx;  %1x=0.2
dy = ly/nely;  %1y=1
%  单元数目
elecount=nelx*nely;
%  节点数目
nodecount=(nely+1)*(nelx+1);
%  节点自由度 
nodedofdegree=2;
nodedofdegree_T=1;
% 模型自由度数目
dofcount=nodedofdegree*nodecount;     %其实就是 dofcount=2*（nely+1）*（nelx+1)
dofcount_T=nodedofdegree_T*nodecount; 
%% 网格划分 （和之前所理解的意义相同）
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
clear edofVec
%% 热问题网格划分
edofVec_T = reshape(nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat_T = repmat(edofVec_T,1,4)+repmat([0 nely+[1 0] -1],nelx*nely,1);
clear nodenrs edofVec_T
%% 力问题的周期性边界条件
alldofs     = [1:2*(nely+1)*(nelx+1)];
fixeddofs11 = [1:1:2*(nely+1)];%左边
fixeddofs12 = [2*nelx*(nely+1)+1:1:2*(nelx+1)*(nely+1)];%右边
fixeddofs21 = union([1:2*(nely+1):2*nelx*(nely+1)+1],[2:2*(nely+1):2*nelx*(nely+1)+2]);%上边
fixeddofs22 = union([2*(nely+1)-1:2*(nely+1):2*(nelx+1)*(nely+1)-1],[2*(nely+1):2*(nely+1):2*(nelx+1)*(nely+1)]);%下边
freedofs    = setdiff(alldofs,fixeddofs12);
freedofs    = setdiff(freedofs,fixeddofs22);
freedofs    = setdiff(freedofs,[1 2]);%需要注意的是，自由度1,2被约束住！，原本为[1 2 4]
%% 热问题的周期性边界条件
alldofs_T     = [1:(nely+1)*(nelx+1)];
fixeddofs11_T = [1:1:(nely+1)];                       %左
fixeddofs12_T = [nelx*(nely+1)+1:1:(nelx+1)*(nely+1)];%右
fixeddofs21_T = [1:(nely+1):nelx*(nely+1)+1];         %上
fixeddofs22_T = [(nely+1):(nely+1):(nelx+1)*(nely+1)];%下
freedofs_T    = setdiff(alldofs_T,fixeddofs12_T);
freedofs_T    = setdiff(freedofs_T,fixeddofs22_T);
freedofs_T    = setdiff(freedofs_T,[1]);
%% 热和力问题的0阶计算
KH = sparse(dofcount, dofcount);
FH = sparse(dofcount,3);
H=zeros(dofcount,3);
KH_T = sparse(dofcount_T, dofcount_T);
FH_T = sparse(dofcount_T,2);
H_T=zeros(dofcount_T,2);
%% 组装
for i=1:elecount
    % 获得单元弹性阵！！！！
    DDY=(0.001+x(i).^penal*(1-0.001))*DD{1,1};
    DDY_T=(0.001+x(i).^penal*(1-0.001))*DD_T;
    %节点坐标
    ex=[0 dx dx 0];
    ey=[0 0 dy dy];
    % 计算单元刚度阵和载荷阵
    [KHe,FHe,KHe_T,FHe_T,~]=plani4e(ex,ey,ep,DDY,DDY_T); %高斯型求积
    edof=edofMat(i,:);
    edof_T=edofMat_T(i,:);
    % 组装
    KH(edof,edof) = KH(edof,edof) + KHe;
    FH(edof,:) = FH(edof,:) +FHe;
   
    KH_T( edof_T, edof_T) = KH_T( edof_T, edof_T) + KHe_T;
    FH_T( edof_T,:) = FH_T( edof_T,:) + FHe_T;
end
%% 力问题求解，（目的是将周期性边界条件加入）
KH(fixeddofs11,alldofs)=KH(fixeddofs11,alldofs)+KH(fixeddofs12,alldofs);
KH(alldofs,fixeddofs11)=KH(alldofs,fixeddofs11)+KH(alldofs,fixeddofs12);
KH(fixeddofs21,alldofs)=KH(fixeddofs21,alldofs)+KH(fixeddofs22,alldofs);
KH(alldofs,fixeddofs21)=KH(alldofs,fixeddofs21)+KH(alldofs,fixeddofs22);
FH(fixeddofs11,:)=FH(fixeddofs11,:)+FH(fixeddofs12,:);
FH(fixeddofs21,:)=FH(fixeddofs21,:)+FH(fixeddofs22,:);
H(freedofs,:) = KH(freedofs,freedofs) \ FH(freedofs,:);
H(fixeddofs12,:)= H(fixeddofs11,:);
H(fixeddofs22,:)= H(fixeddofs21,:);
H([1 2],:)=0;%原本为H([1 2 4],:)=0
%% 热问题求解，
KH_T(fixeddofs11_T,alldofs_T)=KH_T(fixeddofs11_T,alldofs_T)+KH_T(fixeddofs12_T,alldofs_T);  %行对应叠加
KH_T(alldofs_T,fixeddofs11_T)=KH_T(alldofs_T,fixeddofs11_T)+KH_T(alldofs_T,fixeddofs12_T);  %列对应叠加
FH_T(fixeddofs11_T,:)=FH_T(fixeddofs11_T,:)+FH_T(fixeddofs12_T,:);  %左右两侧周期性边界, F只对行进行累加
KH_T(fixeddofs21_T,alldofs_T)=KH_T(fixeddofs21_T,alldofs_T)+KH_T(fixeddofs22_T,alldofs_T);  %行对应叠加
KH_T(alldofs_T,fixeddofs21_T)=KH_T(alldofs_T,fixeddofs21_T)+KH_T(alldofs_T,fixeddofs22_T);  %列对应叠加
FH_T(fixeddofs21_T,:)=FH_T(fixeddofs21_T,:)+FH_T(fixeddofs22_T,:);  %
H_T(freedofs_T,:) = KH_T(freedofs_T,freedofs_T) \ FH_T(freedofs_T,:);
H_T(fixeddofs12_T,:)= H_T(fixeddofs11_T,:);
H_T(fixeddofs22_T,:)= H_T(fixeddofs21_T,:);
H_T([1],:)=0;

%% 0阶等效
DH0=zeros(3,3);
DH0_T=zeros(2,2);
for i=1:elecount
    edof=edofMat(i,:);
    edof_T=edofMat_T(i,:);
    He=H(edof,:);
    He_T=H_T(edof_T,:);
    % 获得单元弹性阵！！！！
    DDY=(0.001+x(i).^penal*(1-0.001))*DD{1,1};
    DDY_T=(0.001+x(i).^penal*(1-0.001))*DD_T;
    [C0e,C0e_T]=planeAssemble1(ex,ey,ep,DDY,He,DDY_T,He_T); %专门求解等效弹性矩阵
    DH0_T=DH0_T+C0e_T;
    DH0=DH0+C0e;
end
DH0=DH0/lx/ly;
DH0_T=DH0_T/lx/ly;
DH0DX=cell(nely,nelx);
DH0DX_T=cell(nely,nelx);
for i=1:elecount
    edof=edofMat(i,:);
    edof_T=edofMat_T(i,:);
    He=H(edof,:);
    He_T=H_T(edof_T,:);
    % 获得单元弹性阵！！！！
    DDYDX=penal*x(i).^(penal-1)*(1-0.001)*DD{1,1};
    DDYDX_T=penal*x(i).^(penal-1)*(1-0.001)*DD_T;
    [C0eDX,C0eDX_T]=planeAssemble1(ex,ey,ep,DDYDX,He,DDYDX_T,He_T);%3x3的矩阵
    DH0DX{i}=C0eDX/lx/ly;    
    DH0DX_T{i}=C0eDX_T/lx/ly; %表示DH0对每一个单元设计变量的灵敏度
end
end