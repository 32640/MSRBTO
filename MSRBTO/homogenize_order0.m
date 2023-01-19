function [ DH0,DH0DX,H,DH0_T,DH0DX_T,H_T] = homogenize_order0(lx,ly,x,DD,ep,penal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lx        = Unit cell length in x-direction. x�����ϵ����ĳ���
% ly        = Unit cell length in y-direction. y�����ϵ����ĳ���
% x         = Material indicator matrix. Size used to determine nelx/nely
% DD      =matrial stiffness for both materials
% ep       =analysis type��thickness; integration rule;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%�ȴ����ĵ��Ծ���
k=56;  %�ȴ���ϵ��
DD_T=[k 0;0 k];
%% INITIALIZE
% Deduce discretization �Ƶ�����ɢ��
[nely, nelx] = size(x);
dx = lx/nelx;  %1x=0.2
dy = ly/nely;  %1y=1
%  ��Ԫ��Ŀ
elecount=nelx*nely;
%  �ڵ���Ŀ
nodecount=(nely+1)*(nelx+1);
%  �ڵ����ɶ� 
nodedofdegree=2;
nodedofdegree_T=1;
% ģ�����ɶ���Ŀ
dofcount=nodedofdegree*nodecount;     %��ʵ���� dofcount=2*��nely+1��*��nelx+1)
dofcount_T=nodedofdegree_T*nodecount; 
%% ���񻮷� ����֮ǰ������������ͬ��
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
clear edofVec
%% ���������񻮷�
edofVec_T = reshape(nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat_T = repmat(edofVec_T,1,4)+repmat([0 nely+[1 0] -1],nelx*nely,1);
clear nodenrs edofVec_T
%% ������������Ա߽�����
alldofs     = [1:2*(nely+1)*(nelx+1)];
fixeddofs11 = [1:1:2*(nely+1)];%���
fixeddofs12 = [2*nelx*(nely+1)+1:1:2*(nelx+1)*(nely+1)];%�ұ�
fixeddofs21 = union([1:2*(nely+1):2*nelx*(nely+1)+1],[2:2*(nely+1):2*nelx*(nely+1)+2]);%�ϱ�
fixeddofs22 = union([2*(nely+1)-1:2*(nely+1):2*(nelx+1)*(nely+1)-1],[2*(nely+1):2*(nely+1):2*(nelx+1)*(nely+1)]);%�±�
freedofs    = setdiff(alldofs,fixeddofs12);
freedofs    = setdiff(freedofs,fixeddofs22);
freedofs    = setdiff(freedofs,[1 2]);%��Ҫע����ǣ����ɶ�1,2��Լ��ס����ԭ��Ϊ[1 2 4]
%% ������������Ա߽�����
alldofs_T     = [1:(nely+1)*(nelx+1)];
fixeddofs11_T = [1:1:(nely+1)];                       %��
fixeddofs12_T = [nelx*(nely+1)+1:1:(nelx+1)*(nely+1)];%��
fixeddofs21_T = [1:(nely+1):nelx*(nely+1)+1];         %��
fixeddofs22_T = [(nely+1):(nely+1):(nelx+1)*(nely+1)];%��
freedofs_T    = setdiff(alldofs_T,fixeddofs12_T);
freedofs_T    = setdiff(freedofs_T,fixeddofs22_T);
freedofs_T    = setdiff(freedofs_T,[1]);
%% �Ⱥ��������0�׼���
KH = sparse(dofcount, dofcount);
FH = sparse(dofcount,3);
H=zeros(dofcount,3);
KH_T = sparse(dofcount_T, dofcount_T);
FH_T = sparse(dofcount_T,2);
H_T=zeros(dofcount_T,2);
%% ��װ
for i=1:elecount
    % ��õ�Ԫ�����󣡣�����
    DDY=(0.001+x(i).^penal*(1-0.001))*DD{1,1};
    DDY_T=(0.001+x(i).^penal*(1-0.001))*DD_T;
    %�ڵ�����
    ex=[0 dx dx 0];
    ey=[0 0 dy dy];
    % ���㵥Ԫ�ն�����غ���
    [KHe,FHe,KHe_T,FHe_T,~]=plani4e(ex,ey,ep,DDY,DDY_T); %��˹�����
    edof=edofMat(i,:);
    edof_T=edofMat_T(i,:);
    % ��װ
    KH(edof,edof) = KH(edof,edof) + KHe;
    FH(edof,:) = FH(edof,:) +FHe;
   
    KH_T( edof_T, edof_T) = KH_T( edof_T, edof_T) + KHe_T;
    FH_T( edof_T,:) = FH_T( edof_T,:) + FHe_T;
end
%% ��������⣬��Ŀ���ǽ������Ա߽��������룩
KH(fixeddofs11,alldofs)=KH(fixeddofs11,alldofs)+KH(fixeddofs12,alldofs);
KH(alldofs,fixeddofs11)=KH(alldofs,fixeddofs11)+KH(alldofs,fixeddofs12);
KH(fixeddofs21,alldofs)=KH(fixeddofs21,alldofs)+KH(fixeddofs22,alldofs);
KH(alldofs,fixeddofs21)=KH(alldofs,fixeddofs21)+KH(alldofs,fixeddofs22);
FH(fixeddofs11,:)=FH(fixeddofs11,:)+FH(fixeddofs12,:);
FH(fixeddofs21,:)=FH(fixeddofs21,:)+FH(fixeddofs22,:);
H(freedofs,:) = KH(freedofs,freedofs) \ FH(freedofs,:);
H(fixeddofs12,:)= H(fixeddofs11,:);
H(fixeddofs22,:)= H(fixeddofs21,:);
H([1 2],:)=0;%ԭ��ΪH([1 2 4],:)=0
%% ��������⣬
KH_T(fixeddofs11_T,alldofs_T)=KH_T(fixeddofs11_T,alldofs_T)+KH_T(fixeddofs12_T,alldofs_T);  %�ж�Ӧ����
KH_T(alldofs_T,fixeddofs11_T)=KH_T(alldofs_T,fixeddofs11_T)+KH_T(alldofs_T,fixeddofs12_T);  %�ж�Ӧ����
FH_T(fixeddofs11_T,:)=FH_T(fixeddofs11_T,:)+FH_T(fixeddofs12_T,:);  %�������������Ա߽�, Fֻ���н����ۼ�
KH_T(fixeddofs21_T,alldofs_T)=KH_T(fixeddofs21_T,alldofs_T)+KH_T(fixeddofs22_T,alldofs_T);  %�ж�Ӧ����
KH_T(alldofs_T,fixeddofs21_T)=KH_T(alldofs_T,fixeddofs21_T)+KH_T(alldofs_T,fixeddofs22_T);  %�ж�Ӧ����
FH_T(fixeddofs21_T,:)=FH_T(fixeddofs21_T,:)+FH_T(fixeddofs22_T,:);  %
H_T(freedofs_T,:) = KH_T(freedofs_T,freedofs_T) \ FH_T(freedofs_T,:);
H_T(fixeddofs12_T,:)= H_T(fixeddofs11_T,:);
H_T(fixeddofs22_T,:)= H_T(fixeddofs21_T,:);
H_T([1],:)=0;

%% 0�׵�Ч
DH0=zeros(3,3);
DH0_T=zeros(2,2);
for i=1:elecount
    edof=edofMat(i,:);
    edof_T=edofMat_T(i,:);
    He=H(edof,:);
    He_T=H_T(edof_T,:);
    % ��õ�Ԫ�����󣡣�����
    DDY=(0.001+x(i).^penal*(1-0.001))*DD{1,1};
    DDY_T=(0.001+x(i).^penal*(1-0.001))*DD_T;
    [C0e,C0e_T]=planeAssemble1(ex,ey,ep,DDY,He,DDY_T,He_T); %ר������Ч���Ծ���
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
    % ��õ�Ԫ�����󣡣�����
    DDYDX=penal*x(i).^(penal-1)*(1-0.001)*DD{1,1};
    DDYDX_T=penal*x(i).^(penal-1)*(1-0.001)*DD_T;
    [C0eDX,C0eDX_T]=planeAssemble1(ex,ey,ep,DDYDX,He,DDYDX_T,He_T);%3x3�ľ���
    DH0DX{i}=C0eDX/lx/ly;    
    DH0DX_T{i}=C0eDX_T/lx/ly; %��ʾDH0��ÿһ����Ԫ��Ʊ�����������
end
end