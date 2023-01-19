function [Ta01,Ta04,Ti1,Ti2,T,freedofs]=tem1_field(Kth,nelx_macro,nely_macro)
nelx=nelx_macro;nely=nely_macro;
q=0.01;
topth = [1:nely+1:(nely+1)*nelx+1];
bottomth =[ nely+1:nely+1:(nely+1)*(nelx+1)];
bottomth1 =[ (nely+1)*nelx/4:nely+1:(nely+1)*3*nelx/4];
leftth = [1:nely+1];
rightth = [(nely+1)*nelx+1:(nely+1)*(nelx+1)];
Ta01=[(nely+1)*nelx/4+25+1];%,(nely+1)*nelx/5+nely/2,(nely+1)*(nelx/5-1)+nely/2+1,(nely+1)*(nelx/5-1)+nely/2];
% Ta02=[(nely+1)*2*nelx/5+nely/2+1,(nely+1)*2*nelx/5+nely/2,(nely+1)*(2*nelx/5-1)+nely/2+1,(nely+1)*(2*nelx/5-1)+nely/2];
% Ta03=[(nely+1)*3*nelx/5+nely/2+1,(nely+1)*3*nelx/5+nely/2,(nely+1)*(3*nelx/5-1)+nely/2+1,(nely+1)*(3*nelx/5-1)+nely/2];
Ta04=[(nely+1)*3*nelx/4+25+1];%,(nely+1)*4*nelx/5+nely/2,(nely+1)*(4*nelx/5-1)+nely/2+1,(nely+1)*(4*nelx/5-1)+nely/2];

Ta011=[(nely+1)*nelx/5+nely/2+1,(nely+1)*nelx/5+nely/2,(nely+1)*(nelx/5-1)+nely/2+1,(nely+1)*(nelx/5-1)+nely/2];
Ta044=[(nely+1)*4*nelx/5+nely/2+1,(nely+1)*4*nelx/5+nely/2,(nely+1)*(4*nelx/5-1)+nely/2+1,(nely+1)*(4*nelx/5-1)+nely/2];

Ta=[ Ta01 Ta04 ];
F= sparse(Ta,1,q,(nely+1)*(nelx+1),1);
Fth01= sparse(Ta01,1,1,(nely+1)*(nelx+1),1);
Fth02= sparse(Ta04,1,1,(nely+1)*(nelx+1),1);
U = zeros((nely+1)*(nelx+1),1);
Ti1=zeros((nely+1)*(nelx+1),1);
Ti2=zeros((nely+1)*(nelx+1),1);
fixeddofs =[bottomth1];
alldofs = [1:(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
U(freedofs) = Kth(freedofs,freedofs)\F(freedofs);
Ti1(freedofs)=-Kth(freedofs,freedofs)\Fth01(freedofs);
Ti2(freedofs)=-Kth(freedofs,freedofs)\Fth02(freedofs);
T=U;