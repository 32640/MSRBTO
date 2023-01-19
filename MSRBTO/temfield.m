 

function [Ta01,Ta04, T, Ti1,Ti2, Kth0]=temfield(Kth,Kth1,nelx_macro,nely_macro)
nelx=nelx_macro;nely=nely_macro;
q=0.05; %热源
% 第一列
%  Ta11=[(nely+1)*nelx/5+nely/4+1];%,(nely+1)*nelx/5+nely/4,(nely+1)*(nelx/5-1)+nely/4+1,(nely+1)*(nelx/5-1)+nely/4];
%  Ta12=[(nely+1)*nelx/5+3*nely/4+1];%,(nely+1)*nelx/5+3*nely/4,(nely+1)*(nelx/5-1)+3*nely/4+1,(nely+1)*(nelx/5-1)+3*nely/4];
% %第二列
%  Ta21=[(nely+1)*2*nelx/5+nely/4+1];%,(nely+1)*2*nelx/5+nely/4,(nely+1)*(2*nelx/5-1)+nely/4+1,(nely+1)*(2*nelx/5-1)+nely/4];
%  Ta22=[(nely+1)*2*nelx/5+3*nely/4+1];%,(nely+1)*2*nelx/5+3*nely/4,(nely+1)*(2*nelx/5-1)+3*nely/4+1,(nely+1)*(2*nelx/5-1)+3*nely/4];
%  %第三列
%  Ta31=[(nely+1)*3*nelx/5+nely/4+1];%,(nely+1)*3*nelx/5+nely/4,(nely+1)*(3*nelx/5-1)+nely/4+1,(nely+1)*(3*nelx/5-1)+nely/4];
%  Ta32=[(nely+1)*3*nelx/5+3*nely/4+1];%,(nely+1)*3*nelx/5+3*nely/4,(nely+1)*(3*nelx/5-1)+3*nely/4+1,(nely+1)*(3*nelx/5-1)+3*nely/4];
%  %第四列
%  Ta41=[(nely+1)*4*nelx/5+nely/4+1];%,(nely+1)*4*nelx/5+nely/4,(nely+1)*(4*nelx/5-1)+nely/4+1,(nely+1)*(4*nelx/5-1)+nely/4];
%  Ta42=[(nely+1)*4*nelx/5+3*nely/4+1];%,(nely+1)*4*nelx/5+3*nely/4,(nely+1)*(4*nelx/5-1)+3*nely/4+1,(nely+1)*(4*nelx/5-1)+3*nely/4];
% 
% %中间行
Ta01=[(nely+1)*nelx/4+nely/3];%,(nely+1)*nelx/5+nely/2,(nely+1)*(nelx/5-1)+nely/2+1,(nely+1)*(nelx/5-1)+nely/2];
% Ta02=[(nely+1)*2*nelx/5+nely/2+1,(nely+1)*2*nelx/5+nely/2,(nely+1)*(2*nelx/5-1)+nely/2+1,(nely+1)*(2*nelx/5-1)+nely/2];
% Ta03=[(nely+1)*3*nelx/5+nely/2+1,(nely+1)*3*nelx/5+nely/2,(nely+1)*(3*nelx/5-1)+nely/2+1,(nely+1)*(3*nelx/5-1)+nely/2];
Ta04=[(nely+1)*3*nelx/4+nely/3];%,(nely+1)*4*nelx/5+nely/2,(nely+1)*(4*nelx/5-1)+nely/2+1,(nely+1)*(4*nelx/5-1)+nely/2];
 

% Ta=[ Ta01 Ta02 Ta03 Ta04 Ta11 Ta12 Ta21 Ta22 Ta31 Ta32 Ta41 Ta42];
Ta=[ Ta01 Ta04 ];
Fth= sparse(Ta,1,q,(nely+1)*(nelx+1),1);  %两端固支梁，左右两点
Fth01= sparse(Ta01,1,1,(nely+1)*(nelx+1),1);
Fth02= sparse(Ta04,1,1,(nely+1)*(nelx+1),1);

topth = [1:nely+1:(nely+1)*nelx+1];
bottomth =[ nely+1:nely+1:(nely+1)*(nelx+1)];
leftth = [1:nely+1];
rightth = [(nely+1)*nelx+1:(nely+1)*(nelx+1)];
T1= zeros((nely+1)*(nelx+1),1);
%% 为边界条件而改造热传导矩阵%
      for ii=1:length(leftth)
       dof=leftth(ii);
          Kth(dof,:) =0;
          Kth(:,dof) =0;
          Kth(dof,dof) =1;
      end
%     
    for ii=1:length(topth)
       dof=topth(ii);
          Kth(dof,:) =0;
          Kth(:,dof) =0;
          Kth(dof,dof) =1;
     end
     for ii=1:length(bottomth)
       dof=bottomth(ii);
          Kth(dof,:) =0;
          Kth(:,dof) =0;
          Kth(dof,dof) =1;
     end
% %  % %  
      for ii=1:length(rightth)
     dof=rightth(ii);
        Kth(dof,:) =0;
        Kth(:,dof) =0;
        Kth(dof,dof) =1;
      end
% %  
% % % %% 添加温度边界条件
    thleftth=0;thbottomth=5;thtopth=5;thrightth=0;
      T1(leftth) =thleftth;
      T1(bottomth)=thbottomth;
      T1(topth) =thtopth;
      T1(rightth)=thrightth;
      Kth2=Kth1-Kth;
     Fth1=Kth2*T1;
     Fth2=Fth-Fth1;
     Fth21=Fth01-Fth1;
     Fth22=Fth02-Fth1;
     Fth2(leftth)=thleftth;
     Fth2(rightth)=thrightth;
     Fth2(bottomth)=thbottomth;
     Fth2(topth)=thtopth;
   %% %%%%%%%%%%%%%%%%%%
     Fth21(leftth)=thleftth;
     Fth21(rightth)=thrightth;
     Fth21(bottomth)=thbottomth;
     Fth21(topth)=thtopth;
     Fth22(leftth)=thleftth;
     Fth22(rightth)=thrightth;
     Fth22(bottomth)=thbottomth;
     Fth22(topth)=thtopth;
  %% 求解温度场
     T = Kth\Fth2;
%       Ti1 =-(Kth\Fth21);
%       Ti2 =-(Kth\Fth22);
       Ti1 =-(Kth\Fth01);
       Ti2 =-(Kth\Fth02);
     Kth0=Kth;
     