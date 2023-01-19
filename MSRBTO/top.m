function [c,Tmax01,Tmax02,Temfield,f0val,df0dx,df0dx2,fval,dfdx,dgdx] = top(elecount_macro,nodecount_macro,elecount,nodecount,lx,ly,x,xPhys,xPhys_macro,xran,nelx_macro,nely_macro,penal_macro,nely,nelx,penal,rmin_macro,nu,hou,iK,jK,edofMat,iK_T,jK_T,edofMat_T,ft,xTilde,xTilde_macro,eta,beta,H,Hs,H_macro,Hs_macro,volfrac,volfrac_macro,m)
  leng=xran(1);   %长度
  width=xran(2);  %宽度
  E=xran(3);      %弹性模量
  P1(:,1)=xran(4); %载荷
  P2(:,1)=xran(5);
  P3(:,1)=xran(6);
  %% 材料参数
ptype=1;Emin=1e-9;
thickness=1;
integrationRule=2;
ep=[ptype thickness integrationRule];
materialNum=1;
DE=cell(materialNum,2);  %创建数组储存位置，DE=[] [];
% 材料1:各向同性材料
v1 = 0.3;density1=10;  %%材料密度
[DD1]=E/(1-v1^2)*[1  v1   0;
                  v1  1   0;
                  0  0 (1-v1)/2];  
DE{1,1}=DD1;
DE{1,2}=density1; 
[DH0 , DH0DX , H_order0, DH0_T, DH0DX_T,  H_order0_T] = homogenize_order0(lx,ly,xPhys,DE,ep,penal); 

 dx1 = leng/nelx_macro;  
 dy1 = width/nely_macro;
 ex=[0 dx1 dx1 0];
 ey=[0 0 dy1 dy1];

 [KE,~,KE_T,~,H1]=plani4e(ex,ey,ep,DH0,DH0_T); %H1为等效热荷载Ft=H1*T中的等效矩阵系数
  sK = reshape(KE(:)*(Emin+xPhys_macro(:)'.^penal_macro*(1-Emin)),64*nelx_macro*nely_macro,1);
  K = sparse(iK,jK,sK); K = (K+K')/2; 
% 热传导的单刚
    sKT = reshape(KE_T(:)*(Emin+xPhys_macro(:)'.^penal_macro*(1-Emin)),16*nelx_macro*nely_macro,1);
    Kth = sparse(iK_T,jK_T,sKT); Kth = (Kth+Kth')/2;
    Kth1=Kth;
  %% 求解温度场 temfield
 [Ta01,Ta04,Ti1,Ti2,T,iAA]=tem1_field(Kth,nelx_macro,nely_macro);
  Temfield=reshape(T,(nely_macro+1),(nelx_macro+1));
  
  %% 对各种随机变量的偏导   
  [DKDA,DKDB,dHda,dHdb,DKthDA,DKthDB,DKDE,dHde]=DKDL_DHDL(hou,xPhys_macro,nelx_macro,nely_macro,penal_macro,DH0,DH0_T,E,Emin,iK,jK,iK_T,jK_T,leng,width);
  
  %% 求解等效温度荷载
   Ht=H1;
   HH = zeros(2*(nely_macro+1)*(nelx_macro+1),(nely_macro+1)*(nelx_macro+1));
   DHHDX = zeros(2*(nely_macro+1)*(nelx_macro+1),(nely_macro+1)*(nelx_macro+1));
   DHDA = zeros(2*(nely_macro+1)*(nelx_macro+1),(nely_macro+1)*(nelx_macro+1));
   DHDB = zeros(2*(nely_macro+1)*(nelx_macro+1),(nely_macro+1)*(nelx_macro+1));
   DHDE = zeros(2*(nely_macro+1)*(nelx_macro+1),(nely_macro+1)*(nelx_macro+1));
 for elx = 1:nelx_macro
     for ely = 1:nely_macro
       n1 = (nely_macro+1)*(elx-1)+ely; 
       n2 = (nely_macro+1)* elx   +ely;
       edof4 = [n1+1; n2+1; n2; n1];
       edof8 = [2*n1+1; 2*n1+2; 2*n2+1; 2*n2+2; 2*n2-1; 2*n2; 2*n1-1; 2*n1];
       HH(edof8,edof4) = HH(edof8,edof4) + (Emin+(1-Emin)*xPhys_macro(ely,elx).^penal_macro).*Ht;%%%
       DHHDX(edof8,edof4) = DHHDX(edof8,edof4) + penal_macro*xPhys_macro(ely,elx)^(penal_macro-1)*Ht;
       DHDA(edof8,edof4) = DHDA(edof8,edof4) + (Emin+xPhys_macro(ely,elx).^penal_macro*(1-Emin))*dHda;
       DHDB(edof8,edof4) = DHDB(edof8,edof4) + (Emin+xPhys_macro(ely,elx).^penal_macro*(1-Emin))*dHdb;
       DHDE(edof8,edof4) = DHDE(edof8,edof4) + (Emin+xPhys_macro(ely,elx).^penal_macro*(1-Emin))*dHde;
     end
  end 
  Ft=HH*T;
  
  %% 荷载和约束条件
  iF = [26 (nely_macro+1)*nelx_macro/2+1 (nely_macro+1)*nelx_macro+26];
  Fm1 = sparse(2*26 ,1,-P1,2*(nely_macro+1)*(nelx_macro+1),1);
  Fm2 = sparse(2*((nely_macro+1)*nelx_macro/2+1) ,1,-P2,2*(nely_macro+1)*(nelx_macro+1),1);
  Fm3 = sparse(2*((nely_macro+1)*nelx_macro+26),1,-P3,2*(nely_macro+1)*(nelx_macro+1),1);
  U = zeros(2*(nely_macro+1)*(nelx_macro+1),1); 
  BOTTOM=[ (nely_macro+1)*nelx_macro/4:nely_macro+1:(nely_macro+1)*3*nelx_macro/4];
  fixeddofs =[2*BOTTOM 2*BOTTOM-1];
  alldofs = [1:2*(nely_macro+1)*(nelx_macro+1)];
  freedofs = setdiff(alldofs,fixeddofs);
 Fm=Fm1+Fm2+Fm3;
  %耦合力
F=Fm+Ft;  

U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    % 柔度
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely_macro,nelx_macro); %柔度
    c = sum(sum((Emin+xPhys_macro.^penal_macro*(1-Emin)).*ce));  

    lamda_T=0*T;
    GG=(-2*U'*HH)';
    lamda_T(iAA)=Kth(iAA,iAA)\GG(iAA); % C=U'*K*U+lamda*(K*U-Fm-Ft)+lamda_T*(Kth*T-Q)

   %% 宏观敏度(下面写法速度较快)
%    T1_macro=reshape(sum((Ti1(edofMat_T)*KE_T).*T(edofMat_T),2),nely_macro,nelx_macro); 
%    dTdX1_macro=penal_macro*(1-Emin)*(Emin+xPhys_macro(ely,elx)).^(penal_macro-1).*T1_macro;
%    T2_macro=reshape(sum((Ti2(edofMat_T)*KE_T).*T(edofMat_T),2),nely_macro,nelx_macro); 
%    dTdX2_macro=penal_macro*(1-Emin)*(Emin+xPhys_macro(ely,elx)).^(penal_macro-1).*T2_macro;
%    C1_macro=reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely_macro,nelx_macro); 
%    C2_macro=reshape(sum((U(edofMat)*Ht).*T(edofMat_T),2),nely_macro,nelx_macro); 
%    C3_macro=reshape(sum((lamda_T(edofMat_T)*KE_T).*T(edofMat_T),2),nely_macro,nelx_macro);
%    dTdX2_macro=penal_macro*(1-Emin)*(Emin+xPhys_macro(ely,elx)).^(penal_macro-1).*T2_macro +...
%    2*penal_macro*(1-Emin)*(Emin+xPhys_macro(ely,elx)).^(penal_macro-1)*C2_macro + penal_macro*(1-Emin)*(Emin+xPhys_macro(ely,elx)).^(penal_macro-1)*C3_macro;
  
%% 宏观敏度
   for elx = 1:nelx_macro
 for ely = 1:nely_macro
     n1 = (nely_macro+1)*(elx-1)+ely;
     n2 = (nely_macro+1)* elx +ely;
     edof4 = [n1+1; n2+1; n2; n1];
     edof8 = [2*n1+1; 2*n1+2; 2*n2+1; 2*n2+2; 2*n2-1; 2*n2; 2*n1-1; 2*n1];
     Ue = U(edof8,1);
     Te = T(edof4,1);
     Ti1e=Ti1(edof4,1);
     Ti2e=Ti2(edof4,1);
     dT1_macro(ely,elx)=penal_macro*(1-Emin)*(Emin+xPhys_macro(ely,elx)).^(penal_macro-1).*(Ti1e'*KE_T*Te);
     dT2_macro(ely,elx)=penal_macro*(1-Emin)*(Emin+xPhys_macro(ely,elx)).^(penal_macro-1).*(Ti2e'*KE_T*Te);
     DHHDXe=DHHDX(edof8,edof4);
     lamda_Te = lamda_T(edof4,1);
     dc_macro(ely,elx) = -penal_macro*(1-Emin)*(Emin+xPhys_macro(ely,elx)).^(penal_macro-1)*(Ue'*KE*Ue)+2*penal_macro*(1-Emin)*(Emin+xPhys_macro(ely,elx)).^(penal_macro-1)*(Ue'*Ht*Te)+ ...
     penal_macro*(1-Emin)*(Emin+xPhys_macro(ely,elx)).^(penal_macro-1)*(lamda_Te'*KE_T*Te);
 end
end
     dv_macro = ones(nely_macro,nelx_macro);
  %% 微观敏度
    dc =0*xPhys;  %开辟储存空间
    dT1=0*xPhys;
    dT2=0*xPhys;

    for i=1:elecount
        dQ=DH0DX{i};
        dQ_T=DH0DX_T{i};
        [KEDX,~,KE_TDX,~,H1DX]=plani4e(ex,ey,ep,dQ,dQ_T);   %有问题
        sKTDX = reshape(KE_TDX(:)*(Emin+xPhys_macro(:)'.^penal_macro*(1-Emin)),16*nelx_macro*nely_macro,1);
        KthDX_mic = sparse(iK_T,jK_T,sKTDX); KthDX_mic = (KthDX_mic+KthDX_mic')/2;
        ce1 = reshape(sum((U(edofMat)*KEDX).*U(edofMat),2),nely_macro,nelx_macro); 
        ce2 = reshape(sum((2*U(edofMat)*H1DX).*T(edofMat_T),2),nely_macro,nelx_macro);
        ce3 = reshape(sum((lamda_T(edofMat_T)*hou*KE_TDX).*T(edofMat_T),2),nely_macro,nelx_macro);
        dc(i) = -sum(sum((Emin+xPhys_macro.^penal_macro*(1-Emin)).*ce1))+... 
            sum(sum((Emin+xPhys_macro.^penal_macro*(1-Emin)).*ce2))+sum(sum((Emin+xPhys_macro.^penal_macro*(1-Emin)).*ce3)); 
        dT1(i)=Ti1'*KthDX_mic*T;
        dT2(i)=Ti2'*KthDX_mic*T; %微观节点温度敏度
    end
    dv = ones(nely,nelx);
 %% FILTERING/MODIFICATION OF SENSITIVITIES过滤/修改敏感
    if ft == 1
        dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
        dT1(:) = H*(x(:).*dT1(:))./Hs./max(1e-3,x(:));
        dT2(:) = H*(x(:).*dT2(:))./Hs./max(1e-3,x(:));
        dc_macro(:) = H_macro*(x_macro(:).*dc_macro(:))./Hs_macro./max(1e-3,x_macro(:));
        dT1_macro(:) = H_macro*(x_macro(:).*dT1_macro(:))./Hs_macro./max(1e-3,x_macro(:));
        dT2_macro(:) = H_macro*(x_macro(:).*dT2_macro(:))./Hs_macro./max(1e-3,x_macro(:));
    elseif ft == 2
        dc(:) = H*(dc(:)./Hs);
        dT1(:) = H*(dT1(:)./Hs);
        dT2(:) = H*(dT2(:)./Hs);
        dv(:) = H*(dv(:)./Hs);
        dc_macro(:) = H_macro*(dc_macro(:)./Hs_macro);
        dT1_macro(:) = H_macro*(dT1_macro(:)./Hs_macro);
        dT2_macro(:) = H_macro*(dT2_macro(:)./Hs_macro);
        dv_macro(:) = H_macro*(dv_macro(:)./Hs_macro);
    elseif ft == 3
        dx =HeavisideSns(xTilde,beta,eta);
        dc(:) = H*(dc(:).*dx(:)./Hs);
        dT1(:) = H*(dT1(:).*dx(:)./Hs);
        dT2(:) = H*(dT2(:).*dx(:)./Hs);
        dv(:) = H*(dv(:).*dx(:)./Hs);
        dx_macro =HeavisideSns(xTilde_macro,beta,eta);
        dc_macro(:) = H_macro*(dc_macro(:).*dx_macro(:)./Hs_macro);
        dT1_macro(:) = H_macro*(dT1_macro(:).*dx_macro(:)./Hs_macro);
        dT2_macro(:) = H_macro*(dT2_macro(:).*dx_macro(:)./Hs_macro);
        dv_macro(:) = H_macro*(dv_macro(:).*dx_macro(:)./Hs_macro);
    end  
   n=elecount_macro+elecount;          %变量个数1
   
%% %%%%%%%%%%%%%%约束条件+目标函数%%%%%%%%%%%%%%%%%
      W1=1e0;M=1e0;W2=1e0;W3=1e0;
     v=sum(xPhys(:))-volfrac*nely*nelx;%%
     v_macro=sum(xPhys_macro(:));
      f0val=M*v_macro;
      df0dx =M*[0*dv(:);dv_macro(:)];
      df0dx2 =zeros(n,1);
      obc=c-17;
   
    Tmax=10; %热源处温度允许值
    Tmax01=T(Ta01);
    Tmax02=T(Ta04);
    Tmax1=(Tmax01-Tmax); 
    Tmax2=(Tmax02-Tmax); 
      fval =[W3*v;W1*obc;W2*Tmax1;W2*Tmax2];
      dfdx=zeros(m,n);
      dfdx(1,1:elecount)=W3*dv(:);
      dfdx(2,:)=W1*[dc(:) ; dc_macro(:)];
      dfdx(3,:)=W2*[dT1(:);dT1_macro(:)]; %先微观后宏观
      dfdx(4,:)=W2*[dT2(:);dT2_macro(:)];

%% 微观体积对对宏观随机变量的灵敏度
   dgdx(1,1)=0;
   dgdx(1,2)=0;  
   dgdx(1,3)=0;  
   dgdx(1,4)=0; 
   dgdx(1,5)=0; 
   dgdx(1,6)=0;  
  %% 柔度对宏观随机变量的灵敏度
   dgdx(2,1)=-U'*DKDA*U+2*U'*DHDA*T+lamda_T'*DKthDA*T;
   dgdx(2,2)=-U'*DKDB*U+2*U'*DHDB*T+lamda_T'*DKthDB*T;
   dgdx(2,3)=-U'*DKDE*U+2*U'*DHDE*T;
   dgdx(2,4)=2*U'*(Fm1/P1);
   dgdx(2,5)=2*U'*(Fm2/P2);
   dgdx(2,6)=2*U'*(Fm3/P3);
 %% 两节点温度对宏观随机变量的灵敏度
  dgdx(3,1)=Ti1'*DKthDA*T;
  dgdx(3,2)=Ti1'*DKthDB*T;
  dgdx(3,3)=0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
  dgdx(3,4)=0;
  dgdx(3,5)=0;
  dgdx(3,6)=0;
  dgdx(4,1)=Ti2'*DKthDA*T;
  dgdx(4,2)=Ti2'*DKthDB*T;
  dgdx(4,3)=0;
  dgdx(4,4)=0;
  dgdx(4,5)=0;
  dgdx(4,6)=0;
   fval=-fval;
   dgdx=-dgdx;
   dfdx=-dfdx;
