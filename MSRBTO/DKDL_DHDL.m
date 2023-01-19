function [DKDA,DKDB,dHda,dHdb,DKthDA,DKthDB,DKDE,dHdE]=DKDL_DHDL(hou,xPhys_macro,nelx_macro,nely_macro,penal_macro,DH0,DH0_T,E,Emin,iK,jK,iK_T,jK_T,leng,width)
la=leng/nelx_macro/2;lb=width/nely_macro/2;
syms x1 y1 a b 
N(1)=1/4*(1+x1/a)*(1+y1/b);
N(2)=1/4*(1-x1/a)*(1+y1/b);
N(3)=1/4*(1-x1/a)*(1-y1/b);
N(4)=1/4*(1+x1/a)*(1-y1/b);
a11=diff(N(1),x1,1);
a12=diff(N(2),x1,1);
a13=diff(N(3),x1,1);
a14=diff(N(4),x1,1);
a21=diff(N(1),y1,1);
a22=diff(N(2),y1,1);
a23=diff(N(3),y1,1);
a24=diff(N(4),y1,1);
Na=[N(3) N(4) N(1) N(2)];
  B=[a13 0 a14 0 a11 0 a12 0
      0 a23 0 a24 0 a21 0 a22
      a23 a13 a24 a14 a21 a11 a22 a12];
 Bth=[a13 a14 a11 a12 
      a23 a24 a21 a22]; %�ȴ����ļ��ξ���Bth


Q1=B'*DH0*B;
Q2=int(Q1,x1,-a,a);
Ke=int(Q2,y1,-b,b); %����
dKda=double(subs(diff(hou*Ke,a,1),[a b],[la,lb]));  %���նԳ���a�ĵ���
dKdb=double(subs(diff(hou*Ke,b,1),[a b],[la,lb]));  %���նԿ��b�ĵ���

bb=12*10^(-6); %������ϵ��

alpha=[bb bb 0];   %���Ӻ���plani4eҲ����A
H1=B'*DH0*alpha'*Na;
H2=int(H1,x1,-a,a);
H3=int(H2,y1,-b,b); %Ft=x^p*(H3*T)
% Ht=double(subs(H3,[a b],[la,lb]));
dHda=double(subs(diff(hou*H3,a,1),[a b],[la,lb])); %H3�Գ���a�ĵ���
dHdb=double(subs(diff(hou*H3,b,1),[a b],[la,lb])); %H3�Կ��b�ĵ���


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HtDX=cell(nely,nelx);
%   for id = 1:nelx*nely
%          DH1DX=B'*DH0DX{id}*alpha'*Na;
%          DH2DX=int(DH1DX,x1,-a,a);
%          DH3DX=int(DH2DX,y1,-b,b); 
%          DH4DX=double(subs(DH3DX,[a b],[la,lb])); 
%          HtDX{id}=DH4DX;   %Ht��΢����Ʊ�����������
%   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
Qth1=Bth'*DH0_T*Bth;
Qth2=int(Qth1,x1,-a,a);
Kth=int(Qth2,y1,-b,b); %�ȴ����ĵ���
dKthda=double(subs(diff(hou*Kth,a,1),[a b],[la,lb])); %�ȴ����ĵ��նԳ���a�ĵ���
dKthdb=double(subs(diff(hou*Kth,b,1),[a b],[la,lb]));  %�ȴ����ĵ��նԿ��b�ĵ���

dDH0dE=DH0./E;
QE1=B'*dDH0dE*B;
QE2=int(QE1,x1,-a,a);
dKdE=int(QE2,y1,-b,b); %���նԵ���ģ��E�ĵ���
dKdE=double(subs(hou*dKdE,[a b],[la,lb]));

HE1=B'*dDH0dE*alpha'*Na;
HE2=int(HE1,x1,-a,a);
dHdE=int(HE2,y1,-b,b);%dFtdE=x^p*(dHdE*T)
dHdE=double(subs(hou*dHdE,[a b],[la,lb]));

sKa=reshape(dKda(:)*(Emin+xPhys_macro(:)'.^penal_macro*(1-Emin)),64*nelx_macro*nely_macro,1);
DKDA = sparse(iK,jK,sKa);
DKDA = (DKDA+DKDA')/2; 

sKb=reshape(dKdb(:)*(Emin+xPhys_macro(:)'.^penal_macro*(1-Emin)),64*nelx_macro*nely_macro,1);
DKDB = sparse(iK,jK,sKb); DKDB = (DKDB+DKDB')/2; 

sKe=reshape(dKdE(:)*(Emin+xPhys_macro(:)'.^penal_macro*(1-Emin)),64*nelx_macro*nely_macro,1);
DKDE = sparse(iK,jK,sKe); DKDE = (DKDE+DKDE')/2; 

sKtha=reshape(dKthda(:)*(Emin+xPhys_macro(:)'.^penal_macro*(1-Emin)),16*nelx_macro*nely_macro,1);
DKthDA = sparse(iK_T,jK_T,sKtha); DKthDA = (DKthDA+DKthDA')/2; 

sKthb=reshape(dKthdb(:)*(Emin+xPhys_macro(:)'.^penal_macro*(1-Emin)),16*nelx_macro*nely_macro,1);
DKthDB = sparse(iK_T,jK_T,sKthb); DKthDB = (DKthDB+DKthDB')/2; 








