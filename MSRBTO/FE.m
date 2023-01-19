function [U,K,F,DKDE,DKDNU,DKDA,DKDB,lamda]=FE(nelx,nely,x,penal,length,width,P,E,u,h,iK,jK,edofMat)
[KE] = elementstiff(nelx,nely,length,width,E,u,h); 
[dkde,dkdnu,dkda,dkdb]=DEelementstiff(nelx,nely,length,width,E,u,h);%%%单刚对于弹性模量的灵敏度，单刚对于泊松比的灵敏度，单刚对于长度的灵敏度，单刚对于宽度的灵敏度，
nelt = nelx*nely;      % Total number of elements
%%%%%刚度组装
xx=reshape(x,nely,nelx);
sK=reshape(KE(:)*xx(:)'.^penal,64*nelx*nely,1);
K=sparse(iK,jK,sK);K=(K+K')/2;
%%%%%%%%%
sKE=reshape(dkde(:)*xx(:)'.^penal,64*nelx*nely,1);
DKDE = sparse(iK,jK,sKE); DKDE = (DKDE+DKDE')/2; 
sKnu = reshape(dkdnu(:)*xx(:)'.^penal,64*nelx*nely,1);
DKDNU = sparse(iK,jK,sKnu); DKDNU = (DKDNU+DKDNU')/2;
sKa = reshape(dkda(:)*xx(:)'.^penal,64*nelx*nely,1);
DKDA = sparse(iK,jK,sKa); DKDA = (DKDA+DKDA')/2;
sKb = reshape(dkdb(:)*xx(:)'.^penal,64*nelx*nely,1);
DKDB = sparse(iK,jK,sKb); DKDB = (DKDB+DKDB')/2;
%%%%%%%悬臂梁
F(2*(nelx+1)*(nely+1),1) = -P;
fixeddofs = [1:2*(nely+1)];
alldofs     = [1:2*(nely+1)*(nelx+1)];
freedofs    = setdiff(alldofs,fixeddofs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%伴随向量
%%%%伴随向量
T=0*F;
T(2*(nelx+1)*(nely+1),1)=1;
lamda(freedofs,:)=-K(freedofs,freedofs) \ T(freedofs,:);
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);     
U(fixeddofs,:)= 0;

end

