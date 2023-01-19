% clc
% clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 几何/物理参数 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[edofMat,iK,jK]=Pre_FEA(nelx,nely);
nelt = nelx*nely;      % Total number of elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  优化参数 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
penal=3.0;rmin=3.0;
x = ones(nelt,1); 

%%
m=1;
n=nelt;
 xold1 = reshape(x,n,1);
 xold2 = reshape(x,n,1);
 xmin = 1e-3; % Densities' Lower bound
 xmax = 1; % Densities' Upper bound
 low = xmin;
 upp = xmax;
a0=1;
a=0;
cmma=10000;
d=0;
move=0.2;
%% Setup Filter
% iW = ones(nelt*(2*(ceil(rmin)-1)+1)^2,1);
% jW = ones(size(iW));
% sW = zeros(size(iW));
% k = 0;
% for i1 = 1:nell
%   for j1 = 1:nell
%     if (i1 > nels && j1 <= nelr); continue; end;
%     e1 = (i1<=nels)*((i1-1)*nell+j1)+(i1>nels)*(nels*nell+(i1-nels-1)*nels+j1-nelr);
%     for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nell)
%       for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nell)
%         if (i2 > nels && j2 <= nelr); continue; end;
%         e2 = (i2<=nels)*((i2-1)*nell+j2)+(i2>nels)*(nels*nell+(i2-nels-1)*nels+j2-nelr);
%         k = k+1;
%         iW(k) = e1;
%         jW(k) = e2;
%         sW(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
%       end
%     end
%   end
% end
% w = sparse(iW,jW,sW);
% W = bsxfun(@rdivide,w,sum(w,2));