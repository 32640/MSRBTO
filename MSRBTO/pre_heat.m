function [iK,jK,edofMat]=pre_heat(nelx,nely)
%%¡ª¡ª¡ª¡ª¡ª×óÏÂ½Ç¿ªÊ¼±àÂë£º
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,4)+repmat([0 nely+[1 0] -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(4,1))',16*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,4))',16*nelx*nely,1);