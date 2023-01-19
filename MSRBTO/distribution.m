function [mu,sigma,type]=distribution
% mu=[16]'; %长度，宽度，弹性模量，载荷
% type=[1];
% sigma=[0.1]'.*mu;
F=2000;
mu=[0.1;0.13;2e11;F;2*F;F]'; %长度，宽度，弹性模量，载荷
type=[1;1;1;1;1;1];
sigma=0.1*[0.1;0.13;2e11;F;2*F;F]';  % S G M