function [mu,sigma,type]=distribution
% mu=[16]'; %���ȣ���ȣ�����ģ�����غ�
% type=[1];
% sigma=[0.1]'.*mu;
F=2000;
mu=[0.1;0.13;2e11;F;2*F;F]'; %���ȣ���ȣ�����ģ�����غ�
type=[1;1;1;1;1;1];
sigma=0.1*[0.1;0.13;2e11;F;2*F;F]';  % S G M