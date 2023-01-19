function kk=elementstiff(nelx,nely,length,width,E,u,h)
% E = 2*1e8; 
% u = 0.3;
% h=0.05;
a=length/nelx/2;
b=width/nely/2;
c1=1/3*b/a+(1-u)/6*a/b;
c2=1/3*a/b+(1-u)/6*b/a;
d1=(1+u)/8;
d2=(1-3*u)/8;
e1=-1/3*b/a+(1-u)/12*a/b;
e2=-1/3*a/b+(1-u)/12*a/b;
f1=-1/6*b/a-(1-u)/12*a/b;
f2=-1/6*a/b-(1-u)/12*a/b;
g1=1/6*b/a-(1-u)/6*a/b;
g2=1/6*a/b-(1-u)/6*a/b;
kk=[c1,d1,e1,-d2,f1,-d1,g1,d2;
        d1,c2,d2,g2,-d1,f2,-d2,e2;
        e1,d2,c1,-d1,g1,-d2,f1,d1;
        -d2,g2,-d1,c2,d2,e2,d1,f2;
        f1,-d1,g1,d2,c1,d1,e1,-d2;
        -d1,f2,-d2,e2,d1,c2,d2,g2;
        g1,-d2,f1,d1,e1,d2,c1,-d1;
        d2,e2,d1,f2,-d2,g2,-d1,c2]*E*h/(1-u^2);


