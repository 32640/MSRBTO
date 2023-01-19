function [dkde,dkdnu,dkda,dkdb]=DEelementstiff(nelx,nely,leng,width,E,NU,hou)
% E = 2*1e8; 
% NU = 0.3;
h=hou;
a=leng/nelx/2;
b=width/nely/2;
%%%%%%%%%%%%刚度对于弹性模量E的矩阵
dkde=[ -(h*(a^2 - NU*a^2 + 2*b^2))/(6*a*b*(NU^2 - 1)),                                -h/(8*(NU - 1)), (h*(NU*a^2 - a^2 + 4*b^2))/(12*a*b*(NU^2 - 1)),                 -(h*(3*NU - 1))/(8*(NU^2 - 1)), (h*(a^2 - NU*a^2 + 2*b^2))/(12*a*b*(NU^2 - 1)),                                 h/(8*(NU - 1)),   -(h*(NU*a^2 - a^2 + b^2))/(6*a*b*(NU^2 - 1)),                  (h*(3*NU - 1))/(8*(NU^2 - 1));
                                -h/(8*(NU - 1)), -(h*(2*a^2 - NU*b^2 + b^2))/(6*a*b*(NU^2 - 1)),                  (h*(3*NU - 1))/(8*(NU^2 - 1)),   -(h*(NU*b^2 + a^2 - b^2))/(6*a*b*(NU^2 - 1)),                                 h/(8*(NU - 1)), (h*(2*a^2 - NU*b^2 + b^2))/(12*a*b*(NU^2 - 1)),                 -(h*(3*NU - 1))/(8*(NU^2 - 1)), (h*(NU*b^2 + 4*a^2 - b^2))/(12*a*b*(NU^2 - 1));
 (h*(NU*a^2 - a^2 + 4*b^2))/(12*a*b*(NU^2 - 1)),                  (h*(3*NU - 1))/(8*(NU^2 - 1)), -(h*(a^2 - NU*a^2 + 2*b^2))/(6*a*b*(NU^2 - 1)),                                 h/(8*(NU - 1)),   -(h*(NU*a^2 - a^2 + b^2))/(6*a*b*(NU^2 - 1)),                 -(h*(3*NU - 1))/(8*(NU^2 - 1)), (h*(a^2 - NU*a^2 + 2*b^2))/(12*a*b*(NU^2 - 1)),                                -h/(8*(NU - 1));
                 -(h*(3*NU - 1))/(8*(NU^2 - 1)),   -(h*(NU*b^2 + a^2 - b^2))/(6*a*b*(NU^2 - 1)),                                 h/(8*(NU - 1)), -(h*(2*a^2 - NU*b^2 + b^2))/(6*a*b*(NU^2 - 1)),                  (h*(3*NU - 1))/(8*(NU^2 - 1)), (h*(NU*b^2 + 4*a^2 - b^2))/(12*a*b*(NU^2 - 1)),                                -h/(8*(NU - 1)), (h*(2*a^2 - NU*b^2 + b^2))/(12*a*b*(NU^2 - 1));
 (h*(a^2 - NU*a^2 + 2*b^2))/(12*a*b*(NU^2 - 1)),                                 h/(8*(NU - 1)),   -(h*(NU*a^2 - a^2 + b^2))/(6*a*b*(NU^2 - 1)),                  (h*(3*NU - 1))/(8*(NU^2 - 1)), -(h*(a^2 - NU*a^2 + 2*b^2))/(6*a*b*(NU^2 - 1)),                                -h/(8*(NU - 1)), (h*(NU*a^2 - a^2 + 4*b^2))/(12*a*b*(NU^2 - 1)),                 -(h*(3*NU - 1))/(8*(NU^2 - 1));
                                 h/(8*(NU - 1)), (h*(2*a^2 - NU*b^2 + b^2))/(12*a*b*(NU^2 - 1)),                 -(h*(3*NU - 1))/(8*(NU^2 - 1)), (h*(NU*b^2 + 4*a^2 - b^2))/(12*a*b*(NU^2 - 1)),                                -h/(8*(NU - 1)), -(h*(2*a^2 - NU*b^2 + b^2))/(6*a*b*(NU^2 - 1)),                  (h*(3*NU - 1))/(8*(NU^2 - 1)),   -(h*(NU*b^2 + a^2 - b^2))/(6*a*b*(NU^2 - 1));
   -(h*(NU*a^2 - a^2 + b^2))/(6*a*b*(NU^2 - 1)),                 -(h*(3*NU - 1))/(8*(NU^2 - 1)), (h*(a^2 - NU*a^2 + 2*b^2))/(12*a*b*(NU^2 - 1)),                                -h/(8*(NU - 1)), (h*(NU*a^2 - a^2 + 4*b^2))/(12*a*b*(NU^2 - 1)),                  (h*(3*NU - 1))/(8*(NU^2 - 1)), -(h*(a^2 - NU*a^2 + 2*b^2))/(6*a*b*(NU^2 - 1)),                                 h/(8*(NU - 1));
                  (h*(3*NU - 1))/(8*(NU^2 - 1)), (h*(NU*b^2 + 4*a^2 - b^2))/(12*a*b*(NU^2 - 1)),                                -h/(8*(NU - 1)), (h*(2*a^2 - NU*b^2 + b^2))/(12*a*b*(NU^2 - 1)),                 -(h*(3*NU - 1))/(8*(NU^2 - 1)),   -(h*(NU*b^2 + a^2 - b^2))/(6*a*b*(NU^2 - 1)),                                 h/(8*(NU - 1)), -(h*(2*a^2 - NU*b^2 + b^2))/(6*a*b*(NU^2 - 1))];
%%%%%%%%%%%%刚度对于泊松比NU的矩阵
dkdnu=[    (E*a*h)/(6*b*(NU^2 - 1)) + (E*NU*h*(a^2 - NU*a^2 + 2*b^2))/(3*a*b*(NU^2 - 1)^2),                                                               (E*h)/(8*(NU - 1)^2),   (E*a*h)/(12*b*(NU^2 - 1)) - (E*NU*h*(NU*a^2 - a^2 + 4*b^2))/(6*a*b*(NU^2 - 1)^2),                      (E*NU*h*(3*NU - 1))/(4*(NU^2 - 1)^2) - (3*E*h)/(8*(NU^2 - 1)), - (E*a*h)/(12*b*(NU^2 - 1)) - (E*NU*h*(a^2 - NU*a^2 + 2*b^2))/(6*a*b*(NU^2 - 1)^2),                                                              -(E*h)/(8*(NU - 1)^2),      (E*NU*h*(NU*a^2 - a^2 + b^2))/(3*a*b*(NU^2 - 1)^2) - (E*a*h)/(6*b*(NU^2 - 1)),                      (3*E*h)/(8*(NU^2 - 1)) - (E*NU*h*(3*NU - 1))/(4*(NU^2 - 1)^2);
                                                               (E*h)/(8*(NU - 1)^2),    (E*b*h)/(6*a*(NU^2 - 1)) + (E*NU*h*(2*a^2 - NU*b^2 + b^2))/(3*a*b*(NU^2 - 1)^2),                      (3*E*h)/(8*(NU^2 - 1)) - (E*NU*h*(3*NU - 1))/(4*(NU^2 - 1)^2),      (E*NU*h*(NU*b^2 + a^2 - b^2))/(3*a*b*(NU^2 - 1)^2) - (E*b*h)/(6*a*(NU^2 - 1)),                                                              -(E*h)/(8*(NU - 1)^2), - (E*b*h)/(12*a*(NU^2 - 1)) - (E*NU*h*(2*a^2 - NU*b^2 + b^2))/(6*a*b*(NU^2 - 1)^2),                      (E*NU*h*(3*NU - 1))/(4*(NU^2 - 1)^2) - (3*E*h)/(8*(NU^2 - 1)),   (E*b*h)/(12*a*(NU^2 - 1)) - (E*NU*h*(NU*b^2 + 4*a^2 - b^2))/(6*a*b*(NU^2 - 1)^2);
   (E*a*h)/(12*b*(NU^2 - 1)) - (E*NU*h*(NU*a^2 - a^2 + 4*b^2))/(6*a*b*(NU^2 - 1)^2),                      (3*E*h)/(8*(NU^2 - 1)) - (E*NU*h*(3*NU - 1))/(4*(NU^2 - 1)^2),    (E*a*h)/(6*b*(NU^2 - 1)) + (E*NU*h*(a^2 - NU*a^2 + 2*b^2))/(3*a*b*(NU^2 - 1)^2),                                                              -(E*h)/(8*(NU - 1)^2),      (E*NU*h*(NU*a^2 - a^2 + b^2))/(3*a*b*(NU^2 - 1)^2) - (E*a*h)/(6*b*(NU^2 - 1)),                      (E*NU*h*(3*NU - 1))/(4*(NU^2 - 1)^2) - (3*E*h)/(8*(NU^2 - 1)), - (E*a*h)/(12*b*(NU^2 - 1)) - (E*NU*h*(a^2 - NU*a^2 + 2*b^2))/(6*a*b*(NU^2 - 1)^2),                                                               (E*h)/(8*(NU - 1)^2);
                      (E*NU*h*(3*NU - 1))/(4*(NU^2 - 1)^2) - (3*E*h)/(8*(NU^2 - 1)),      (E*NU*h*(NU*b^2 + a^2 - b^2))/(3*a*b*(NU^2 - 1)^2) - (E*b*h)/(6*a*(NU^2 - 1)),                                                              -(E*h)/(8*(NU - 1)^2),    (E*b*h)/(6*a*(NU^2 - 1)) + (E*NU*h*(2*a^2 - NU*b^2 + b^2))/(3*a*b*(NU^2 - 1)^2),                      (3*E*h)/(8*(NU^2 - 1)) - (E*NU*h*(3*NU - 1))/(4*(NU^2 - 1)^2),   (E*b*h)/(12*a*(NU^2 - 1)) - (E*NU*h*(NU*b^2 + 4*a^2 - b^2))/(6*a*b*(NU^2 - 1)^2),                                                               (E*h)/(8*(NU - 1)^2), - (E*b*h)/(12*a*(NU^2 - 1)) - (E*NU*h*(2*a^2 - NU*b^2 + b^2))/(6*a*b*(NU^2 - 1)^2);
 - (E*a*h)/(12*b*(NU^2 - 1)) - (E*NU*h*(a^2 - NU*a^2 + 2*b^2))/(6*a*b*(NU^2 - 1)^2),                                                              -(E*h)/(8*(NU - 1)^2),      (E*NU*h*(NU*a^2 - a^2 + b^2))/(3*a*b*(NU^2 - 1)^2) - (E*a*h)/(6*b*(NU^2 - 1)),                      (3*E*h)/(8*(NU^2 - 1)) - (E*NU*h*(3*NU - 1))/(4*(NU^2 - 1)^2),    (E*a*h)/(6*b*(NU^2 - 1)) + (E*NU*h*(a^2 - NU*a^2 + 2*b^2))/(3*a*b*(NU^2 - 1)^2),                                                               (E*h)/(8*(NU - 1)^2),   (E*a*h)/(12*b*(NU^2 - 1)) - (E*NU*h*(NU*a^2 - a^2 + 4*b^2))/(6*a*b*(NU^2 - 1)^2),                      (E*NU*h*(3*NU - 1))/(4*(NU^2 - 1)^2) - (3*E*h)/(8*(NU^2 - 1));
                                                              -(E*h)/(8*(NU - 1)^2), - (E*b*h)/(12*a*(NU^2 - 1)) - (E*NU*h*(2*a^2 - NU*b^2 + b^2))/(6*a*b*(NU^2 - 1)^2),                      (E*NU*h*(3*NU - 1))/(4*(NU^2 - 1)^2) - (3*E*h)/(8*(NU^2 - 1)),   (E*b*h)/(12*a*(NU^2 - 1)) - (E*NU*h*(NU*b^2 + 4*a^2 - b^2))/(6*a*b*(NU^2 - 1)^2),                                                               (E*h)/(8*(NU - 1)^2),    (E*b*h)/(6*a*(NU^2 - 1)) + (E*NU*h*(2*a^2 - NU*b^2 + b^2))/(3*a*b*(NU^2 - 1)^2),                      (3*E*h)/(8*(NU^2 - 1)) - (E*NU*h*(3*NU - 1))/(4*(NU^2 - 1)^2),      (E*NU*h*(NU*b^2 + a^2 - b^2))/(3*a*b*(NU^2 - 1)^2) - (E*b*h)/(6*a*(NU^2 - 1));
      (E*NU*h*(NU*a^2 - a^2 + b^2))/(3*a*b*(NU^2 - 1)^2) - (E*a*h)/(6*b*(NU^2 - 1)),                      (E*NU*h*(3*NU - 1))/(4*(NU^2 - 1)^2) - (3*E*h)/(8*(NU^2 - 1)), - (E*a*h)/(12*b*(NU^2 - 1)) - (E*NU*h*(a^2 - NU*a^2 + 2*b^2))/(6*a*b*(NU^2 - 1)^2),                                                               (E*h)/(8*(NU - 1)^2),   (E*a*h)/(12*b*(NU^2 - 1)) - (E*NU*h*(NU*a^2 - a^2 + 4*b^2))/(6*a*b*(NU^2 - 1)^2),                      (3*E*h)/(8*(NU^2 - 1)) - (E*NU*h*(3*NU - 1))/(4*(NU^2 - 1)^2),    (E*a*h)/(6*b*(NU^2 - 1)) + (E*NU*h*(a^2 - NU*a^2 + 2*b^2))/(3*a*b*(NU^2 - 1)^2),                                                              -(E*h)/(8*(NU - 1)^2);
                      (3*E*h)/(8*(NU^2 - 1)) - (E*NU*h*(3*NU - 1))/(4*(NU^2 - 1)^2),   (E*b*h)/(12*a*(NU^2 - 1)) - (E*NU*h*(NU*b^2 + 4*a^2 - b^2))/(6*a*b*(NU^2 - 1)^2),                                                               (E*h)/(8*(NU - 1)^2), - (E*b*h)/(12*a*(NU^2 - 1)) - (E*NU*h*(2*a^2 - NU*b^2 + b^2))/(6*a*b*(NU^2 - 1)^2),                      (E*NU*h*(3*NU - 1))/(4*(NU^2 - 1)^2) - (3*E*h)/(8*(NU^2 - 1)),      (E*NU*h*(NU*b^2 + a^2 - b^2))/(3*a*b*(NU^2 - 1)^2) - (E*b*h)/(6*a*(NU^2 - 1)),                                                              -(E*h)/(8*(NU - 1)^2),    (E*b*h)/(6*a*(NU^2 - 1)) + (E*NU*h*(2*a^2 - NU*b^2 + b^2))/(3*a*b*(NU^2 - 1)^2)];
%%%%%%%%%%%%刚度对于长度a的矩阵
dkda=1/nelx/2*[     (E*h*(a^2 - NU*a^2 + 2*b^2))/(6*a^2*b*(NU^2 - 1)) - (E*h*(2*a - 2*NU*a))/(6*a*b*(NU^2 - 1)),                                                                             0, - (E*h*(NU*a^2 - a^2 + 4*b^2))/(12*a^2*b*(NU^2 - 1)) - (E*h*(2*a - 2*NU*a))/(12*a*b*(NU^2 - 1)),                                                                             0,   (E*h*(2*a - 2*NU*a))/(12*a*b*(NU^2 - 1)) - (E*h*(a^2 - NU*a^2 + 2*b^2))/(12*a^2*b*(NU^2 - 1)),                                                                             0,       (E*h*(2*a - 2*NU*a))/(6*a*b*(NU^2 - 1)) + (E*h*(NU*a^2 - a^2 + b^2))/(6*a^2*b*(NU^2 - 1)),                                                                             0;
                                                                                               0,  (E*h*(2*a^2 - NU*b^2 + b^2))/(6*a^2*b*(NU^2 - 1)) - (2*E*h)/(3*b*(NU^2 - 1)),                                                                                               0,      (E*h*(NU*b^2 + a^2 - b^2))/(6*a^2*b*(NU^2 - 1)) - (E*h)/(3*b*(NU^2 - 1)),                                                                                               0,   (E*h)/(3*b*(NU^2 - 1)) - (E*h*(2*a^2 - NU*b^2 + b^2))/(12*a^2*b*(NU^2 - 1)),                                                                                               0, (2*E*h)/(3*b*(NU^2 - 1)) - (E*h*(NU*b^2 + 4*a^2 - b^2))/(12*a^2*b*(NU^2 - 1));
 - (E*h*(NU*a^2 - a^2 + 4*b^2))/(12*a^2*b*(NU^2 - 1)) - (E*h*(2*a - 2*NU*a))/(12*a*b*(NU^2 - 1)),                                                                             0,     (E*h*(a^2 - NU*a^2 + 2*b^2))/(6*a^2*b*(NU^2 - 1)) - (E*h*(2*a - 2*NU*a))/(6*a*b*(NU^2 - 1)),                                                                             0,       (E*h*(2*a - 2*NU*a))/(6*a*b*(NU^2 - 1)) + (E*h*(NU*a^2 - a^2 + b^2))/(6*a^2*b*(NU^2 - 1)),                                                                             0,   (E*h*(2*a - 2*NU*a))/(12*a*b*(NU^2 - 1)) - (E*h*(a^2 - NU*a^2 + 2*b^2))/(12*a^2*b*(NU^2 - 1)),                                                                             0;
                                                                                               0,      (E*h*(NU*b^2 + a^2 - b^2))/(6*a^2*b*(NU^2 - 1)) - (E*h)/(3*b*(NU^2 - 1)),                                                                                               0,  (E*h*(2*a^2 - NU*b^2 + b^2))/(6*a^2*b*(NU^2 - 1)) - (2*E*h)/(3*b*(NU^2 - 1)),                                                                                               0, (2*E*h)/(3*b*(NU^2 - 1)) - (E*h*(NU*b^2 + 4*a^2 - b^2))/(12*a^2*b*(NU^2 - 1)),                                                                                               0,   (E*h)/(3*b*(NU^2 - 1)) - (E*h*(2*a^2 - NU*b^2 + b^2))/(12*a^2*b*(NU^2 - 1));
   (E*h*(2*a - 2*NU*a))/(12*a*b*(NU^2 - 1)) - (E*h*(a^2 - NU*a^2 + 2*b^2))/(12*a^2*b*(NU^2 - 1)),                                                                             0,       (E*h*(2*a - 2*NU*a))/(6*a*b*(NU^2 - 1)) + (E*h*(NU*a^2 - a^2 + b^2))/(6*a^2*b*(NU^2 - 1)),                                                                             0,     (E*h*(a^2 - NU*a^2 + 2*b^2))/(6*a^2*b*(NU^2 - 1)) - (E*h*(2*a - 2*NU*a))/(6*a*b*(NU^2 - 1)),                                                                             0, - (E*h*(NU*a^2 - a^2 + 4*b^2))/(12*a^2*b*(NU^2 - 1)) - (E*h*(2*a - 2*NU*a))/(12*a*b*(NU^2 - 1)),                                                                             0;
                                                                                               0,   (E*h)/(3*b*(NU^2 - 1)) - (E*h*(2*a^2 - NU*b^2 + b^2))/(12*a^2*b*(NU^2 - 1)),                                                                                               0, (2*E*h)/(3*b*(NU^2 - 1)) - (E*h*(NU*b^2 + 4*a^2 - b^2))/(12*a^2*b*(NU^2 - 1)),                                                                                               0,  (E*h*(2*a^2 - NU*b^2 + b^2))/(6*a^2*b*(NU^2 - 1)) - (2*E*h)/(3*b*(NU^2 - 1)),                                                                                               0,      (E*h*(NU*b^2 + a^2 - b^2))/(6*a^2*b*(NU^2 - 1)) - (E*h)/(3*b*(NU^2 - 1));
       (E*h*(2*a - 2*NU*a))/(6*a*b*(NU^2 - 1)) + (E*h*(NU*a^2 - a^2 + b^2))/(6*a^2*b*(NU^2 - 1)),                                                                             0,   (E*h*(2*a - 2*NU*a))/(12*a*b*(NU^2 - 1)) - (E*h*(a^2 - NU*a^2 + 2*b^2))/(12*a^2*b*(NU^2 - 1)),                                                                             0, - (E*h*(NU*a^2 - a^2 + 4*b^2))/(12*a^2*b*(NU^2 - 1)) - (E*h*(2*a - 2*NU*a))/(12*a*b*(NU^2 - 1)),                                                                             0,     (E*h*(a^2 - NU*a^2 + 2*b^2))/(6*a^2*b*(NU^2 - 1)) - (E*h*(2*a - 2*NU*a))/(6*a*b*(NU^2 - 1)),                                                                             0;
                                                                                               0, (2*E*h)/(3*b*(NU^2 - 1)) - (E*h*(NU*b^2 + 4*a^2 - b^2))/(12*a^2*b*(NU^2 - 1)),                                                                                               0,   (E*h)/(3*b*(NU^2 - 1)) - (E*h*(2*a^2 - NU*b^2 + b^2))/(12*a^2*b*(NU^2 - 1)),                                                                                               0,      (E*h*(NU*b^2 + a^2 - b^2))/(6*a^2*b*(NU^2 - 1)) - (E*h)/(3*b*(NU^2 - 1)),                                                                                               0,  (E*h*(2*a^2 - NU*b^2 + b^2))/(6*a^2*b*(NU^2 - 1)) - (2*E*h)/(3*b*(NU^2 - 1))];
%%%%%%%%%%%%刚度对于宽度b的矩阵    
dkdb=1/nely/2*[  (E*h*(a^2 - NU*a^2 + 2*b^2))/(6*a*b^2*(NU^2 - 1)) - (2*E*h)/(3*a*(NU^2 - 1)),                                                                                               0, (2*E*h)/(3*a*(NU^2 - 1)) - (E*h*(NU*a^2 - a^2 + 4*b^2))/(12*a*b^2*(NU^2 - 1)),                                                                                               0,   (E*h)/(3*a*(NU^2 - 1)) - (E*h*(a^2 - NU*a^2 + 2*b^2))/(12*a*b^2*(NU^2 - 1)),                                                                                               0,      (E*h*(NU*a^2 - a^2 + b^2))/(6*a*b^2*(NU^2 - 1)) - (E*h)/(3*a*(NU^2 - 1)),                                                                                               0;
                                                                             0,     (E*h*(2*a^2 - NU*b^2 + b^2))/(6*a*b^2*(NU^2 - 1)) - (E*h*(2*b - 2*NU*b))/(6*a*b*(NU^2 - 1)),                                                                             0,       (E*h*(2*b - 2*NU*b))/(6*a*b*(NU^2 - 1)) + (E*h*(NU*b^2 + a^2 - b^2))/(6*a*b^2*(NU^2 - 1)),                                                                             0,   (E*h*(2*b - 2*NU*b))/(12*a*b*(NU^2 - 1)) - (E*h*(2*a^2 - NU*b^2 + b^2))/(12*a*b^2*(NU^2 - 1)),                                                                             0, - (E*h*(NU*b^2 + 4*a^2 - b^2))/(12*a*b^2*(NU^2 - 1)) - (E*h*(2*b - 2*NU*b))/(12*a*b*(NU^2 - 1));
 (2*E*h)/(3*a*(NU^2 - 1)) - (E*h*(NU*a^2 - a^2 + 4*b^2))/(12*a*b^2*(NU^2 - 1)),                                                                                               0,  (E*h*(a^2 - NU*a^2 + 2*b^2))/(6*a*b^2*(NU^2 - 1)) - (2*E*h)/(3*a*(NU^2 - 1)),                                                                                               0,      (E*h*(NU*a^2 - a^2 + b^2))/(6*a*b^2*(NU^2 - 1)) - (E*h)/(3*a*(NU^2 - 1)),                                                                                               0,   (E*h)/(3*a*(NU^2 - 1)) - (E*h*(a^2 - NU*a^2 + 2*b^2))/(12*a*b^2*(NU^2 - 1)),                                                                                               0;
                                                                             0,       (E*h*(2*b - 2*NU*b))/(6*a*b*(NU^2 - 1)) + (E*h*(NU*b^2 + a^2 - b^2))/(6*a*b^2*(NU^2 - 1)),                                                                             0,     (E*h*(2*a^2 - NU*b^2 + b^2))/(6*a*b^2*(NU^2 - 1)) - (E*h*(2*b - 2*NU*b))/(6*a*b*(NU^2 - 1)),                                                                             0, - (E*h*(NU*b^2 + 4*a^2 - b^2))/(12*a*b^2*(NU^2 - 1)) - (E*h*(2*b - 2*NU*b))/(12*a*b*(NU^2 - 1)),                                                                             0,   (E*h*(2*b - 2*NU*b))/(12*a*b*(NU^2 - 1)) - (E*h*(2*a^2 - NU*b^2 + b^2))/(12*a*b^2*(NU^2 - 1));
   (E*h)/(3*a*(NU^2 - 1)) - (E*h*(a^2 - NU*a^2 + 2*b^2))/(12*a*b^2*(NU^2 - 1)),                                                                                               0,      (E*h*(NU*a^2 - a^2 + b^2))/(6*a*b^2*(NU^2 - 1)) - (E*h)/(3*a*(NU^2 - 1)),                                                                                               0,  (E*h*(a^2 - NU*a^2 + 2*b^2))/(6*a*b^2*(NU^2 - 1)) - (2*E*h)/(3*a*(NU^2 - 1)),                                                                                               0, (2*E*h)/(3*a*(NU^2 - 1)) - (E*h*(NU*a^2 - a^2 + 4*b^2))/(12*a*b^2*(NU^2 - 1)),                                                                                               0;
                                                                             0,   (E*h*(2*b - 2*NU*b))/(12*a*b*(NU^2 - 1)) - (E*h*(2*a^2 - NU*b^2 + b^2))/(12*a*b^2*(NU^2 - 1)),                                                                             0, - (E*h*(NU*b^2 + 4*a^2 - b^2))/(12*a*b^2*(NU^2 - 1)) - (E*h*(2*b - 2*NU*b))/(12*a*b*(NU^2 - 1)),                                                                             0,     (E*h*(2*a^2 - NU*b^2 + b^2))/(6*a*b^2*(NU^2 - 1)) - (E*h*(2*b - 2*NU*b))/(6*a*b*(NU^2 - 1)),                                                                             0,       (E*h*(2*b - 2*NU*b))/(6*a*b*(NU^2 - 1)) + (E*h*(NU*b^2 + a^2 - b^2))/(6*a*b^2*(NU^2 - 1));
      (E*h*(NU*a^2 - a^2 + b^2))/(6*a*b^2*(NU^2 - 1)) - (E*h)/(3*a*(NU^2 - 1)),                                                                                               0,   (E*h)/(3*a*(NU^2 - 1)) - (E*h*(a^2 - NU*a^2 + 2*b^2))/(12*a*b^2*(NU^2 - 1)),                                                                                               0, (2*E*h)/(3*a*(NU^2 - 1)) - (E*h*(NU*a^2 - a^2 + 4*b^2))/(12*a*b^2*(NU^2 - 1)),                                                                                               0,  (E*h*(a^2 - NU*a^2 + 2*b^2))/(6*a*b^2*(NU^2 - 1)) - (2*E*h)/(3*a*(NU^2 - 1)),                                                                                               0;
                                                                             0, - (E*h*(NU*b^2 + 4*a^2 - b^2))/(12*a*b^2*(NU^2 - 1)) - (E*h*(2*b - 2*NU*b))/(12*a*b*(NU^2 - 1)),                                                                             0,   (E*h*(2*b - 2*NU*b))/(12*a*b*(NU^2 - 1)) - (E*h*(2*a^2 - NU*b^2 + b^2))/(12*a*b^2*(NU^2 - 1)),                                                                             0,       (E*h*(2*b - 2*NU*b))/(6*a*b*(NU^2 - 1)) + (E*h*(NU*b^2 + a^2 - b^2))/(6*a*b^2*(NU^2 - 1)),                                                                             0,     (E*h*(2*a^2 - NU*b^2 + b^2))/(6*a*b^2*(NU^2 - 1)) - (E*h*(2*b - 2*NU*b))/(6*a*b*(NU^2 - 1))];
                                                                         
end