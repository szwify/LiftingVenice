function [uhat]=Cylinder_Poro_LaplaceSolution_Mode1_ur(propObject,s,r)
%

B=propObject.B;
nu_u=propObject.nu_u;
nu=propObject.nu;
c=propObject.c_dif;
g=propObject.g;

beta=sqrt(s/c);
xi=beta*r;

uhat=-(r/(s*2*g))*...
    ((1-2*nu_u)*(1-nu)*besseli(0,beta)+2*(nu_u-nu)*besseli(1,xi)/xi)./...
    ((1-nu)*besseli(0,beta)-2*(nu_u-nu)*besseli(1,beta)/beta);
 



