function [phat]=Cylinder_Poro_LaplaceSolution_Mode1_p(propObject,s,r)
%

B=propObject.B;
nu_u=propObject.nu_u;
nu=propObject.nu;
c=propObject.c_dif;

beta=sqrt(s/c);
xi=beta*r;

phat=(1/s)*(2./3.)*B*(1+nu_u)*...
    ((1-nu)*(besseli(0,beta)-besseli(0,xi))./...
    ((1-nu)*besseli(0,beta)-2*(nu_u-nu)*besseli(1,beta)/beta));


