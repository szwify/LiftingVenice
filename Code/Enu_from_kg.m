function [E,nu]=Enu_from_kg(k,g)

    nu=(1-2*(g)/(3*k))/(2+2*(g)/(3*k));
    E=(9*(k))/(1+3*(k)/(g));
 
end
