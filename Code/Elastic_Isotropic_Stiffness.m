function L_elas=Elastic_Isotropic_Stiffness(k,g,Geo)

% create L_elas

La=k+(4./3.)*g;
Lb=k-(2./3.)*g;

switch Geo
    case 'PlaneStrain'
        
        L_elas=[La Lb 0. ;...
            Lb La 0.;...
            0 0  g];
   
    otherwise
        
        L_elas=[La Lb Lb 0 0 0 ;...
            Lb La Lb 0 0 0; ...
            Lb Lb La 0 0 0 ; ...
            0 0 0 2*g 0 0 ;...
            0 0 0 0 2*g 0 ;...
            0 0 0 0 0 2*g];
        
end
end