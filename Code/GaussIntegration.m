% Gaussian integration over unit elements
function [Kel]=GaussIntegration(fun,dim,IntegrationOrder)

switch dim
    case 1
       
        %one dimensional integration
        [xeta,Wl]=GaussQuadratureCoefs(IntegrationOrder);
        
        for i=1:length(Wl)
            
            Ra=Wl(i)*fun(xeta(i));
            
            if i==1
                Kel=Ra;
            else
                Kel=Kel+Ra;
            end
            
        end
        
    case 2
        % 2 D integration ...
        
        [xeta,Wl]=GaussQuadratureCoefs(IntegrationOrder);
        
        k=0;
        for i=1:IntegrationOrder
            for j=1:IntegrationOrder
                k=k+1;
                
                Ra = Wl(i)*Wl(j)*fun([xeta(i), xeta(j)]);
                
                if k==1
                    Kel=Ra;
                else
                    Kel=Kel+Ra;
                end
            end
        end
        
end


end

