function [eq_free,fix_nonZero,eq_fix]=PrepareDirichletBC(Imp_dirichlet,ID_array)

% Dirichlet BC .....
% Imp_dirichlet is a  matrix with columms:  node, dof, imposed value

if (isempty(Imp_dirichlet))
    eq_free=ID_array(:);
    fix_nonZero=[];
    eq_fix=[];
else
    eq_fix=[];
    for imp=1:length(Imp_dirichlet(:,1))
        eq_fix=[eq_fix ;  ID_array(Imp_dirichlet(imp,1),Imp_dirichlet(imp,2)) ];
    end
    
    %disp( 'fix eq.'); disp(length(eq_fix(:))); disp(' --');
    
    eq_free=setdiff(ID_array(:),eq_fix(:));
    
    fix_nonZero=unique(find(Imp_dirichlet(:,3)~=0));
    
end

end