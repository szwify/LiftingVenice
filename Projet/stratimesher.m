function [vert,etri,tria,tnum] = stratimesher(L, Hs)

    function [HH] = szfunc(PP, Hs)
        H = PP(:,2);
        HH = ones(size(PP,1),1);
        for k=1:size(PP,1)
            for j=1:(size(Hs,1)-1)
                if H(k)<= Hs(j,1) && H(k)>= Hs(j+1,1)
                    HH(k) = ((Hs(j,1)-H(k))*Hs(j+1,2) - (Hs(j+1,1)-H(k))*Hs(j,2))/(Hs(j,1)-Hs(j+1,1));
                    break;
                end
            end
        end
    end
    
    hfun = @szfunc;
    
    H1 = Hs(1,1);
    H2 = Hs(2,1);
    
    nod1 = [
       0,H1;
         L, H1;
         L, H2;
         0, H2
        ];
    edg1 = [
         1 ,  2 ;  2 ,  3
         3 ,  4 ;  4 ,  1
        ] ;
    edg1(:,3) = +1;
    
    edge = [edg1];
    node = [nod1];
    for i=2:(size(Hs,1)-1)
        
        H1 = Hs(i,1);
        H2 = Hs(i+1,1);

        nodtemp = [
           0,H1;
             L, H1;
             L, H2;
             0, H2
            ];
        edgtemp = [
             1 ,  2 ;  2 ,  3
             3 ,  4 ;  4 ,  1
            ] ;

        edgtemp(:,3) = i;

        edgtemp(:,1:2) = ...
        edgtemp(:,1:2)+size(node,1);
        edge = [edge; edgtemp];
        node = [node; nodtemp];
    end
    
    part{1} = [...
            find(edge(:,3) == 1)
            ];
    for i=2:(size(Hs,1)-1)
    part{1} = [ ...
        part{1} 
        find(edge(:,3) == i)
        ] ;
    part{i} = [ ...
        find(edge(:,3) == i)
        ] ;
    
    
    end
    edge = edge(:,1:2) ;
    
[vert,etri, tria,tnum] = refine2(node,edge,part,[], hfun, Hs) ; % do not touch
 
end
    