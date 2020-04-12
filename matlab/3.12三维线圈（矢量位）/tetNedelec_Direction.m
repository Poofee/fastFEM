function [si] = tetNedelec_Direction(x,y,z)
% Function: determine direction of Nedelec by coordinates
tet_n = [ 0 , 1 , 2 , 3  ;
    0 , 2 , 3 , 1 ;
    0 , 3 , 1 , 2 ;
    1 , 2 , 0 , 3 ;
    1 , 3 , 2 , 0 ;
    2 , 3 , 0 , 1 ];
tet_n = tet_n + 1;
si = zeros(6,1);

for i=1:6
    % try x-define
    ii = tet_n(i,1+1) ;
    jj = tet_n(i,0+1) ;
    pp = x(ii) - x(jj) ;
    
    if pp > 0.0
        si(i) =  1 ;
    elseif pp < 0.0
        si(i) = -1 ;
    else
        % try y-define
        pp = y(ii) - y(jj) ;
        if pp > 0.0
            si(i) =  1 ;
        elseif pp < 0.0
            si(i) = -1 ;
        else
            % try z-define
            pp = z(ii) - z(jj) ;
            if pp > 0.0
                si(i) =  1 ;
            elseif pp < 0.0
                si(i) = -1 ;
            else
                % error
                error('tetNedelec_Direction');
            end
        end
    end
end