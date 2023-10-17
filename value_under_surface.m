function [val,val_u] = value_under_surface(s, map)
    loc = s.coord-s.normal+1;
    val = interp3(map,loc(2,:),loc(1,:),loc(3,:));
    
    loc_u = s.coord-(s.normal*2)+1;
    val_u = interp3(map,loc_u(2,:),loc_u(1,:),loc_u(3,:));
end

%Interpolate outside white matter surface
%function [val] = value_under_surface(s, map)
    %loc = s.coord+s.normal+1;
    %val = interp3(map,loc(1,:),loc(2,:),loc(3,:));
%end


%Interpolate inside white matter surface- multiple steps
%function [val] = value_under_surface(s, map)
    %loc = s.coord+s.normal+1;
    %val = interp3(map,loc(1,:),loc(2,:),loc(3,:));
    %val2 = interp3(map,val(1,:),val(2,:),val(3,:));
%end