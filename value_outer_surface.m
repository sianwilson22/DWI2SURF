%Interpolate outside white matter surface
function [val,val_o] = value_outer_surface(s, map)
    loc = s.coord+s.normal+1;
    val = interp3(map,loc(2,:),loc(1,:),loc(3,:));
    
    loc_o = s.coord+(s.normal*2)+1;
    val_o = interp3(map,loc_o(2,:),loc_o(1,:),loc_o(3,:));
end