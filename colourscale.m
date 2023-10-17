%% Make colourscale bar from negative (blue) to positive (red)
c = hot;
c = flipud(c);
s = parula;
full = [s;c];
full(240:282,:) = ones(43,3);

full_smooth = smoothdata(full);
%%