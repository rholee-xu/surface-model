function beads = get_beads(raw_beads)
    
    beads = table2array(raw_beads(:,[1,2,3]))';
    %beads_array = raw_beads;
    %tip_index = find(A3(1,:) == max(A3(1,:)));
    [~, order] = sort(beads(1,:));
    beads = beads(:,order);

end