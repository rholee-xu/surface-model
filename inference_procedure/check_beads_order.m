function [] = check_beads_order(subset_before,subset_after)
c = 1:size(subset_after,2);

quickscatter(subset_after);
text(subset_after(1,:),subset_after(2,:),subset_after(3,:),string(c))
hold on;
c1 = 1:size(subset_before,2);

quickscatter(subset_before);
text(subset_before(1,:),subset_before(2,:),subset_before(3,:),string(c1))

view([87.9 5.1])
end