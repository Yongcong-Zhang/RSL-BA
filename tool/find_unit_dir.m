function d = find_unit_dir(state_vector)
d1 = state_vector(:, 4: 6);
d2 = state_vector(:, 16: 18);
%d3 = state_vector(:, 28: 30);
d = (d2-d1)/norm(d2-d1);
%d = (d3-d1)/norm(d3-d1);
end