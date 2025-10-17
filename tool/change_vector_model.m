function new_vector = change_vector_model(vector, param, unit_dir, is_RT)
if is_RT == 1
    new_vector = zeros(1, 9+2*param.num_camera+3*param.num_point);
    new_vector(1,1:6) = vector(1,1:6);
    new_vector(1,7:9) = unit_dir;
    t1 = vector(:, 4: 6);
    for i = 1: param.num_camera
        ti = vector(:, 12*(i-1)+4: 12*(i-1)+6);
        distance = (ti-t1) * unit_dir';
        new_vector(8+2*i) = distance;

        di = vector(:, 12*(i-1)+10: 12*(i-1)+12);
        new_vector(8+2*i+1) = di * unit_dir';
    end
    new_vector(1, 10+2*param.num_camera:end) = vector(1, 1+12*param.num_camera: end);
else
    new_vector = zeros(1, 12*param.num_camera+3*param.num_point);
    for i = 1: param.num_camera
        new_vector(1, 12*(i-1)+1: 12*(i-1)+3) = vector(1,1:3);
        new_vector(1, 12*(i-1)+4: 12*(i-1)+6) = vector(1,4:6) + vector(8+2*i)*vector(1,7:9);
        new_vector(1, 12*(i-1)+10: 12*(i-1)+12) = vector(9+2*i)*vector(1,7:9);
    end
    new_vector(1, 1+12*param.num_camera: end) = vector(1, 10+2*param.num_camera:end);
end