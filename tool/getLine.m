function corrected_line = getLine(refined_param_middle_jacobian, param)
    num_cameras = param.num_camera;
    num_points = param.num_point;
    num_lines = param.num_line;
    corrected_line = [];
    for i = 1:num_lines
        id = 12 * num_cameras + 3 * num_points + (i - 1) * 4 + 1;
        corrected_line = [corrected_line; refined_param_middle_jacobian(1, id:id + 3)];
    end
end