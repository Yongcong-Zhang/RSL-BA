function corrected_point = getPoint(refined_param_middle_jacobian, param)
    num_cameras = param.num_camera;
    num_points = param.num_point;
    corrected_point = [];
    for i = 1:num_points
        id = 12 * num_cameras + (i - 1) * 3 + 1;
        corrected_point = [corrected_point; refined_param_middle_jacobian(1, id:id + 2)];
    end
end