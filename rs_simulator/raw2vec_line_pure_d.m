function cam_point_pram = raw2vec_line_pure_d(camera_set, point_set, line_set, param, is_gt)
    cam_point_pram = zeros(1, 12 * param.num_camera + 3 * param.num_point + 4 * param.num_line);
    for i = 1:param.num_camera
        id = (i - 1) * 12 + 1;
        if is_gt
            axang = rotm2axang(camera_set{i}.gt_oritation);
        else
            axang = rotm2axang(camera_set{i}.oritation);
        end
        if axang(4) == 0
            axang = [1 0 0 2 * pi];
        end
        axis_vec = axang(1:3) * axang(4);
        if is_gt
            cam_point_pram(id:id + 2) = axis_vec;
            cam_point_pram(id + 3: id + 5) = camera_set{i}.gt_translation';
            cam_point_pram(id + 6:id + 8) = camera_set{i}.w_gt';
            cam_point_pram(id + 9:id + 11) = camera_set{i}.d_gt';
        else
            cam_point_pram(id:id + 2) = axis_vec;
            cam_point_pram(id + 3: id + 5) = camera_set{i}.translation';
            cam_point_pram(id + 6:id + 8) = 1e-7 * ones(1,3);
            cam_point_pram(id + 9:id + 11) = 1e-7 * ones(1,3);
        end
    end
    for i = 1:param.num_point
        id = param.num_camera * 12 + (i - 1) * 3 + 1;
        if is_gt
            cam_point_pram(id: id + 2) = point_set.point_gt(i,:);
        else
            cam_point_pram(id: id + 2) = point_set.point(i,:);
        end
    end
    for i = 1:param.num_line
        id = param.num_camera * 12 + param.num_point * 3 + (i - 1) * 4 + 1;
        if is_gt
            cam_point_pram(id: id + 3) = line_set.line_gt(i,:);
        else
            cam_point_pram(id: id + 3) = line_set.line(i,:);
        end
    end

end