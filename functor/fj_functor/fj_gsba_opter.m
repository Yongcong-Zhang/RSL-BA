function [F,J] = fj_gsba_opter(x, camera_set, param)
    K = param.K_matrix;
    num_camera = param.num_camera;
    num_point = param.num_point;
    F = zeros(2 * num_camera * num_point, 1);
    J = zeros(2 * num_camera * num_point, num_camera * 6 + num_point * 3);
    f_x = K(1, 1);
    f_y = K(2, 2);
    u_0 = K(1, 3);
    v_0 = K(2, 3);
    for camera_id = 1:param.num_camera
        for point_id = 1:param.num_point
            
            pose_state_id = (camera_id - 1) * 6 + 1;
            point_state_id = (num_camera  * 6) + (point_id - 1) * 3 + 1;   
            measurement_id = (camera_id - 1) * num_point + point_id;

            camera_pose = x(1, pose_state_id:pose_state_id + 5);
            point = x(1, point_state_id:point_state_id + 2);

            obs = camera_set{camera_id}.feature_point(point_id,:)';

            r = camera_pose(1, 1:3);
            t = camera_pose(1, 4:6);
            R = axang2rotm([r / norm(r) norm(r)]);
            p_c = R * point' + t';
            project_obs = [f_x;f_y] .* p_c(1:2) / p_c(3) + [u_0; v_0];
            
            
            F((measurement_id - 1) * 2 + 1:(measurement_id - 1) * 2 + 2, 1) = project_obs - obs;


            j_e_pc = - [f_x / p_c(3), 0           , - f_x * p_c(1) / (p_c(3) * p_c(3));
                        0           , f_y / p_c(3), - f_y * p_c(1) / (p_c(3) * p_c(3))];

            j_pc_pw = R;
            j_pc_r = - X_(R * point');
            j_pc_t = eye(3);
            j_id = (measurement_id - 1) * 2 + 1;

            J(j_id : j_id + 1,pose_state_id : pose_state_id + 5) = j_e_pc * [j_pc_r, j_pc_t];
            J(j_id : j_id + 1, point_state_id : point_state_id + 2) = j_e_pc * j_pc_pw;

        end
    end
end