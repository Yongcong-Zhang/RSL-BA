function [F,J] = fj_nmrsba_opter(cam_oth, camera_set, param)
    K = param.K_matrix;
    num_camera = param.num_camera;
    num_point = param.num_point;
    F = zeros(2 * num_camera * num_point,1);
    J = zeros(2 * num_camera * num_point, num_camera * 12 + num_point * 3);

    f_x = K(1, 1);
    f_y = K(2, 2);
    u_0 = K(1, 3);
    v_0 = K(2, 3);
    for camera_id = 1:param.num_camera
        for point_id = 1:param.num_point
            
            pose_state_id = (camera_id - 1) * 12 + 1;
            point_state_id = (num_camera  * 12) + (point_id - 1) * 3 + 1;   
            measurement_id = (camera_id - 1) * num_point + point_id;

            camera_pose = cam_oth(1, pose_state_id:pose_state_id + 5);

            rs_camera_pose = cam_oth(1, pose_state_id + 6:pose_state_id + 11);
            point = cam_oth(1, point_state_id:point_state_id + 2);

            obs = camera_set{1, camera_id}.feature_point(point_id,:)';
            obs_normal = (obs - [u_0; v_0]) ./ [f_x;f_y];

            r = camera_pose(1, 1:3);
            t = camera_pose(1, 4:6);
            w = rs_camera_pose(1, 1:3);
            d = rs_camera_pose(1, 4:6);
            R = axang2rotm([r / norm(r) norm(r)]);
            p_c = (eye(3) + X_(w) * obs_normal(2)) * R * point' + t' + d' * obs_normal(2);
            project_obs = p_c(1:2) / p_c(3);
            
            j_id = (measurement_id - 1) * 2 + 1;
            F(j_id:j_id+ 1, 1) = obs_normal - project_obs;


                j_e_pc = - [1 / p_c(3), 0           , - 1 * p_c(1) / (p_c(3) * p_c(3));
                            0           , 1 / p_c(3), - 1 * p_c(2) / (p_c(3) * p_c(3))];
    
                j_pc_pw = (eye(3) + X_(w) * obs_normal(2)) * R;
                j_pc_r = -(eye(3) + X_(w) * obs_normal(2)) * X_(R * point');
                j_pc_t = eye(3);
                j_pc_w = - obs_normal(2) * X_(R * point');
                j_pc_d = obs_normal(2) * eye(3);
                

                    J(j_id : j_id + 1,pose_state_id : pose_state_id + 11) = j_e_pc * [j_pc_r, j_pc_t, j_pc_w, j_pc_d];
                    J(j_id : j_id + 1, point_state_id : point_state_id + 2) = j_e_pc * j_pc_pw;


            
        end
    end

end