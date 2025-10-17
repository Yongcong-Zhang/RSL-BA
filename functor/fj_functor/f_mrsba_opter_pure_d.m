function [F] = f_mrsba_opter_pure_d(x, camera_set, param)
    cam_fir = x(1:10);
    cam_oth = x(11:end);
    K = param.K_matrix;
    num_camera = param.num_camera;
    num_point = param.num_point;
    F = zeros(2, num_camera * num_point);

    f_x = K(1, 1);
    f_y = K(2, 2);
    u_0 = K(1, 3);
    v_0 = K(2, 3);
    for camera_id = 1:param.num_camera
        for point_id = 1:param.num_point
            
            pose_state_id = camera_id;
            point_state_id = num_camera*2 + (point_id - 1) * 3;   
            measurement_id = (camera_id - 1) * num_point + point_id;
            point = cam_oth(1, point_state_id:point_state_id + 2);
            obs = camera_set{1, camera_id}.feature_point(point_id,:)';
            r = cam_fir(1, 1:3);
            if camera_id ~= 1
                t = cam_fir(1, 4:6) + cam_fir(1, 7:9)* cam_oth((pose_state_id-1)*2);
            else
                t = cam_fir(1, 4:6);
            end
            w = [0 0 0];
            d = cam_fir(1, 7:9)* cam_oth((pose_state_id-1)*2+1);
            R = axang2rotm([r / norm(r) norm(r)]);
            p_c = (eye(3) + X_(w) * obs(2)) * R * point' + t' + d' * obs(2);
            project_obs = [f_x;f_y] .* p_c(1:2) / p_c(3) + [u_0; v_0];
            F(1:2, measurement_id) = obs - project_obs;

            
        end
    end
F = reshape(F, 2 * size(F, 2), 1);
end