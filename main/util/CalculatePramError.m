function [point_error_norm, r_error_norm, t_error_norm, w_error_norm, d_error_norm, l_r_error_norm, l_d_error_norm] = CalculatePramError(refined_state_vector, gt_state_vector, param)
    num_camera = param.num_camera;
    num_point = param.num_point;
    num_line = param.num_line;
    r_error = 0;
    t_error = 0;
    w_error = 0;
    d_error = 0;
    point_error = 0;
    line_r_error = 0;
    line_d_error = 0;

    for i = 1:num_camera
        camera_id = (i - 1) * 12 + 1;
        camera_rota = refined_state_vector(1, camera_id:camera_id + 2);
        gt_camera_rota = gt_state_vector(1, camera_id:camera_id + 2);
        rs_camera_rota = refined_state_vector(1, camera_id + 6:camera_id + 8);
        gt_rs_camera_rota = gt_state_vector(1, camera_id + 6:camera_id + 8);
        camera_tran = refined_state_vector(1, camera_id + 3:camera_id + 5);
        gt_camera_tran = gt_state_vector(1, camera_id + 3:camera_id + 5);
        rs_camera_tran = refined_state_vector(1, camera_id + 9:camera_id + 11);
        gt_rs_camera_tran = gt_state_vector(1, camera_id + 9:camera_id + 11);


        tmp_error_r = axang2rotm([camera_rota / norm(camera_rota) norm(camera_rota)]);
        tmp_gt_error_r = axang2rotm([gt_camera_rota / norm(gt_camera_rota) norm(gt_camera_rota)]);
        tmp_rs_error_r = eye(3) + X_(rs_camera_rota);
        tmp_gt_rs_error_r = eye(3) + X_(gt_rs_camera_rota);
        tmp_error_r = rotm2axang((tmp_error_r)' * (tmp_gt_error_r));
        tmp_error_r = tmp_error_r(1:3) * tmp_error_r(4);
        tmp_error_w = rotm2axang((tmp_rs_error_r)' * (tmp_gt_rs_error_r));
        tmp_error_w = tmp_error_w(1:3) * tmp_error_w(4);
        tmp_error_w = abs(tmp_error_w) / norm(gt_rs_camera_rota);
        tmp_error_t = gt_camera_tran - camera_tran;
        tmp_error_d = abs(gt_rs_camera_tran - rs_camera_tran) / norm(gt_rs_camera_tran) ;
        r_error = r_error + tmp_error_r * tmp_error_r';
        t_error = t_error + tmp_error_t * tmp_error_t';
        w_error = w_error + tmp_error_w * tmp_error_w';
        d_error = d_error + tmp_error_d * tmp_error_d';
    end
    if param.is_error_p == 1
        for i = 1:num_point
            point_id = num_camera * 12 + (i - 1) * 3 + 1;
            point = refined_state_vector(1, point_id:point_id + 2);
            gt_point = gt_state_vector(1, point_id:point_id + 2);
            point_error = point_error + (point - gt_point) * (point - gt_point)';
        end
    end

    if param.is_error_l == 1
        for i = 1: num_line
            line_id = num_camera * 12 + num_point * 3 + (i - 1) * 4 + 1;
            line = refined_state_vector(1, line_id:line_id + 3);
            gt_line = gt_state_vector(1, line_id:line_id + 3);
            % line_rotm = eul2rotm(line(1:3), 'XYZ');
            % gt_line_rotm = eul2rotm(gt_line(1:3), 'XYZ');
            %line_r_error = line_r_error + acos((trace(line_rotm'*gt_line_rotm)-1)/2);
            gt_plucker = orth_to_plucker(gt_line);
            line_plucker = orth_to_plucker(line);
            gt_n = [gt_plucker(2,3), gt_plucker(1,3), gt_plucker(1,2)];
            gt_d = [gt_plucker(1,4), gt_plucker(2,4), gt_plucker(3,4)];
            line_n = [line_plucker(2,3), line_plucker(1,3), line_plucker(1,2)];
            line_d = [line_plucker(1,4), line_plucker(2,4), line_plucker(3,4)];
            %line_d_error = line_d_error + abs( norm(gt_n)/norm(gt_d) -  norm(line_n)/norm(line_d));
            
            e_s = gt_d*line_d'/(norm(gt_d)*norm(line_d ));
            if e_s < 0
                e_s = -e_s;
            end
            if e_s > 1
                e_s = 1;
            end
            line_r_error = line_r_error + acos(e_s);

            if abs(gt_plucker(1,4)) > abs(gt_plucker(2,4)) && abs(gt_plucker(1,4)) > abs(gt_plucker(3,4))
                P_gt = gt_plucker*[1; 0; 0; 0];
                P_line = line_plucker*[1; 0; 0; 0];
            elseif abs(gt_plucker(2,4)) > abs(gt_plucker(1,4)) && abs(gt_plucker(2,4)) > abs(gt_plucker(3,4))
                P_gt = gt_plucker*[0; 1; 0; 0];
                P_line = line_plucker*[0; 1; 0; 0];
            elseif abs(gt_plucker(3,4)) > abs(gt_plucker(1,4)) && abs(gt_plucker(3,4)) > abs(gt_plucker(2,4))
                P_gt = gt_plucker*[0; 0; 1; 0];
                P_line = line_plucker*[0; 0; 1; 0];  
            else
                P_gt = gt_plucker*[0; 1; 0; 0];
                P_line = line_plucker*[0; 1; 0; 0];
            end
            P_gt = P_gt./P_gt(4);
            P_gt = P_gt(1:3);
            P_line = P_line./P_line(4);
            P_line = P_line(1:3);    
            d_p1p2 = P_gt - P_line;
            proj_v1_d = dot(d_p1p2, gt_d) / norm(gt_d)^2 * gt_d;
            proj_v2_d = dot(d_p1p2, line_d) / norm(line_d)^2 * line_d;
            shortest_distance = norm(proj_v1_d - proj_v2_d);
            line_d_error = line_d_error + shortest_distance;
            if isnan(shortest_distance)
                disp('sdfsf')
            end
            %line_d_error = line_d_error + abs(line(4) - gt_line(4));
        end
    end

    point_error_norm = point_error / num_point;
    r_error_norm = r_error / num_camera;
    t_error_norm = t_error / num_camera;
    w_error_norm = w_error / num_camera;
    d_error_norm = d_error / num_camera;
    l_r_error_norm = line_r_error./num_camera;
    l_d_error_norm = line_d_error./num_camera;

end