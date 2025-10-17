function [F, J]  = fj_nslrsba_opter(x, camera_set, param)
K = param.K_matrix;
num_camera = param.num_camera;
num_point = param.num_point;
num_line = param.num_line;
f_x = K(1, 1);
f_y = K(2, 2);
u_0 = K(1, 3);
v_0 = K(2, 3);


num_error_of_line = param.num_error_of_line;
F_lines = zeros(num_error_of_line, num_camera * num_line);
F_slope = zeros(num_error_of_line, num_camera * num_line);
J_lines = zeros(num_camera * num_line * num_error_of_line, num_camera * 12 + num_point * 3 + num_line * 4);
J_slope = zeros(num_camera * num_line * num_error_of_line, num_camera * 12 + num_point * 3 + num_line * 4);
line_error_id = 1;
for camera_id = 1:param.num_camera
    for line_id = 1:num_line

        pose_state_id = (camera_id - 1) * 12 + 1;
        line_state_id = (num_camera  * 12) + (num_point) * 3 + (line_id - 1) * 4 + 1;
        measurement_id = (camera_id - 1) * num_line + line_id;


        camera_pose = x(1, pose_state_id:pose_state_id + 5);
        rs_camera_pose = x(1, pose_state_id + 6:pose_state_id + 11);
        line_orth = x(1,line_state_id:line_state_id + 3);
        line_plucker = orth_to_plucker(line_orth);
        line_plucker = line_plucker./norm(line_plucker);


        r = camera_pose(1, 1:3);
        t = camera_pose(1, 4:6);
        w = rs_camera_pose(1, 1:3);
        d = rs_camera_pose(1, 4:6);
        R = axang2rotm([r / norm(r) norm(r)]);

        P0 = [R, t'];
        Q = [X_(w)*R, d'];
        A1 = P0*line_plucker*P0';
        A2 = P0*line_plucker*Q' + Q*line_plucker*P0';
        A3 = Q*line_plucker*Q';


        l = [A3(1,3), A3(3,2), A2(1,3), A3(2,1), A2(3,2), A1(1,3), A2(2,1), A1(3,2), A1(2,1)];
        n = [-line_plucker(2,3); line_plucker(1,3); -line_plucker(1,2)];
        a = [line_plucker(1,4); line_plucker(2,4); line_plucker(3,4)];

        for s = 1: num_error_of_line
            pol_id = (line_id-1)*num_error_of_line+s;
            u_measure = camera_set{camera_id}.feature_points_on_line(pol_id, 1);
            v_measure = camera_set{camera_id}.feature_points_on_line(pol_id, 2);
            u_measure = (u_measure - u_0)/f_x;
            v_measure = (v_measure - v_0)/f_y;
            v_measures = v_measure;
            u_measures = u_measure;
            ld1 = A3(1,3)*v_measure^3+A3(3,2)*u_measure*v_measure^2+(A2(1,3)+A3(2,1))*v_measure^2+A2(3,2)*u_measure*v_measure+(A1(1,3)+A2(2,1))*v_measure+A1(3,2)*u_measure+A1(2,1);
            ld2 = ((A1(3,2)+A2(3,2)*v_measure+A3(3,2)*v_measure^2)^2+(A1(1,3)+A2(1,3)*v_measure+A3(1,3)*v_measure^2)^2)^0.5;
            F_lines(s, measurement_id) = -ld1/ld2;

            slope_measure_u = camera_set{camera_id}.feature_points_on_line(pol_id, 3);
            slope_measure_v = camera_set{camera_id}.feature_points_on_line(pol_id, 4);
            %u_projects = -(A3s(1,3)*v_measures^3+(A2s(1,3)+A3s(2,1))*v_measures^2+(A1s(1,3)+A2s(2,1))*v_measures+A1s(2,1))/(A3s(3,2)*v_measures^2+A2s(3,2)*v_measures+A1s(3,2));
            u_projects = u_measures;
            slope_project1 = 3*A3(1,3)*v_measures^2+2*A3(3,2)*u_projects*v_measures+2*(A2(1,3)+A3(2,1))*v_measures+A2(3,2)*u_projects+(A2(2,1)+A1(1,3));
            slope_project2 = A3(3,2)*v_measures^2 + A2(3,2)*v_measures + A1(3,2);
            slope_project_u = slope_project1/slope_project2;
            slope_project_v = slope_project2/slope_project1;
            sm = [slope_measure_u; slope_measure_v];
            sp = [slope_project_u; slope_project_v];
            
            e_s = sm'*sp/(norm(sm)*norm(sp));

            if e_s < 0
                aaa = slope_measure_u;
                slope_measure_u = -slope_measure_u;
                slope_measure_v = -slope_measure_v ;
                e_s = -e_s;
            end

            if e_s > 1
                e_s = 1;
            end
            F_slope(s, measurement_id) = acos(e_s);
            [J_0, J_1] = get_line_J_norm_nsl1(u_measure, v_measure ,l, R, t, w, d, n, a, line_orth,slope_measure_u, slope_measure_v);
            J_lines(line_error_id,pose_state_id : pose_state_id + 11) = J_0(:,1:12);
            J_lines(line_error_id, line_state_id : line_state_id + 3) = J_0(:,13:16);
            J_slope(line_error_id,pose_state_id : pose_state_id + 11) = J_1(:,1:12);
            J_slope(line_error_id, line_state_id : line_state_id + 3) = J_1(:,13:16);

            line_error_id = line_error_id + 1;

        end
    end
end
F_lines = reshape(F_lines, param.num_error_of_line * size(F_lines, 2), 1);
F_slope = reshape(F_slope, param.num_error_of_line * size(F_slope, 2), 1);
F = [F_lines; F_slope];
J = [J_lines; J_slope];

end