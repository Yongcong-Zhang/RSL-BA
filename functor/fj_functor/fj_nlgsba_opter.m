function [F, J]  = fj_nlgsba_opter(x, camera_set, param)
K = param.K_matrix;
num_camera = param.num_camera;
num_point = param.num_point;
num_line = param.num_line;
num_measurement = num_point * num_camera;

f_x = K(1, 1);
f_y = K(2, 2);
u_0 = K(1, 3);
v_0 = K(2, 3);


num_error_of_line = param.num_error_of_line;
F_lines = zeros(1, num_camera * num_line);
J_lines = zeros(num_camera * num_line , num_camera * 6 + num_line * 4);

line_error_id = 1;
for camera_id = 1:param.num_camera
    for line_id = 1:num_line


        pose_state_id = (camera_id - 1) * 6 + 1;
        line_state_id = (num_camera  * 6) + (num_point) * 3 + (line_id - 1) * 4 + 1;
        measurement_id = (camera_id - 1) * num_line + line_id;

        camera_pose = x(1, pose_state_id:pose_state_id + 5);

        line_orth = x(1,line_state_id:line_state_id + 3);
        line_plucker = orth_to_plucker(line_orth);
        line_plucker = line_plucker./norm(line_plucker);

        r = camera_pose(1, 1:3);
        t = camera_pose(1, 4:6);
        R = axang2rotm([r / norm(r) norm(r)]);


        
        pol_id1 = (line_id-1)*num_error_of_line+1;
        u_measure1 = camera_set{camera_id}.feature_points_on_line(pol_id1, 1);
        v_measure1 = camera_set{camera_id}.feature_points_on_line(pol_id1, 2);

        pol_id2 = (line_id-1)*num_error_of_line+num_error_of_line;
        u_measure2 = camera_set{camera_id}.feature_points_on_line(pol_id2, 1);
        v_measure2 = camera_set{camera_id}.feature_points_on_line(pol_id2, 2);


        u_measure1 = (u_measure1 - u_0)/f_x;
        v_measure1 = (v_measure1 - v_0)/f_y;

        u_measure2 = (u_measure2 - u_0)/f_x;
        v_measure2 = (v_measure2 - v_0)/f_y;

        P0 = [R, t'];
        m_l = P0 * line_plucker * P0';
        mx = -m_l(2,3);
        my = m_l(1,3);
        mz = -m_l(1,2);
        ml = [mx; my; mz];
        A = [[u_measure1, v_measure1, 1]; [u_measure2, v_measure2, 1]];
        B = (1/(3*(mx^2+my^2))).*[[1 0.5];[0.5 1]];
        error = ml'*(A'*B*A)*ml;

        F_lines(line_error_id) = error;


        orth1 = line_orth(1); orth2 = line_orth(2); orth3 = line_orth(3); orth4 = line_orth(4);tx=t(1);ty=t(2);tz=t(3);
        r11 = R(1,1); r12 = R(1,2); r13 = R(1,3);
        r21 = R(2,1); r22 = R(2,2); r23 = R(2,3);
        r31 = R(3,1); r32 = R(3,2); r33 = R(3,3);
        J_0 = J_full_nlgsba_func(orth1,orth2,orth3,orth4,r11,r12,r13,r21,r22,r23,r31,r32,r33,tx,ty,tz,u_measure1,u_measure2,v_measure1,v_measure2);
        
        J_lines(line_error_id,pose_state_id : pose_state_id + 5) = J_0(:,1:6);
        J_lines(line_error_id, line_state_id : line_state_id + 3) = J_0(:,7:10);

                line_error_id = line_error_id + 1;
        


    end
end

F_lines = reshape(F_lines, size(F_lines, 2), 1);
F = F_lines;
J = J_lines;
end