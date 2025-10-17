function [F,J] = fj_nwrsba_opter(x,camera_set,pram)

    num_points = pram.num_point;
    num_camera = pram.num_camera;
    num_measurement = num_points * num_points;
    K = pram.K_matrix;
    F = zeros(2 * num_measurement, 1);
    J = zeros(2 * num_measurement, num_camera * 12 + num_points * 3 - 6);
    f_x = K(1, 1);
    f_y = K(2, 2);
    u_0 = K(1, 3);
    v_0 = K(2, 3);
 %--------------------------------------------------------------------------   
    for i=1:num_camera 
        for j = 1:num_points
            
            pose_id = (i - 1) * 12 + 1;
            point_id = (num_camera  * 12) + (j - 1) * 3 + 1;   
            measurement_id = (i - 1) * num_points + j;

            camera_pose = x(1, pose_id:pose_id + 5);
            rs_camera_pose = x(1, pose_id + 6:pose_id + 11);
            point_pose = x(1, point_id:point_id + 2);
            obs = camera_set{i}.feature_point(j,:)';



            vi = obs(2);
            X = point_pose(1,1);
            Y = point_pose(1,2);
            Z = point_pose(1,3);
            fx = f_x;
            fy = f_y;
            rx = camera_pose(1,1);
            ry = camera_pose(1,2);
            rz = camera_pose(1,3);
            tx = camera_pose(1,4);
            ty = camera_pose(1,5);
            tz = camera_pose(1,6);
            w_x = rs_camera_pose(1,1);
            w_y = rs_camera_pose(1,2);
            w_z = rs_camera_pose(1,3);
            d_x = rs_camera_pose(1,4);
            d_y = rs_camera_pose(1,5);
            d_z = rs_camera_pose(1,6);
            u0 = u_0;
            v0 = v_0;
            u_obs = obs(1);
            v_obs = obs(2);
            j_id = (measurement_id - 1) * 2 + 1;

            F(j_id:j_id + 1, 1) = [1 / fx 0;0 1 / fy] * error_obs_weighted_normal(X,Y,Z,d_x,d_y,d_z,fx,fy,rx,ry,rz,tx,ty,tz,u0,u_obs,v0,v_obs,vi,w_x,w_y,w_z);

                tmp_j = [1 / fx 0;0 1 / fy] * J_weighted_normal(X,Y,Z,d_x,d_y,d_z,fx,fy,rx,ry,rz,tx,ty,tz,v0,v_obs,vi,w_x,w_y,w_z);
                J(j_id:j_id + 1,pose_id:pose_id + 11) = tmp_j(:,1:12);
                J(j_id:j_id + 1, point_id:point_id + 2) = tmp_j(:,13:15);

            
         
        end
    end
end