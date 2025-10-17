function [camera_set, point_set,line_set, param, flag] = rs_simulate_line_pure_d(param)
cameraParams = cameraParameters('IntrinsicMatrix', param.K_matrix');
sigma = [param.measurement_noise 0; 0 param.measurement_noise];
R_sig = chol(sigma);
point_set.point_gt = generate3DCube(1.5);

% point_set.point_gt = point_set.point_gt + ones(size(point_set.point_gt))*0.5;
r_p3d = [1 1 1, pi/18];
R_p3d = axang2rotm(r_p3d);
point_set.point_gt = R_p3d*point_set.point_gt';
point_set.point_gt = point_set.point_gt';

param.num_point = size(point_set.point_gt, 1);
% ------------- camera and feature generation -------------
camera_set = cell(1, param.num_camera);

% ------ first camera ------
w = 2 * (rand(3,1) - 0.5 * ones(3, 1));
d = 2 * (rand(3,1) - 0.5 * ones(3, 1));

a=0;b=1;c=0;
d = [a; b; c];

w_gt = [0; 0; 0];
d_gt = randn(1) * d / norm(d) *(param.d_norm * param.t_v);

camera_pose_1.w_gt = w_gt;
camera_pose_1.d_gt = d_gt;
r = [1 0 0, pi];
R = axang2rotm(r);
t = [0; 0; 13];


camera_pose_1.gt_oritation = R;
camera_pose_1.gt_translation = t;
[Rot_Rows, trans_Rows] = linearEgoMotion(R,t,w_gt,d_gt,cameraParams);
[p2d_1,p3d_RS, flag] = RSWorld2Image(point_set.point_gt, Rot_Rows, trans_Rows, cameraParams);
camera_pose_1.gt_feature_point = p2d_1;
noise = randn(param.num_point, 2) * R_sig;
p2d_1 = p2d_1 + noise;
camera_pose_1.feature_point = p2d_1;
camera_pose_1.oritation = camera_pose_1.gt_oritation;
camera_pose_1.translation = camera_pose_1.gt_translation;
camera_set{1, 1} = camera_pose_1;


axis_rotate = zeros(param.num_camera - 1, 2);
for j = 1:param.num_camera - 1
    r = [0 0 1 2 * pi * (j) / (param.num_camera - 1)];
    R = axang2rotm(r);
    new_axis = R * [1; 0; 0];
    axis_rotate(j, :) = new_axis(1:2)';
end

AbsolutePose = [rigid3d(camera_pose_1.oritation, - camera_pose_1.translation' * camera_pose_1.oritation)];
for i = 2:param.num_camera

    r = [1 0 0, pi];
    R = axang2rotm(r);
    %t = [0 ; 0 ; (13 - (i - 1) * randn(1) * 3)];
    t_scale = 0.7;
    t = abs(randn(1))*[(i-1)*a ; (i-1)*b ; (i-1)*c].*t_scale + [0; 0; 13];


    w_gt = [0; 0; 0];
    d_gt = randn(1) * d / norm(d) * (param.d_norm * param.t_v);


    camera_pose.w_gt = w_gt;
    camera_pose.d_gt = d_gt;

    camera_pose.gt_oritation = R;
    camera_pose.gt_translation = t;
    [Rot_Rows, trans_Rows] = linearEgoMotion(R,t,w_gt,d_gt,cameraParams);
    [p2d,p3d_RS, flag_rs] = RSWorld2Image(point_set.point_gt, Rot_Rows, trans_Rows, cameraParams);
    camera_pose.gt_feature_point = p2d;
    noise = randn(param.num_point, 2) * R_sig;
    p2d = p2d + noise;
    camera_pose.feature_point = p2d;
    flag = flag || flag_rs;
    camera_set{1, i} = camera_pose;


    E = estimateEssentialMatrix(camera_pose_1.feature_point, camera_set{1, i}.feature_point, cameraParams);
    [relativeOrientation,relativeLocation] = relativeCameraPose(E,cameraParams,camera_pose_1.feature_point,camera_set{1, i}.feature_point);

    scale = norm(camera_pose_1.gt_oritation' * camera_pose_1.gt_translation - camera_set{1, i}.gt_oritation' * camera_set{1, i}.gt_translation);
    camera_set{1, i}.oritation = relativeOrientation * camera_pose_1.gt_oritation;
    camera_set{1, i}.translation = ((camera_pose_1.translation' - scale * relativeLocation) * relativeOrientation')';
    absolute_pose = rigid3d(camera_set{1, i}.oritation, - camera_set{1, i}.translation' * camera_set{1, i}.oritation);
    AbsolutePose = [AbsolutePose;absolute_pose];
end


tracks = [];
for i = 1:param.num_point
    feature_uv = [];
    for j = 1:param.num_camera
        feature_uv = [feature_uv; camera_set{1, j}.feature_point(i, :)];
    end
    track_d = pointTrack((1:param.num_camera), feature_uv);
    tracks = [tracks, track_d];
end

ViewId = uint32(1:param.num_camera)';
camPoses = table(ViewId, AbsolutePose);

focalLength = [param.fx, param.fy];
principalPoint = [param.cx, param.cy];
imageSize = [param.H, param.W];
intrinsics_matrix = cameraIntrinsics(focalLength, principalPoint, imageSize);

noise_sigma = 0.1;
point_set.point = point_set.point_gt + noise_sigma * randn(size(point_set.point_gt));



line_set.line_pair = [[1,13];[4,16];[41,53];[44,56];[1,4];[13,16];[41,44];[53,56];[1,41];[13,53];[4,44];[16,56]];

line_set.line_gt = [];
line_set.line = [];
for j = 1: size(camera_set,2)
    camera_set{j}.gt_feature_line = [];
    camera_set{j}.feature_line = [];
    camera_set{j}.gt_feature_points_on_line = [];
    camera_set{j}.feature_points_on_line = [];
end

for i = 1: size(line_set.line_pair,1)
    j1 = line_set.line_pair(i,1);
    j2 = line_set.line_pair(i,2);
    P1 = [point_set.point_gt(j1, :), 1];
    P2 = [point_set.point_gt(j2, :), 1];
    L = P1'*P2-P2'*P1;
    L = L./norm(L);
    orth = plucker_to_orth(L);
    line_set.line_gt = [line_set.line_gt; orth];

    for j = 1: size(camera_set,2)
        R_line = camera_set{j}.gt_oritation;
        T_line = camera_set{j}.gt_translation;
        w_line = camera_set{j}.w_gt;
        d_line = camera_set{j}.d_gt;
        K_line = param.K_matrix;
        line_2d = get_project_line(K_line, L, R_line, T_line, w_line, d_line,camera_set{j}.gt_feature_point(line_set.line_pair(i,1), 2),camera_set{j}.gt_feature_point(line_set.line_pair(i,2), 2));
        camera_set{j}.gt_feature_line = [camera_set{j}.gt_feature_line; line_2d];
    end

end


for i = 1: size(line_set.line_pair,1)
    j1 = line_set.line_pair(i,1);
    j2 = line_set.line_pair(i,2);
    P1 = [point_set.point(j1, :), 1];
    P2 = [point_set.point(j2, :), 1];
    L = P1'*P2-P2'*P1;
    L = L./norm(L);
    orth = plucker_to_orth(L);
    line_set.line = [line_set.line; orth];

    for j = 1: size(camera_set,2)
        R_line = camera_set{j}.gt_oritation;
        T_line = camera_set{j}.gt_translation;
        w_line = camera_set{j}.w_gt;
        d_line = camera_set{j}.d_gt;
        K_line = param.K_matrix;
        line_2d = get_project_line(K_line, L, R_line, T_line, w_line, d_line,camera_set{j}.gt_feature_point(line_set.line_pair(i,1), 2),camera_set{j}.gt_feature_point(line_set.line_pair(i,2), 2));
        camera_set{j}.feature_line = [camera_set{j}.feature_line; line_2d];
    end

end

num_error_of_line = param.num_error_of_line;


sigma = [param.measurement_noise 0; 0 param.measurement_noise];
%sigma = [0.1 0; 0 0.1];
R_sig = chol(sigma);
R_sig(2,2) = 0;


for camera_id = 1:param.num_camera
    for line_id = 1:size(line_set.line_pair,1)


        obs = camera_set{1, camera_id}.gt_feature_line(line_id,:);
        v1 = obs(8);
        v2 = obs(9);
        d_v = (v2- v1)/(num_error_of_line-1);
        %line_2d = [A3(1,3), A3(3,2), A2(1,3)+A3(2,1), A2(3,2), A1(1,3)+A2(2,1), A1(3,2), A1(2,1), v1, v2];


        for s = 1: num_error_of_line
            v_measure = obs(8) + (s-1) * d_v;
            u_measure = -(obs(1)*v_measure^3+obs(3)*v_measure^2+obs(5)*v_measure+obs(7))/(obs(2)*v_measure^2+obs(4)*v_measure+obs(6));
            slope_measure1 = 3*obs(1)*v_measure^2+2*obs(2)*u_measure*v_measure+2*obs(3)*v_measure+obs(4)*u_measure+obs(5);
            slope_measure2 = obs(2)*v_measure^2 + obs(4)*v_measure + obs(6);
            slope_measure_u = slope_measure1/slope_measure2;
            slope_measure_v = slope_measure2/slope_measure1;
            camera_set{camera_id}.gt_feature_points_on_line = [camera_set{camera_id}.gt_feature_points_on_line; [u_measure, v_measure, slope_measure_u, slope_measure_v]];
        end

    end

    points_on_line = camera_set{camera_id}.gt_feature_points_on_line(:,1:2);
    noise = randn(size(camera_set{camera_id}.gt_feature_points_on_line, 1),2) * R_sig;
    camera_set{camera_id}.feature_points_on_line = camera_set{camera_id}.gt_feature_points_on_line;
    camera_set{camera_id}.feature_points_on_line(:,1:2)= points_on_line + noise;
end

param.num_line = size(line_set.line_pair, 1);


end