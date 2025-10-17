function [camera_set, point_set, param, flag] = rs_simulate_degeneracy(param)
cameraParams = cameraParameters('IntrinsicMatrix', param.K_matrix');
sigma = [param.measurement_noise 0; 0 param.measurement_noise];
R_sig = chol(sigma);
point_set.point_gt = generate3DCube(1.5);
param.num_point = size(point_set.point_gt, 1);
camera_set = cell(1, param.num_camera);

% ------------- first camera -------------
w = 2 * (rand(3,1) - 0.5 * ones(3, 1));
d = 2 * (rand(3,1) - 0.5 * ones(3, 1));
w_gt = w / norm(w) * (pi * param.w_norm  / 180) * param.t_v;
d_gt = d / norm(d) * param.d_norm * param.t_v;
camera_pose_1.w_gt = w_gt;
camera_pose_1.d_gt = d_gt;
r = [1 0 0, pi];
R = axang2rotm(r);
t = [0; 0; param.scale]; 
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










AbsolutePose = [rigid3d(camera_pose_1.oritation, - camera_pose_1.translation' * camera_pose_1.oritation)];
for i = 2:param.num_camera
    w = 2 * (rand(3,1) - 0.5 * ones(3, 1));
    d = 2 * (rand(3,1) - 0.5 * ones(3, 1));
    w_gt = w / norm(w) * (pi * param.w_norm / 180) * param.t_v;
    d_gt = d / norm(d) * param.d_norm * param.t_v;
    camera_pose.w_gt = w_gt;
    camera_pose.d_gt = d_gt;
    
    r = [1 0 0, pi];
    R = axang2rotm(r);
    r = [0 0 1 param.pitch * pi * i / param.num_camera];
    R = axang2rotm(r) * R;
    r = [0 1 0 (rand(1)) * 2 * pi];
    R = axang2rotm(r) * R;

    t = [0; 0; param.scale]; 
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
xyzPoints = triangulateMultiview(tracks, camPoses, intrinsics_matrix);
point_set.point = xyzPoints;



end