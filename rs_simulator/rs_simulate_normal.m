function [camera_set_mid, param] = rs_simulate_normal(camera_set, param)

camera_set_mid = camera_set;
cameraParams = cameraParameters('IntrinsicMatrix',param.K_matrix');


for i = 1:param.num_camera
R = camera_set_mid{i}.gt_oritation;
t = camera_set_mid{i}.gt_translation;
w_gt = camera_set_mid{i}.w_gt;
d_gt = camera_set_mid{i}.d_gt;
[Rot_Rows, trans_Rows] = linearEgoMotion(R,t,w_gt,d_gt,cameraParams);
camera_set_mid{i}.gt_oritation = Rot_Rows{param.cy};
camera_set_mid{i}.gt_translation = trans_Rows{param.cy};
camera_set_mid{i}.w_gt = w_gt * param.fy;
camera_set_mid{i}.d_gt = d_gt * param.fy;
end
camera_set_mid{1}.oritation = camera_set_mid{1}.gt_oritation;
camera_set_mid{1}.translation = camera_set_mid{1}.gt_translation;

end