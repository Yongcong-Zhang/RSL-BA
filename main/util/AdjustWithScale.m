function [camera_param_adjusted, scale] = AdjustWithScale(camera_param, gt_camera_param, param)
    camera_param_adjusted = camera_param;
    num_camera = param.num_camera;
    num_point = param.num_point;
    scale_gt = 0;
    scale = 0;
    gt_camera_t_fir = gt_camera_param(1, 4:6);
    gt_camera_r_fir = gt_camera_param(1, 1:3);
    gt_camera_r_fir = axang2rotm([gt_camera_r_fir / norm(gt_camera_r_fir) norm(gt_camera_r_fir)]);
    for i = 2:num_camera
        id = (i - 1) * 12 + 1;
        gt_camera_rota = gt_camera_param(1, id:id + 2);
        gt_camera_r = axang2rotm([gt_camera_rota / norm(gt_camera_rota) norm(gt_camera_rota)]);
        gt_camera_tran = gt_camera_param(1, id + 3:id + 5);
        gt_camera_t = gt_camera_r' * gt_camera_tran' - gt_camera_r_fir' * gt_camera_t_fir';
        scale_gt = scale_gt + norm(gt_camera_t);
    end
    for i = 2:num_camera
        id = (i - 1) * 12 + 1;
        camera_rota = camera_param(1, id:id + 2);
        camera_r = axang2rotm([camera_rota / norm(camera_rota) norm(camera_rota)]);
        camera_tran = camera_param(1, id + 3:id + 5);
        camera_t = camera_r' * camera_tran' - gt_camera_r_fir' * gt_camera_t_fir';
        scale = scale + norm(camera_t);
    end
    scale = scale_gt / scale;
    for i = 2:num_camera
        id = (i - 1) * 12 + 1;
        camera_rota = camera_param(1, id:id + 2);
        camera_r = axang2rotm([camera_rota / norm(camera_rota) norm(camera_rota)]);
        camera_tran = camera_param(1, id + 3:id + 5);
        camera_t = camera_r' * camera_tran' - gt_camera_r_fir' * gt_camera_t_fir';
        camera_t = camera_t * scale;
        camera_tran = camera_r * (camera_t + gt_camera_r_fir' * gt_camera_t_fir');
        camera_param_adjusted(1, id + 3:id + 5) = camera_tran';
        camera_param_adjusted(1, id + 9:id + 11) = camera_param_adjusted(1, id + 9:id + 11) * scale;
    end
    for i = 1:num_point
        id = num_camera * 12 + (i - 1) * 3 + 1;
        point = camera_param(1, id:id + 2);
        
        point = point + (gt_camera_r_fir' * gt_camera_t_fir')';
        camera_param_adjusted(1, id:id + 2) = point * scale - (gt_camera_r_fir' * gt_camera_t_fir')';
    end

end