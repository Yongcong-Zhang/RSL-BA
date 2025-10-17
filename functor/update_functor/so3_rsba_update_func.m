function x = so3_rsba_update_func(x, dx, param)
    for i = 1:param.num_camera
        r_id = 12 * (i - 1) + 1;
        t_id = 12 * (i - 1) + 4;
        w_id = 12 * (i - 1) + 7;
        d_id = 12 * (i - 1) + 10;

            r_x_vector = x(1, r_id:r_id + 2);
            r_dx_vector = dx(1, r_id:r_id + 2);
            if norm(r_x_vector) <= 1e-10
                r_x = eye(3);
            else
                r_x = axang2rotm([r_x_vector / norm(r_x_vector) norm(r_x_vector)]);
            end

            if norm(r_dx_vector) <= 1e-10
                r_dx = eye(3);
            else
                r_dx = axang2rotm([r_dx_vector / norm(r_dx_vector) norm(r_dx_vector)]);
            end

            
            
            axang = rotm2axang(r_dx * r_x);
            axis_vec = axang(1:3) * axang(4);
            x(1, r_id:r_id + 2) = axis_vec;
            x(1, t_id:t_id + 2) = x(1, t_id:t_id + 2) + dx(1, t_id:t_id + 2);

        x(1, w_id:w_id + 2) = x(1, w_id:w_id + 2) + dx(1, w_id:w_id + 2);
        x(1, d_id:d_id + 2) = x(1, d_id:d_id + 2) + dx(1, d_id:d_id + 2);
    end
    x(1, 12 * param.num_camera + 1:end) = x(1, 12 * param.num_camera + 1:end) + dx(1, 12 * param.num_camera + 1:end);
end