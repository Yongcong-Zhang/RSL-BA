function line_2d = get_project_line(K, L, R, T, w, d, v1, v2)
    P0 = K*[R, T];
    Q = K*[X_(w)*R, d];
    A1 = P0*L*P0';
    A2 = P0*L*Q' + Q*L*P0';
    A3 = Q*L*Q';
    line_2d = [A3(1,3), A3(3,2), A2(1,3)+A3(2,1), A2(3,2), A1(1,3)+A2(2,1), A1(3,2), A1(2,1), v1, v2];
end