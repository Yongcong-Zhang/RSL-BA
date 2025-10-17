function plucker = orth_to_plucker(orth)
angles = orth(1:3);
U = eul2rotm(angles, 'XYZ');
u1 = U(:,1);
u2 = U(:,2);
theta = orth(4);
w1 = cos(theta);
w2 = sin(theta);

L = [w1*u1', w2*u2']';

plucker = [[0 -L(3) L(2) L(4)]; [L(3) 0 -L(1) L(5)]; [-L(2) L(1) 0 L(6)]; [-L(4) -L(5) -L(6) 0]];
end