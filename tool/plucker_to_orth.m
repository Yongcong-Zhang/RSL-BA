function orth = plucker_to_orth(plucker)
d = [plucker(1,4); plucker(2,4); plucker(3,4)];
n = [-plucker(2,3);  plucker(1,3); -plucker(1,2)];

U = [n./norm(n), d./norm(d), X_(n)*d./(norm(X_(n)*d))];
angles = rotm2eul(U, 'XYZ');
W = (1/(norm(n)^2+norm(d)^2)^0.5).*[[norm(n), -norm(d)]; [norm(d), norm(n)]];
theta = acos(W(1));

orth = [angles(1), angles(2), angles(3), theta];
end