function draw_line(a1,a2,a3,a4,a5,a6,a7,v1,v2)
%figure;
if v1 < v2
    y_values = v1:0.1:v2;
else
    y_values = v2:0.1:v1;
end
x_values = y_values;

for i = 1:length(y_values)
    v = y_values(i);
    x_values(i) = -(a1*v*v*v+a3*v*v+a5*v+a7)/(a2*v*v+a4*v+a6);
end


plot(x_values, y_values);
end