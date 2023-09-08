function Y = rotateCw(X,theta)
% rotate clockwise
rotMat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
Y = rotMat * X;
end