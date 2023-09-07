function [xOffset,yOffset] = offsetIn(x,y,R)
% Generate 2 vectors holding coordinates of offset curve by R from a closed
% curve.
% If R is positive, the offset curve will be inside the original curve.
% If R is negative, the offset curve will be outside the original curve.

L = length(x); % the length of both vectors
% Check whether the input is closed
%if (x(1) ~= x(L)) %the same point could be shifted because of calculation
%precision. Although the difference is extremely small, this condition may
%not hold true, hence the condition in the next line.
if (abs(x(1) - x(L)) > 10^-10)
    disp("The input is not a closed curve. The function will terminate.")
    return
end

xOffset = zeros(size(x));
yOffset = zeros(size(x));

% Boundary. Note that the first and the last points on pitch curve are the
% same
[xOffset(1),yOffset(1)] = normalIn([x(L-1) x(1) x(2)],[y(length(y)-1) y(1) y(2)],R);
[xOffset(L),yOffset(L)] = normalIn([x(L-1) x(L) x(2)],[y(L-1) y(L) y(2)],R);

for k = 1:1:L-2
    X = x(k:1:k+2);
    Y = y(k:1:k+2);
    [xOffset(k+1),yOffset(k+1)] = normalIn(X,Y,R);
end

end