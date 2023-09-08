function R = circumscribedR(X,Y)
% Calculate radius of circumscribed circle through three point in vectors X
% and Y of size 1*3

A = [X(1) Y(1)]; B = [X(2) Y(2)]; C = [X(3) Y(3)];

a = norm(B-C);
b = norm(C-A);
c = norm(A-B);

s = (a+b+c)/2;
temp = sqrt(s*(s-a)*(s-b)*(s-c));

R = a*b*c/4/temp;
end