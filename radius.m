function radius = radius(p1,p2,p3)

% Calculate midpoint of p1-p2 and its slope and y-intercept
x_mid = (p1(1) + p2(1))/2;
y_mid = (p1(2) + p2(2))/2;
slope = -(p2(1) - p1(1)) / (p2(2) - p1(2));
y_intercept = y_mid - slope * x_mid;

% Calculate midpoint of p2-p3 and its slope and y-intercept
x_mid = (p2(1) + p3(1))/2;
y_mid = (p2(2) + p3(2))/2;
slope = -(p3(1) - p2(1)) / (p3(2) - p2(2));
y_intercept = y_mid - slope * x_mid;

% Calculate midpoint of p3-p1 and its slope and y-intercept
x_mid = (p3(1) + p1(1))/2;
y_mid = (p3(2) + p1(2))/2;
slope = -(p1(1) - p3(1)) / (p1(2) - p3(2));
y_intercept = y_mid - slope * x_mid;

% Find the point of intersection of the three perpendicular bisectors
x_center = (y_intercept(2) - y_intercept(1)) / (slope(1) - slope(2));
y_center = slope(1) * x_center + y_intercept(1);
center = [x_center, y_center];

% Use the distance formula to find the radius
radius = sqrt((center(1) - p1(1))^2 + (center(2) - p1(2))^2);
