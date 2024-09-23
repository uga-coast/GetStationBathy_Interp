function [z] = GetStationBathy_Interp(meshname,coordinate)

% Interpolates a coordinate's elevation from the nodes associated with the
% element it is located in

% Meshname:   File name of mesh, e.g. 'fort.14'. Must be a string.
% Coordinate: [longitude latitude] of coordinate to interpolate to

% REQUIRES adcirc_util: https://github.com/BrianOBlanton/adcirc_util


%%
clc;
grd = grd_to_opnml(meshname);

numElem = length(grd.e(:,1));
index = 0;

for j = 1:numElem
    if rem(j,1000) == 0
        clc; fprintf('%i%%\n',floor(100*j/numElem));
    else
    end
    
    elem = grd.e(j,:);
    triangle = [grd.x(elem(1)) grd.y(elem(1)); grd.x(elem(2)) grd.y(elem(2)); grd.x(elem(3)) grd.y(elem(3))];
    
    v0 = triangle(3, :) - triangle(1, :);
    v1 = triangle(2, :) - triangle(1, :);
    v2 = coordinate - triangle(1, :);
    
    dot00 = dot(v0, v0);
    dot01 = dot(v0, v1);
    dot02 = dot(v0, v2);
    dot11 = dot(v1, v1);
    dot12 = dot(v1, v2);
    
    invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    v = (dot00 * dot12 - dot01 * dot02) * invDenom;
    
    % Boolean indicating if the coordinate is inside the element
    if (u >= 0) && (v >= 0) && (u + v <= 1)
        index = j; clc;
        break
    else
    end
end

if index == 0
    error 'Coordinate is not within mesh domain'
else
end

elem_nodes = grd.e(index,:);

x1 = grd.x(elem_nodes(1));   y1 = grd.y(elem_nodes(1));   z1 = grd.z(elem_nodes(1));
x2 = grd.x(elem_nodes(2));   y2 = grd.y(elem_nodes(2));   z2 = grd.z(elem_nodes(2));
x3 = grd.x(elem_nodes(3));   y3 = grd.y(elem_nodes(3));   z3 = grd.z(elem_nodes(3));

gamma1 = grd.z(elem_nodes(1));
gamma2 = grd.z(elem_nodes(2));
gamma3 = grd.z(elem_nodes(3));

a1 = x3 - x2;   a2 = x1 - x3;   a3 = x2 - x1;
b1 = y2 - y3;   b2 = y3 - y1;   b3 = y1 - y2;

A = ((b1*a2) - (b2*a1))/2;

phi1  = @(x,y)((x2*y3) - (x3*y2) + (b1*x) + (a1*y))/(2*A);
phi2  = @(x,y)((x3*y1) - (x1*y3) + (b2*x) + (a2*y))/(2*A);
phi3  = @(x,y)((x1*y2) - (x2*y1) + (b3*x) + (a3*y))/(2*A);

gamma = @(x,y) gamma1*phi1(x,y) + gamma2*phi2(x,y) + gamma3*phi3(x,y);

z = gamma(coordinate(1),coordinate(2));
fprintf('Finished\n')
