%% Output the statistic of the stl file.
clear all; close all; clc;

% Reading STL using stlread
GEO  = stlread('../stl/Perdigao_rotated_clipped_8m_repaired.stl');
% sphere = = stlread('./stl/sphere_gmshv2.stl');
% perdigao = stlread('./stl/dragon_allFix.stl');
% bunny = = stlread('./stl/bunny_allFix.stl');

patch(GEO,'FaceColor', [0.8 0.8 1.0], ...
    'FaceLighting', 'gouraud', ...
    'AmbientStrength', 0.15);

% Fix the axes scaling, and set a nice view angle
axis('image');
view([-135 35]);
grid on;
movegui('northeast')

% triangle is made up of 3 verticies, and each verticies have an X,Y,Z
% coordinate
face = numel(GEO.faces)/3;
vertex = numel(GEO.vertices)/3;
length = zeros(face,3);

j = 1;
for i=1:3:vertex
    v1 = GEO.vertices(i,:);
    v2 = GEO.vertices(i+1,:);
    v3 = GEO.vertices(i+2,:);
    length(j,1) = sqrt((v2(1)-v1(1))^2+((v2(2)-v1(2)))^2+(v2(3)-v1(3))^2);
    length(j,2) = sqrt((v3(1)-v2(1))^2+((v3(2)-v2(2)))^2+(v3(3)-v2(3))^2);
    length(j,3) = sqrt((v1(1)-v3(1))^2+((v1(2)-v3(2)))^2+(v1(3)-v3(3))^2);
    j = j + 1;
end

%finding maximum
Lmax = zeros(1,3);
Lmin = zeros(1,3);
for i = 1:3
    Lmax(i) = max(length(i,:)); %X,Y,Z max
    Lmin(i) = min(length(i,:)); %X,Y,Z min
end
globalMax = max(Lmax);
globalMin = min(Lmin);

% the average length in x
Avg = zeros(1,3);
for i = 1:3
    Avg(i) = sum(length(:,i))/(face);
end
fprintf("**********STL Info**********\nTriangle=%d, Vertex=%d\n",face, vertex)
fprintf("Average[x y z] = [%f %f %f]\nglobal[min max] = [%f %f]\n",globalMin,globalMax, Avg(1), Avg(2), Avg(3))
%% Carteisan mesh quality
Cart = [16000 8000 3300];
CartMesh = [129 65 65]; %449 257 257 was previous
delta = [Cart(1)/(CartMesh(1)-2), Cart(2)/(CartMesh(2)-2), Cart(3)/(CartMesh(3)-2)];

Diag = sqrt(delta(1).^2+delta(2).^2+delta(3).^2);
ratioMin = Diag / globalMax;
ratioMax = Diag / globalMin;

for i = 1:3
    ratiodirection(i) = delta(i) / Avg(i);
end

fprintf("**********Cartesian Info**********\ndiagonal = %f\n", Diag)
fprintf("[dx dy dz] = [%f %f %f]\n", delta(1),delta(2),delta(3))
fprintf("ratio[min max] =  [%f %f]\nratioAvg[x y z] = [%f %f %f]\n", ratioMin, ratioMax, ratiodirection(1),ratiodirection(2),ratiodirection(3))
