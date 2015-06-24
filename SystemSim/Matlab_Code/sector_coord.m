function [ model ] = sector_coord( model )
N = model.config.NumCells;
n = model.config.sector_size;
x = model.bs.position(:,1);
y = model.bs.position(:,2);
R = model.sector_info.R;
theta = model.sector_info.theta;
model.sector.position = zeros(3,2,N);
model.sector.coord = zeros(3*N,2);
c1 = 1; c2 = n;

for i = 1:N
    model.sector.position(:,:,i) = [x(i)-R/2 y(i)+R*sin(theta);x(i)-R/2 y(i)-R*sin(theta);x(i)+R y(i)+0];
    model.sector.coord(c1:c2,:) = [x(i)-R/2 y(i)+R*sin(theta);x(i)-R/2 y(i)-R*sin(theta);x(i)+R y(i)+0];
    c1 = c1+n; c2 = c2+n;
end

