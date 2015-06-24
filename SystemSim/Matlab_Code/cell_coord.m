function [ model ] = cell_coord( model )
% Cell parameters
R = model.sector_info.R;
ISD = 3*R;
% Inital cell site coord = x = 0, y = 0
model.bs.position = zeros(model.config.NumCells,2);

theta = pi/3*[0:5]';

theta_N = 0;
model.bs.position(2:7,:) = [ISD*cos(theta + ...
    theta_N) ISD*sin(theta + theta_N)];

model.bs.position(8:2:end,:) = [2*ISD*cos(theta + ...
    theta_N) 2*ISD*sin(theta + theta_N)];

theta_N = pi/6;
model.bs.position(9:2:end,:) = [sqrt(3)*ISD*cos(theta + ...
    theta_N) sqrt(3)*ISD*sin(theta + theta_N)];

end

