function [ model ] = getCellinfo( model )
R = model.sector_info.R;
a = sqrt(3)*R;
A = 2*R;
r = a/2;
theta = pi/3;
ISD = 3*R;
MaxRadius = ISD*(model.config.Num_tier-1)+A;
model.sector_info.R = R;
model.sector_info.A = A;
model.sector_info.a = a;
model.sector_info.r = r;
model.sector_info.theta = theta;
model.sector_info.ISD = ISD;
model.sector_info.MaxRadius = MaxRadius;
end

