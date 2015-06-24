function [ model ] = BSradiation_pattern( model )
% Base station antenna radiation pattern
switch model.config.environment
    case 'Urban'
     model.bs.Dhb  = 30;
     model.bs.antennaGain = 15 ; %2000MHz
    case 'Rural'
       Dhb  = 30;
end
    
theta(:,1) = -180:180;
Ntheta = length(theta);
g_rx = zeros(Ntheta,1);
Am = 20; %dB
theta_3dB = 65;
for i = 1:Ntheta
    g_rx(i)= -min(12*(theta(i)/theta_3dB).^2,Am);
end
antenna3 = [g_rx(181:end);g_rx(1:180)];
antenna1 = [antenna3(242:end);antenna3(1:241)];
antenna2 = [antenna3(122:end);antenna3(1:121)];
lin1 = 10.^(antenna1/10);
lin2 = 10.^(antenna2/10);
lin3 = 10.^(antenna3/10);
newTheta(:,1) = (0:360)*pi/180;
figure
polar(newTheta,lin1)
hold on
polar(newTheta,lin2)
polar(newTheta,lin3)

model.bs.data.radiation_pattern = [antenna1.'; antenna2.'; antenna3.';]+model.bs.antennaGain;
end

