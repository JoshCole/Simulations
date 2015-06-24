function [ model ] = getAlpaha( model )
j = model.config.j;
if (j == 0) || (j == 1)
    alpha = [1];
%     alpha = [0 0.4 0.5 0.6 0.7 0.8 0.9 1];
elseif j == 2
    alpha = 1;
end
model.ue.alpha = alpha;

end

