function [ eyeD ] = eyeDiagram( val , x )

y=x(val+(1:640));
y = reshape(y,16,40);
eyeD = real(y);
end

