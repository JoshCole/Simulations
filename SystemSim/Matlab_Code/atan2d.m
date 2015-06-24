function f = atan2d(y,x)
% Function calculates the arc tangent of y/x and placees the result in
% range of [0..360]
f = 0;
if x == 0
    if y == 0
        elseif y > 0
           f= 90;
        else
   f = 270;
    end
elseif x > 0
            if y >= 0
                f = atand(y/x);
            else
                    f = atand(y/x) + 360;
            end
elseif x < 0
if y == 0
f = 180;
else
f = atand(y/x) + 180;
end
end
end