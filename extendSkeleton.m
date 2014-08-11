function S = extendSkeleton(S, perimPixels)
% extend skeleton to end of worm by linear extrapolation using the last 5 points
% after each new point is added, a new window of the 5 most last points is used
S = extend(S, perimPixels);
S = flipud(S);
S = extend(S, perimPixels);
S = flipud(S);

function S = extend(S, perimPixels)

window = 5;

while 1
    lastPoint = [S(1,1),S(1,2)];
    p = polyfit(S(1:window,1), S(1:window,2), 1);
    if all(S(1:window,1)==S(1,1)) %horizontal
        y = round(S(1,2) - mean(diff(S(1:window,2))));
        x = (y - p(2)) / p(1);
    else
        x = round(S(1,1) - mean(diff(S(1:window,1))));
        y = round(p(1) * x + p(2));
    end
    if all(lastPoint == [x,y])
        x=x-1;
        y = round(p(1) * x + p(2));
        warning('Check the skeleton extension, may not be good')
    end
    S=[[x,y];S];
    if perimPixels(x,y)
        break
    end
end