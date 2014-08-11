function S = extendSkeleton(S, perimPixels)

if all(S(1:5,1)==S(1,1)) %horizontal
    while 1
        y = (S(1,2) - S(2,2)) + S(1,2);
        x = S(1,1);
        S=[[x,y];S];
        if perimPixels(x,y)
            break
        end
    end
else %not horizontal
    p = polyfit(S(1:5,1), S(1:5,2), 1);
    while 1
        x = (S(1,1) - S(2,1)) + S(1,1);
        y = round(p(1) * x + p(2));
        S=[[x,y];S];
        if perimPixels(x,y)
            break
        end
    end
end