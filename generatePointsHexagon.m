function points = generatePointsHexagon(dimensions,radius,minDistance)

whichRhombus = randi(3,dimensions);

xDim = radius*rand(dimensions);
yDim = radius*rand(dimensions);


points = -1i*xDim + (sqrt(3)/2+1i/2)*yDim;

points = points .* exp(1i*whichRhombus*2*pi/3);


if nargin>2

    notReady = abs(points) < minDistance;

    if ~isempty(notReady)
    
        points(notReady) = generatePointsHexagon([sum(notReady(:)) 1],radius,minDistance);
    
    end

end
