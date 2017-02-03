%This Matlab function is used in the simulations of the following article:
%
%Emil Bjornson, Elisabeth de Carvalho, Jesper H. Sorensen, Erik G. Larsson,
%Petar Popovski, "A Random Access Protocol for Pilot Allocation in Crowded
%Massive MIMO Systems," IEEE Transactions on Wireless Communications,
%To appear.
%
%Download article: http://arxiv.org/pdf/1604.04248
%
%This is version 1.0 (Last edited: 2017-02-03)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
%
%This function selects "nbrOfPoints" uniformly at random in a hexagonal
%cell centered in the origin. The radius of the hexagon is "radius" and the
%(optional) minimum distance from the origin is "minDistance".


function points = generatePointsHexagon(nbrOfPoints,radius,minDistance)

%The hexagon is divided into three rhombus. Each point is uniformly
%distributed among these rhombus
whichRhombus = randi(3,nbrOfPoints);

%Place the points uniformly distributed in a square
xDim = radius*rand(nbrOfPoints);
yDim = radius*rand(nbrOfPoints);

%Rotate the points in the squares to make them uniformly distributed in
%the right rhombus
points = -1i*xDim + (sqrt(3)/2+1i/2)*yDim;
points = points .* exp(1i*whichRhombus*2*pi/3);

%Remove points that are closer to the origin than "minDistance" and
%generate new points in the hexagon to replace these ones.
if nargin>2

    notReady = abs(points) < minDistance;

    if ~isempty(notReady)
    
        points(notReady) = generatePointsHexagon([sum(notReady(:)) 1],radius,minDistance);
    
    end

end
