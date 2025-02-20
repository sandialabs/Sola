function points = sampleUnitCircle(numPoints)
    angles = linspace(0, 2 * pi, numPoints + 1);
    angles(end) = [];
    x = cos(angles);
    y = sin(angles);
    points = [x', y'];
end
