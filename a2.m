% Use sphere fitting to calibrate pointer 8700340
warning('off','all')

stylusID = '8700340';
markerID = '8700339';
dataFile = 'pivot_calibration_0.csv';

setupDrawing();

% Read the raw data into 'pos' (a translation 3-vector) and 'orient'
% (a quaternion 4-vector).

[pos, orient] = read_NDI_data( dataFile, stylusID );

pos = pos(200:size(pos,1),:);
orient = orient(200:size(orient,1),:);

% fit the sphere to get centre c and radius r

[c, r] = fitSphere( pos );

[c_ransac,r_ransac, bestHoes] = fitSphereWithRANSAC(pos);

% Show the fit
%drawPoints( pos );
drawCoordSystems( pos, orient );
drawSphere( c, r );

% Transform c into the coordinate system of each pose

% [YOUR CODE HERE]
allCenters = zeros(length(pos),3);
for i=1:size(pos,1)

    m = quaternion_to_matrix( orient(i,:) );
    t = pos(i,:);
    
    allCenters(i,:) = m.' * (c-t).';
end


% Find the average transformed c, which should be the same in all of
% the stylus coordinate systems.  Also find the standard deviation.

c_average = mean(allCenters, 1);
c_stdev = std(allCenters,1,1);

% Report the results
%
% 'c_average' is the average tip position in the stylus coordinate
% system.  'c_stdev' is its standard deviation in the stylus
% coordinate system.

disp( sprintf( 'tip position in stylus CS: (%g, %g, %g)', c_average(1), c_average(2), c_average(3) ) );
disp( sprintf( 'tip stdev in stylus CS:    (%g, %g, %g)', c_stdev(1),   c_stdev(2),   c_stdev(3) ) );

% Show vectors to tips in world coordinate system
%
% This is for debugging, so that you can see that the vector touches
% the same pivot point from all stylus coordinate systems.

drawLocalVectorInCoordSystems( c_average, pos, orient );

% Show tip points in global system, along with 95% confidence
% interval as an ellipsoid.
%
% 'c_world' are the tip points in the world coordinate system.
% They should all be very near the pivot point.

c_world = zeros(length(pos),3);
for i=1:size(pos,1)

    m = quaternion_to_matrix( orient(i,:) );
    t = pos(i,:);
    
    c_world(i,:) = m*c_average.' + t.';
end
c_world_average = mean(c_world, 1);
c_world_stdv = std(c_world, 1);



drawPointsWithEllipsoid( c_world, c_stdev );



% ---------------- END OF MAIN CODE ----------------



% Fit a sphere to a set of positions
%
% See http://watkins.cs.queensu.ca/~jstewart/472/notes/08-matrix-applications/08-matrix-applications.html

function [c, r] = fitSphere(X)
    
    A=[mean(X(:,1).*(X(:,1)-mean(X(:,1)))), ...
    2*mean(X(:,1).*(X(:,2)-mean(X(:,2)))), ...
    2*mean(X(:,1).*(X(:,3)-mean(X(:,3)))); ...
    0, ...
    mean(X(:,2).*(X(:,2)-mean(X(:,2)))), ...
    2*mean(X(:,2).*(X(:,3)-mean(X(:,3)))); ...
    0, ...
    0, ...
    mean(X(:,3).*(X(:,3)-mean(X(:,3))))];
    A=A+A.';
    B=[mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,1)-mean(X(:,1))));...
        mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,2)-mean(X(:,2))));...
        mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,3)-mean(X(:,3))))];
    c=(A\B).';
    r=sqrt(mean(sum([X(:,1)-c(1),X(:,2)-c(2),X(:,3)-c(3)].^2,2)));
end
  

% Fit a sphere to a set of positions using RANSAC.
%
% ALSO RETURN THE INDICES OF THE BEST INLIERS.  THE CALLING CODE
% SHOULD RESTRICT ITSELF TO THOSE INLIERS.
%
% See https://en.wikipedia.org/wiki/Random_sample_consensus

function [c, r, bestInlierIndices] = fitSphereWithRANSAC( pos )
    
    %[ YOUR CODE HERE ]
    num_points = size(pos,1);
    max_dist_from_sphere = 5;
    min_points_needed = .1*num_points;
    iters = 0;
    
    while iters < 500
        %randomly pick 4 points
        rand_indices = randsample(size(pos,1), 4);
        rand_points = [pos(rand_indices(1),:); pos(rand_indices(2),:); pos(rand_indices(3),:); pos(rand_indices(4),:)];

        %calculate sphere defined by 4 points
        [c, r] = calculateSphereFromPoints(rand_points);

        %find distance between all points and sphere
        distances = zeros(num_points,1);
        for i = 1:num_points
            distances(i) = abs(norm(pos(i,:) - c) - r);
        end
        
        %If ENOUGH points (% of points defined by min_points_needed) 
        %are CLOSE (within range defined by max_dist_from_sphere) stop
        
        %Get indices of points which are within max distance from sphere
        bestInlierIndices=find(distances<max_dist_from_sphere);
        
        %If count of indices is greater than set % of points
        if size(bestInlierIndices,1) > min_points_needed
            break
        end
        iters = iters + 1;
    end
    
end

function [c, r] = calculateSphereFromPoints( pts )
    p1 = pts(1,:);
    p2 = pts(2,:);
    p3 = pts(4,:);
    p4 = pts(4,:);
    
    sympref('FloatingPointOutput',true);
    syms a b c r;
    eqn1 = (p1(1) - a).^2 + (p1(2) - b).^2 + (p1(3) - c).^2 == r.^2;
    eqn2 = (p2(1) - a).^2 + (p2(2) - b).^2 + (p2(3) - c).^2 == r.^2;
    eqn3 = (p3(1) - a).^2 + (p3(2) - b).^2 + (p3(3) - c).^2 == r.^2;
    eqn4 = (p4(1) - a).^2 + (p4(2) - b).^2 + (p4(3) - c).^2 == r.^2;
    
    eqn1 = eqn1 - eqn4;
    eqn2 = eqn2 - eqn4;
    eqn3 = eqn3 - eqn4;
    
    
    [A, B] = equationsToMatrix([eqn1, eqn2, eqn3], [a,b,c]);
    
    Unknowns = linsolve(A,B);
    
    rSquared = (p4(1) - Unknowns(1)).^2 + (p4(2) - Unknowns(2)).^2 + (p4(3) - Unknowns(3)).^2;
    
    r = double(sqrt(rSquared));
    
    center_x = double(Unknowns(1));
    center_y = double(Unknowns(2));
    center_z = double(Unknowns(3));
    c = [center_x, center_y, center_z];
    
    % TESTING5
%     figure; hold on;
%     scatter3(pts(:,1), pts(:,2), pts(:,3), 'blue');
%     scatter3(center_x, center_y, center_z, 'red');
%     drawSphere([center_x, center_y, center_z], r);
%     axis equal;
%     figure;

end
  

% From https://www.mathworks.com/matlabcentral/fileexchange/35475-quaternions
%
% Convert quaternion to 3x3 matrix.

function R = quaternion_to_matrix( Qrotation )

  w = Qrotation( 1 );
  x = Qrotation( 2 );
  y = Qrotation( 3 );
  z = Qrotation( 4 );

  Rxx = 1 - 2*(y^2 + z^2);
  Rxy = 2*(x*y - z*w);
  Rxz = 2*(x*z + y*w);
  Ryx = 2*(x*y + z*w);
  Ryy = 1 - 2*(x^2 + z^2);
  Ryz = 2*(y*z - x*w );
  Rzx = 2*(x*z - y*w );
  Rzy = 2*(y*z + x*w );
  Rzz = 1 - 2 *(x^2 + y^2);

  R = [ Rxx,    Rxy,    Rxz;
        Ryx,    Ryy,    Ryz;
        Rzx,    Rzy,    Rzz  ];
end


% Read a CSV file from the NDI tracker software and extract the
% position and orientation columns of the named marker.
%
% pos    = n x 3 of translations
% orient = n x 4 of quaternion orientations

function [pos, quat] = read_NDI_data( dataFile, markerID )

  t = readtable( dataFile, 'PreserveVariableNames', true );

  % Find the column 

  colIndex = find(contains( t.Properties.VariableNames, markerID ));

  if colIndex == []
    disp( sprintf( "In %s: Could not find a column header containing '%s'.", dataFile, markerID ) );
    exit;
  end

  % From the rows with state == 'OK' (state is at +3 offset from the ID
  % column), extract the Tx, Ty, Tz columns.

  status = t{:,colIndex+3};
  n = size( t( strcmp(status,'OK'), 1 ), 1 );

  pos = zeros( n, 3 );
  quat = zeros( n, 4 );
  k = 1;

  for i = 1:size(t,1)
    if strcmp( t{i,colIndex+3}, 'OK' )

      % Coerce the columns to 'double', since 'readtable' sometimes
      % records numbers as strings.  MATLAB BUG!

      % Extract pose's rotation as a quaternion.  This is in columns
      % offset by +4, +5, +6, +7 from the ID column.

      for j=1:4
	if iscell( t{i,colIndex+3+j} )
	  quat(k,j) = str2double( t{i,colIndex+3+j}{1} );
	else
	  quat(k,j) = t{i,colIndex+3+j};
	end
      end

      % Extract pose's translation as a vector.  This is in columns
      % offset by +8, +9, +10 from the ID column.

      for j=1:3
	if iscell( t{i,colIndex+7+j} )
	  pos(k,j) = str2double( t{i,colIndex+7+j}{1} );
	else
	  pos(k,j) = t{i,colIndex+7+j};
	end
      end

      k = k + 1;
    end
  end
  
  disp( sprintf( '%d points collected', size(pos,1) ) );
end



% Set up the drawing

function setupDrawing() 

    f = figure(1);
    clf(f);
    view(3);
    daspect( [1 1 1] );
    pbaspect manual;
    hold on;
end


% Draw a set of coordinate systems
%
% pos and orient store their vectors in rows.

function drawCoordSystems( pos, orient )
    
    colours = [ 'r' 'g' 'b' ];

    scale = 0.005 * norm(max(pos) - min(pos));
    for i=1:size(pos,1)
        m = quaternion_to_matrix( orient(i,:) );
        t = pos(i,:);
        for j=1:3
            head = t + m(:,j)' .* scale;
            x = [t(1) head(1)];
            plot3( [t(1) head(1)], [t(2) head(2)], [t(3) head(3)], colours(j) );
        end
    end
end

function drawPoints( pos )
    %figure;hold on;
    scatter3(pos(:,1),pos(:,2),pos(:,3),50,'filled');
end



% Draw a sphere

function drawSphere( c, r )
    
    [x,y,z] = sphere;
    surf( x*r+c(1), y*r+c(2), z*r+c(3), 'FaceAlpha', 0.05, 'FaceColor', [0.6 0.3 0.3] );
end


% Draw a local vectors in different coordinate systems
%
% v is a row vector.  pos and orient store their vectors in rows, too.

function drawLocalVectorInCoordSystems( v, pos, orient )
    
    for i=1:size(pos,1)
        m = quaternion_to_matrix( orient(i,:) );
        t = pos(i,:);
        for j=1:3
            head = t + (m * v')';
            plot3( [t(1) head(1)], [t(2) head(2)], [t(3) head(3)] );
        end
    end
end




% Draw points with a 95% CI ellipsoid.
%
% Use matlab's 'ellipsoid' and 'surf' functions.
%
% [YOUR CODE HERE]

function drawPointsWithEllipsoid( points, stdev )

    mean_point = mean(points,1);
    std_dev_95_constant = 1.960;
    num_points = size(points,1);
    step_maybe = stdev(1) / sqrt(num_points);
    x_radius = std_dev_95_constant * stdev(1);
    y_radius = std_dev_95_constant * stdev(2);
    z_radius = std_dev_95_constant * stdev(3)
    
    [x,y,z] = ellipsoid(mean_point(1), mean_point(2), mean_point(3), x_radius, y_radius, z_radius);
    
    figure; hold on;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    scatter3(points(:,1),points(:,2),points(:,3));
    scatter3(mean_point(1), mean_point(2), mean_point(3), 50, 'red');
    surf(x,y,z, 'FaceAlpha', 0.05, 'FaceColor', [0.6 0.3 0.3])
    axis equal;
end


