% Use sphere fitting to calibrate pointer 8700340
warning('off','all')

stylusID = '8700340';
markerID = '8700339';

%From James Stewart provided for assignment
dataFile = 'pivot_calibration_0.csv';
%Collected by KB
%dataFile = 'pivot_calibration_1.csv';
%Colltected by GP
%dataFile = 'pivot_calibration_2.csv';

setupDrawing();

% Read the raw data into 'pos' (a translation 3-vector) and 'orient'
% (a quaternion 4-vector).

[pos, orient] = read_NDI_data( dataFile, stylusID );

% %TEST
% pos = pos(200:size(pos,1),:);
% orient = orient(200:size(orient,1),:);

% fit the sphere to get centre c and radius r

[c, r] = fitSphere( pos );

% % RANSAC outlier removal
% [c, r, bestInliers] = fitSphereWithRANSAC(pos);
% pos_best = zeros(length(bestInliers), 3);
% orient_best = zeros(length(bestInliers), 4);
% for ix=1:length(bestInliers)
%     pos_best(ix,:) = pos(bestInliers(ix),:);
%     orient_best(ix,:) = orient(bestInliers(ix),:);
% end
% pos = pos_best;
% orient = orient_best;

% Show the fit
drawCoordSystems( pos, orient );
drawSphere( c, r );

% Transform c into the coordinate system of each pose
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

%drawPointsWithEllipsoid( c_world, c_world_stdv );

drawPointsWithEllipsoidWithArrow(c_world, c_world_stdv );

% ---------------- END OF MAIN CODE ----------------


% Fit a sphere to a set of positions
%
% See http://watkins.cs.queensu.ca/~jstewart/472/notes/08-matrix-applications/08-matrix-applications.html

function [c, r] = fitSphere( pos )

    %Use least squares solutions to solve Ax=b
    
    Ax = [ mean(pos(:,1).*(pos(:,1)-mean(pos(:,1)))) 2*mean(pos(:,1).*(pos(:,2)-mean(pos(:,2))))  2*mean(pos(:,1).*(pos(:,3)-mean(pos(:,3))))];
    Ay = [ 0  mean(pos(:,2).*(pos(:,2)-mean(pos(:,2)))) 2*mean(pos(:,2).*(pos(:,3)-mean(pos(:,3))))];
    Az = [ 0 0 mean(pos(:,3).*(pos(:,3)-mean(pos(:,3)))) ];
    A = [Ax; Ay; Az];
    A = A + A.';
    
    Bx = mean((pos(:,1).^2+pos(:,2).^2+pos(:,3).^2).*(pos(:,1)-mean(pos(:,1))));
    By = mean((pos(:,1).^2+pos(:,2).^2+pos(:,3).^2).*(pos(:,2)-mean(pos(:,2))));
    Bz = mean((pos(:,1).^2+pos(:,2).^2+pos(:,3).^2).*(pos(:,3)-mean(pos(:,3))));
    B=[Bx; By; Bz];
    
    c = (A\B).';
    r = sqrt(mean(sum([pos(:,1)-c(1), pos(:,2)-c(2), pos(:,3)-c(3)].^2,2)));
end
  

% Fit a sphere to a set of positions using RANSAC.
%
% ALSO RETURN THE INDICES OF THE BEST INLIERS.  THE CALLING CODE
% SHOULD RESTRICT ITSELF TO THOSE INLIERS.
%
% See https://en.wikipedia.org/wiki/Random_sample_consensus

function [c, r, bestInlierIndices] = fitSphereWithRANSAC( pos )
    
    iters = 0;
    num_points = size(pos,1);
    
    %Number of points required to use sphere aka "ENOUGH"
    min_points_needed = .9*num_points;
    
    %Attempt 500 iterations of RANSAC
    while iters < 500
        
        %randomly pick 4 points
        rand_indices = randsample(size(pos,1), 4);
        rand_points = [pos(rand_indices(1),:); pos(rand_indices(2),:); pos(rand_indices(3),:); pos(rand_indices(4),:)];

        %calculate sphere defined by 4 points
        [c_temp, r_temp] = fitSphere(rand_points);
        
        %Distance points must be from sphere to count as "CLOSE"
        max_dist_from_sphere = .05*r_temp;

        %find distance (abs val) between all points and sphere
        distances = zeros(num_points,1);
        for i = 1:num_points
            distances(i) = abs(norm(pos(i,:) - c_temp) - r_temp);
        end
        
        %If ENOUGH points (% of points defined by min_points_needed) 
        %are CLOSE (within range defined by max_dist_from_sphere) stop
        
        %Get indices of points which are within max distance from sphere
        currentInlierIndices=find(distances<max_dist_from_sphere);
        bestInlierIndices = [];
        
        %If count of indices is greater than set % of points return sphere
        if length(currentInlierIndices) > min_points_needed
            bestInlierIndices = currentInlierIndices;
            c = c_temp;
            r = r_temp;
            disp("FOUND BEST");
            break
        end
        
        %Keep track of best sphere (if not returning sphere)
        if size(currentInlierIndices,1) > size(bestInlierIndices,1)
            bestInlierIndices = currentInlierIndices;
            c = c_temp;
            r = r_temp;
        end
                
        iters = iters + 1;
    end
    
end

% From https://www.mathworks.com/matlabcentral/fileexchange/35475-quaternions
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

% Draw all points

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

function drawPointsWithEllipsoid( points, stdev )

    %Find mean aka center of ellipsoid
    mean_point = mean(points,1);
    
    %Define constant for use with 95% confidence interval
    std_dev_95_constant = 1.960;

    %Calculate the radius in each direction
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

function drawPointsWithEllipsoidWithArrow(points, stdev)
    [X,Y,Z] = ellipsoid(0,0,0,1,1,1);

    [eigvec,eigval] = eig(cov(points));

    XYZ = [X(:),Y(:),Z(:)] * (1.96*sqrt(eigval)) * eigvec';

    mu = mean( points );

    X(:) = XYZ(:,1)+mu(1);
    Y(:) = XYZ(:,2)+mu(2);
    Z(:) = XYZ(:,3)+mu(3);

    figure; hold on;
    xlabel("X");
    ylabel("Y");
    zlabel("Z");
    scatter3(points(:,1),points(:,2),points(:,3));
    scatter3(mu(1), mu(2), mu(3), 50, 'red');
    surf( X, Y, Z, 'FaceAlpha', 0.1);
    axis equal;

    mn = min(points);
    mx = max(points);
    axis( [mn(1),mx(1),mn(2),mx(2),mn(3),mx(3)] );

    %This draws the vector along the stylus axis:
    
    tip = mean( points); 
    arrow3( [0, 0, 0], [tip(1), tip(2), tip(3)], 'f', 1 );

end

