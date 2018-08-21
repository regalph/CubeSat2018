% Robin Pflager
% June 2018
% UW CubeSat Orbit Propagator 

close all
clear vars
%% intial conditions and constants
G = 6.675452 * 10^-11; % N per kg^2 per m^2
au = 149597870700; % m
year = 31557600; % s, Julian Year
earth_radius = 6400000; % m (overestimate by ~30km)
sun_radius = 695700000; % m (overestimate by ~4300km)
m_sun = 1.98855 * 10^30; % kg
p_sun = [-449311, 0, 0]; % m
v_sun = [0, -2*pi*449311/ year, 0]; % m per s
m_earth = 5.97237 * 10^24; % kg
p_earth = [149597421389, 0, 0]; % m
v_earth = [0, 2*pi*149597421389/year, 0]; % m per s
a_0 = 5 * 10^-5; % m per s

%% define time resolution and length
TR = 1/480; % samples per second
L = year*4; % seconds per simulation

%% new object exiting Moon SOI
days_since_full = 0 ; % 0 corresponds to collinear with Earth and Sun, outside orbit of Earth
theta_time_of_month = days_since_full*(2*pi/29.5);
p_object = p_earth + 2847176085*[cos(theta_time_of_month),sin(theta_time_of_month),0];
v_object = v_earth + 1415.34*[cos(theta_time_of_month+1.23884),sin(theta_time_of_month+1.23884),0];

%% store
sun_path = zeros(TR*L, 3);
earth_path = zeros(TR*L, 3);
object_path = zeros(TR*L, 3);

sun_path_rot = zeros(TR*L, 3);
earth_path_rot = zeros(TR*L, 3);
object_path_rot = zeros(TR*L, 3);

theta_store = zeros(TR*L,1);
L4_p_difference_store = zeros(TR*L,3);
L4_p_difference_magnitude_store = zeros(TR*L,1);
L4_v_difference_store = zeros(TR*L,3);
L4_v_difference_magnitude_store = zeros(TR*L,1);
L5_p_difference_store = zeros(TR*L,3);
L5_p_difference_magnitude_store = zeros(TR*L, 1);
L5_v_difference_store = zeros(TR*L,3);
L5_v_difference_magnitude_store = zeros(TR*L, 1);
%% define thrust profile
thrust_profile = zeros(TR*L, 4);

% initial burn
burnlength = year; % s
for step = ((0)/135)*burnlength*TR+1:((10)/135)*burnlength*TR    
    thrust_profile(step, 1) = 0; % anti-velocity (forward)
    thrust_profile(step, 2) = sqrt(2)/2; % velocity (backward) 
    thrust_profile(step, 3) = sqrt(2)/2; % left velocity (rightward)
    thrust_profile(step, 4) = 0; % right velocity (leftward)
end
% for step = (80/135)*burnlength*TR:(92/135)*burnlength*TR    
%     thrust_profile(step, 1) = sqrt(2)/2; % anti-velocity (forward)
%     thrust_profile(step, 2) = 0; % velocity (backward)
%     thrust_profile(step, 3) = 0; % left velocity (rightward)
%     thrust_profile(step, 4) = sqrt(2)/2; % right velocity (leftward)

%%
for step = 1:(TR*L)
    % collision detection
    if norm(p_object - p_sun) < sun_radius || norm(p_object - p_earth) < earth_radius
        break
    end
    
    % define rotation matrix for step
    frame_theta = atan2(p_earth(2), p_earth(1));
    rotation = [cos(frame_theta), -sin(frame_theta), 0; sin(frame_theta), cos(frame_theta), 0; 0, 0, 1];
    p_object_rot = p_object * rotation;
    p_sun_rot = p_sun * rotation;
    p_earth_rot = p_earth * rotation;
    
    % store positions in 'path' matrices
    sun_path(step, 1) = p_sun(1);
    sun_path(step, 2) = p_sun(2);
    sun_path_rot(step, 1) = p_sun_rot(1);
    sun_path_rot(step, 2) = p_sun_rot(2);
    
    earth_path(step, 1) = p_earth(1);
    earth_path(step, 2) = p_earth(2);
    earth_path_rot(step, 1) = p_earth_rot(1);
    earth_path_rot(step, 2) = p_earth_rot(2);
    
    object_path(step, 1) = p_object(1);
    object_path(step, 2) = p_object(2);
    object_path_rot(step, 1) = p_object_rot(1);
    object_path_rot(step, 2) = p_object_rot(2);
    
    theta_store(step,1) = frame_theta;
    %compare craft to L4
    L4_p_difference_store(step, 1) = (au-449311)*cos(pi/3+frame_theta)-p_object(1);
    L4_p_difference_store(step, 2) = (au-449311)*sin(pi/3+frame_theta)-p_object(2);
    L4_p_difference_magnitude_store(step) = sqrt(L4_p_difference_store(step, 1)^2+L4_p_difference_store(step, 2)^2);
    L4_v_difference_store(step, 1) = 2*pi*(au-449311)/year*cos(5*pi/6+frame_theta)-v_object(1);
    L4_v_difference_store(step, 2) = 2*pi*(au-449311)/year*sin(5*pi/6+frame_theta)-v_object(2);
    L4_v_difference_magnitude_store(step) = sqrt(L4_v_difference_store(step, 1)^2+L4_v_difference_store(step, 2)^2);
    %compare craft to L5
    L5_p_difference_store(step, 1) = (au-449311)*cos(-pi/3+frame_theta)-p_object(1);
    L5_p_difference_store(step, 2) = (au-449311)*sin(-pi/3+frame_theta)-p_object(2);
    L5_p_difference_magnitude_store(step) = sqrt(L5_p_difference_store(step, 1)^2+L5_p_difference_store(step, 2)^2);
    L5_v_difference_store(step, 1) = 2*pi*(au-449311)/year*cos(pi/6+frame_theta)-v_object(1);
    L5_v_difference_store(step, 2) = 2*pi*(au-449311)/year*sin(pi/6+frame_theta)-v_object(2);
    L5_v_difference_magnitude_store(step) = sqrt(L5_v_difference_store(step, 1)^2+L5_v_difference_store(step, 2)^2);
    
    % accelerations on earth
    a_SunDueToEarth = G * m_earth * (p_earth - p_sun) / norm(p_earth - p_sun)^3;
    
    % accelerations on sun
    a_EarthDueToSun = G * m_sun * (p_sun - p_earth) / norm(p_sun - p_earth)^3;
    
    % accelerations on object
    a_ObjectDueToSun = G * m_sun * (p_sun - p_object) / norm(p_sun - p_object)^3;
    a_ObjectDueToEarth = G * m_earth * (p_earth - p_object) / norm(p_earth - p_object)^3;
    % calculate thrust based on profile
    v_hat = v_object / norm(v_object);
    v_rot = [cos(pi/2), -sin(pi/2), 0; sin(pi/2), cos(pi/2), 0; 1, 0, 0];
    a_ObjectThrust = a_0 * thrust_profile(step, :) * [v_hat; -v_hat; v_hat * v_rot; -v_hat * v_rot];
    a_ObjectSum = a_ObjectDueToSun + a_ObjectDueToEarth + a_ObjectThrust;
    
    % integrate accelerations
    v_object = v_object + a_ObjectSum * (1/TR);
    v_sun = v_sun + a_SunDueToEarth * (1/TR);
    v_earth = v_earth + a_EarthDueToSun * (1/TR);
    
    % integrate velocities
    p_object = p_object + v_object * (1/TR);
    p_sun = p_sun + v_sun * (1/TR);
    p_earth = p_earth + v_earth * (1/TR);
    
end

%% Display Minimums
[MinL4_separation,I_MinL4_separation] = min(L4_p_difference_magnitude_store(L4_p_difference_magnitude_store > 0));
[MinL4_rel_vel,I_MinL4_rel_vel] = min(L4_v_difference_magnitude_store(L4_v_difference_magnitude_store > 0));
[MinL5_separation,I_MinL5_separation] = min(L5_p_difference_magnitude_store(L5_p_difference_magnitude_store > 0));
[MinL5_rel_vel,I_MinL5_rel_vel] = min(L5_v_difference_magnitude_store(L5_v_difference_magnitude_store > 0));
fprintf(['\nMinimum L4 separation = ' num2str(MinL4_separation/au, '%10.2e') ' AU @ t = ' num2str((I_MinL4_separation-1)/TR) ' s going ' num2str(L4_v_difference_magnitude_store(I_MinL4_separation), '%.1f') ' m/s\n']);
fprintf(['Minimum L4 relative velocity = ' num2str(MinL4_rel_vel, '%.1f') ' m/s @ t = ' num2str((I_MinL4_rel_vel-1)/TR) ' s @ separation = ' num2str(L4_p_difference_magnitude_store(I_MinL4_rel_vel)/au, '%10.2e') ' AU\n\n']);
fprintf(['Minimum L5 separation = ' num2str(MinL5_separation/au, '%10.2e') ' AU @ t = ' num2str((I_MinL5_separation-1)/TR) ' s going ' num2str(L5_v_difference_magnitude_store(I_MinL5_separation), '%.1f') ' m/s\n']);
fprintf(['Minimum L5 relative velocity = ' num2str(MinL5_rel_vel, '%.1f') ' m/s @ t = ' num2str((I_MinL5_rel_vel-1)/TR) ' s @ separation = ' num2str(L5_p_difference_magnitude_store(I_MinL5_rel_vel)/au, '%10.2e') ' AU\n']);

% %% inertial reference frame
% figure(1);
% scatter(sun_path(:, 1), sun_path(:, 2), 100, 'y', 'filled');
% xlim([-2 * 10^11, 2 * 10^11])
% ylim([-2 * 10^11, 2 * 10^11])
% hold on;
% scatter(earth_path(:, 1), earth_path(:, 2), 9, 'b', 'filled');
% %plot path
% scatter(object_path(:, 1), object_path(:, 2), 2, 'k', 'filled')

%% rotating Sun-Earth reference frame
figure(2);
scatter(sun_path_rot(:, 1), sun_path_rot(:, 2) ,100, 'y', 'filled');
xlim([-2 * 10^11, 2 * 10^11])
ylim([-2 * 10^11, 2 * 10^11])
pbaspect([1, 1, 1]);
hold on;

%target ellipse
ellipse_xCenter = cos(pi/3)*(au-449311);
ellipse_yCenter = sin(pi/3)*(au-449311);
ellipse_xRadius = au*5*10^-3;
ellipse_yRadius = au*5*10^-2;
ellipse_theta = 0 : 0.01 : 2*pi;
ellipse_rot = [cos(pi/3),-sin(pi/3);sin(pi/3),cos(pi/3)]*[ellipse_xRadius * cos(ellipse_theta); ellipse_yRadius * sin(ellipse_theta)];
ellipse_rot(1,:) = ellipse_rot(1,:) + ellipse_xCenter;
ellipse_rot(2,:) = ellipse_rot(2,:) + ellipse_yCenter;
plot(ellipse_rot(1,:), ellipse_rot(2,:), 'LineWidth', 1);
plot(ellipse_rot(1,:), -ellipse_rot(2,:), 'LineWidth', 1);
scatter(earth_path_rot(:, 1), earth_path_rot(:, 2), 9, 'b', 'filled');
%plot path
scatter(object_path_rot(:, 1), object_path_rot(:, 2), 2, 'k', 'filled')

% %% animated plot
% 
% figure(3);
% 
% h = animatedline;
% h.Color = 'r';
% h.LineWidth = 1;
% i = animatedline;
% i.Color = 'k';
% i.LineWidth = 3;
% j = animatedline;
% j.Color = 'b';
% j.LineWidth = 2;
% 
% axis([-2 * 10^11, 2 * 10^11, -2 * 10^11, 2 * 10^11])
% 
% numpoints = TR*L;
% for k = 1:1:numpoints
%     addpoints(i,cos(theta_store(k, 1)+pi/3)*(au-449311), sin(theta_store(k, 1)+pi/3)*(au-449311)); % target zone
%     addpoints(h, object_path(k, 1), object_path(k, 2));
%     addpoints(j,earth_path(k, 1) ,earth_path(k, 2));
%     drawnow limitrate;
% end
