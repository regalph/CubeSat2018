
%% define intial conditions and constants
G = 6.674 * 10^-11; % N per kg^2 per m^2
m_sun = 1.989 * 10^30; % kg
p_sun = [-447.5, 0, 0]; % m, m
m_earth = 5.972 * 10^24; % kg
p_earth = [149597870252.5, 0, 0]; % m, m
rotating_frame = [0, 0, (2*pi)/31557600];

%% define time resolution and length
TR = .001; % samples per second
L = 600000000; % seconds per simulation

%% call new object function -> object(xpos,ypos,xvel,yvel)

p_object = [149597870252.5/2,149597870252.5/2*sqrt(3), 0];
v_object = [0, 3000, 0];


%% propogate and store
object_path = zeros(TR*L,3);
for step = 1:(TR*L);
    % object acceleration due to large bodies
    object_path(step,1)=p_object(1);
    object_path(step,2)=p_object(2);
    a_DueToSun = G * m_sun * (p_sun - p_object) / norm(p_object - p_sun)^3;
    a_DueToEarth = G * m_earth * (p_earth - p_object) / norm(p_object - p_earth)^3;
    a_Centripetal =  (2*pi*149597870252.5/31557600)^2*p_object;
    a_Coriolis = 2*cross(v_object, rotating_frame);
    a_Sum = a_DueToSun + a_DueToEarth + a_Centripetal + a_Coriolis;
    
    % update velocity with acceleration
    v_object = v_object + a_Sum*(1/TR);
    % update position with velocity
    p_object = p_object + v_object*(1/TR);
end
%% plot
figure(1);
scatter(p_sun(1),p_sun(2),100,'y','filled');
xlim([-2*10^11, 2*10^11])
ylim([-2*10^11, 2*10^11])
hold on;
scatter(p_earth(1),p_earth(2),9,'b','filled');
%plot path
scatter(object_path(:,1),object_path(:,2),2,'k','filled')