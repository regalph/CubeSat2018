%% define intial conditions and constants
G = 6.674 * 10^-11; % N per kg^2 per m^2
m_sun = 1.989 * 10^30; % kg
p_sun = [447.5, pi, 0]; % m, rad
v_sun = [0,(2*pi)/31557600,0]; % rad per sec
m_earth = 5.972 * 10^24; % kg
p_earth = [149597870252.5, 0, 0]; % m, rad
v_earth = [0,(2*pi)/31557600,0]; % rad per sec


%% define time resolution and length
TR = .1; % samples per second
L = 600000; % seconds per simulation

%% call new object function -> object(xpos,ypos,xvel,yvel)

p_object = [149597870252.5,pi/3, 0];
v_object = [0, (2*pi)/31557600, 0];

sun_path = zeros(TR*L,3);
earth_path = zeros(TR*L,3);
object_path = zeros(TR*L,3);

for step = 1:(TR*L);
    sun_path(step,1)=p_sun(1);
    sun_path(step,2)=p_sun(2);
    earth_path(step,1)=p_earth(1);
    earth_path(step,2)=p_earth(2);
    object_path(step,1)=p_object(1);
    object_path(step,2)=p_object(2);
    
    a_DueToSun = G * m_sun * (p_sun - p_object) / norm(p_object - p_sun)^3;
    a_DueToEarth = G * m_earth * (p_earth - p_object) / norm(p_object - p_earth)^3;
    a_Sum = 0;
    
    v_object = v_object + a_Sum*(1/TR);
    
    p_object = p_object + v_object*(1/TR);
    p_sun = p_sun + v_sun*(1/TR);
    p_earth = p_earth + v_earth*(1/TR);
        
end

%% plot
figure(1);
polarscatter(sun_path(:,2),sun_path(:,1),100,'y','filled');
rlim([0, 2*10^11]);
hold on;
polarscatter(earth_path(:,2),earth_path(:,1),9,'b','filled');
%plot path
polarscatter(object_path(:,2),object_path(:,1),2,'k','filled')