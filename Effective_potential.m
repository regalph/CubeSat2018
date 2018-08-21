close all
clear vars
G = 6.675452 * 10^-11; % N*m^2 per kg^2
au = 149597870700;
year = 31557600; % s, Julian Year

m_sun = 1.98855 * 10^30; % kg
m_earth = 5.97237 * 10^28; % kg

p_sun = [-449311, 0, 0]; % m
p_earth = [149597421389, 0, 0]; % m
p_object = [0, 0, 0];
range = 2*10^11;
step_size = 1*10^9;
U_per_mass_store = zeros(range/step_size+1, range/step_size+1);
 



for i = -range:step_size:range
    p_object(1)=146759062430+i;
    for j = -range:step_size:range
        p_object(2)=-223897479+j;
%        T = 2*pi*sqrt((norm(p_object)^3)/(G*m_sun));
%        L_per_mass = norm(p_object)*(2*pi*norm(p_sun-p_earth)/T);
%        U_per_mass_effective = (1/2)*(L_per_mass^2)/(norm(p_object)^2)-G*m_sun/(norm(p_object-p_sun))-G*m_earth/(norm(p_object-p_earth));
        U_per_mass_effective = (-G/2)*(m_sun+m_earth)*(norm(p_object)^2)/(norm(p_sun-p_earth)^3)-G*m_sun/(norm(p_object-p_sun))-G*m_earth/(norm(p_object-p_earth));
        U_per_mass_store((j+range+step_size)/step_size,(i+range+step_size)/step_size) = U_per_mass_effective;
    end
end

%%
U_per_mass_store(U_per_mass_store > 10^10 |U_per_mass_store < -1.5*10^9 ) = NaN;
figure(1);
surf(U_per_mass_store,'EdgeColor','no ne');
%%
figure(2);
[dx,dy] = gradient(U_per_mass_store);
surf(sqrt(dx.^2+dy.^2),'EdgeColor','none');
%set(gca, 'ZScale', 'log');

%%
figure(3); 
contour(U_per_mass_store,100);
%%
figure(4);
contour(U_per_mass_store.^2,100);