function Semi_Infinite_Convection()

temp = temp_at_time(400);
disp(temp)

end


function temp = temp_at_time(t)
x = 0;
h = 200;
k = 400;
alpha = 10^-4;
t_i = 25;
t_inf = 300;

temp = (erfc(x/2*sqrt(alpha*t)) - (exp((h*x/k)+((h^2*alpha*t)/(k^2))))*(erfc((x/2*sqrt(alpha*t))+(h*sqrt(alpha*t)/k))))*(t_inf-t_i) + t_i;

end