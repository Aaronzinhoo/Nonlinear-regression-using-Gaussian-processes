load temperature_data;
%%
%1 Fit Linear regression model to data
t_bar = mean(t);
%t_bar = ones(137,1)*t_bar;
x_bar = mean(x);
%x_bar = ones(137,1)*x_bar;
dist_t = t - t_bar;
dist_x = x - x_bar;
B_1 = sum(dist_t'*dist_x)/sum(dist_x.^2);
B_0 = t_bar - B_1*x_bar;

figure
hold on
scatter(x,t)
plot(x, B_0 + B_1*x)
xlabel('Year')
ylabel('Temp. Diff')
title('Temperature of Years w/ fitted line')
hold off

%COMMENT

%%
%2 Extrapolate the values to 2100
j = 1;
x_2 = zeros(221,1);
for i=1880:2100
    x_2(j) = i;
    t_2(j) = B_0(1) + B_1*x_2(j);
    j=j+1;
end
figure
hold on
scatter(x,t)
plot(x_2, t_2)
xlabel('Year')
ylabel('Temp. Diff')
title('Temperature of Years w/ fitted line')
hold off

% COMMENT

%%
%3 Code the Gaussian Kernal Function
%

%Gaussian Kernal
function k = kernal(x,x_hat,t)
    k = exp(-(x-x_hat)^2/(2*t^2));
end
