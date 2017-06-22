load temperature_data;
%%
%1 Fit Linear regression model to data
t_bar = mean(t);
%t_bar = ones(137,1)*t_bar; for 2016a version
x_bar = mean(x);
%x_bar = ones(137,1)*x_bar; 2016a version
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

%Although the model does fit the data ok it is missing the trends
%in the data that are characteristic of it. The data seems to have a
%curvature to it that the linear regression is not capturing which is the
%problem with this type of model for data of this type. This leads me to 
%assume there is a better model that can capture this relationship that is 
%is not being accounted for by LR.

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

% Extrapolation of points beyond the training data result in predictions
% that seem unreasonable simply because the data doesnt seem to be behaving
% this way at all. Since linear regression results in a predictions based
% on a linear relationship we are assuming the points have a linear trend
% which is simpy just not the case here. This makes it reasonable to assume
% that linear regression is not a good model to fit on this data.

%%
%3 Code the Gaussian Kernal Function
% Refer to coded function at the end

%%
%4 
% Create the C matrix
X = linspace(1880,2100,1000);
X =X';
tao = 10;
sigma = .01;
noise = sigma*eye(1000);
K = kernal(X,X,tao,size(X,1),size(X,1));
%C = K+noise;

mu = zeros(1000,1);
figure
for k=1:4
    hold on
    sample = mvnrnd(mu,K);
    plot(X,sample); 
    
end
xlabel('Year')
ylabel('Gaussian Prior Value')
title(sprintf('Prior Plot Sample w/ Realizations'))
hold off
figure
for k=1:4
    sample = mvnrnd(mu,K);
    subplot(2,2,k);
    plot(X,sample);
    xlabel('Year')
    ylabel('Gaussian Prior Value')
    title(sprintf('Prior Plot Sample (Realizations sep.)'))
end

% Here I have plotted the realizations together and then again seperately.
% We can see that the prior is actually a smooth curve of values which
% matches intuition since we are using a gaussian kernel to set up our
% covariance with. The values oscillate around the mean 0 and at tau^2 = 10
% seem to be very complex. I assume then that we may be overfitting at a
% value of 10 since we shouldnt be so sensitive between years. Plotting the
% data again with a higher tau we see the curve is fluctuating much less
% proving my earlier assumption. The plot is really complex with all of the
% realizations there (4) the other plot is less cluttered and lets us see
% the individual differences more clearly. But the point is to see them
% together exhibiting similar characteristics
%%
%5 Plot Predictive Dist of 1000 points w/ orig points
K = kernal(x,x,tao,size(x,1),size(x,1));
%K_X = zeros(size(K,1),1000);
mu = zeros(1000,1);
sig = zeros(1000,1);
C = K+sigma*eye(size(x,1));
C_inv = inv(C);
for i=1:1000
    k = kernal(X(i),x,tao,1,137);
    c = kernal(X(i),X(i),tao,1,1) + sigma;
    mu(i)= k*C_inv*t;
    sig(i)= c - k*C_inv*k'; 
end
figure
hold on
plot(X,mu);
plot(X,mu - 2*sqrt(sig),'--');
plot(X,mu + 2*sqrt(sig),':');
scatter(x,t,'.');
xlabel('Year')
ylabel('Temp Difference')
title('#5 Mean of the Data with SD +/-2')
legend('\mu' , '\mu - 2*\sigma', '\mu + 2*\sigma')

% The distribution does well fitting the data but doesnt do well when
% extrapolating points since the range where a point lies within the
% extrapolated region is wide. Meaning that we are confident a point lies 
% within a larger range of values than maybe what it can be. The standard
% deviation seems to drift quickly from the mean after we get passed 2050 
% and the mean seems to approach 3 rather quickly. The confidence interval
% being large is similar to us being confident that someones height lies
% between 4ft to 8ft which for the most case is obvious and doesnt give any
% relevant predition of info. Since the modeling of the data could be due
% to overfitting as mentioned earlier this might actully be bad for us and
% the reason why we are getting bad predictions since we overfit the
% training data and now fail on the testing data. Unlike linear regression
% this model captures trends in the data.

%%
%6 Tune width of the gaussian model

figure
range = [10,1000,10000,100000];
j=1;
k=1;
% plot realizations again
for tau=range
    K = kernal(X,X,tau,size(X,1),size(X,1));
    mu = zeros(1000,1);
    subplot(2,2,k);
    hold on
    for l=1:4
        sample = mvnrnd(mu,K);          
        plot(X,sample);
    end
    hold off
    xlabel('Year')
    ylabel('Gaussian Prior Value')
    title(['Prior Plot Sample \tau^2=',num2str(tau)])
    k=k+1;
end
% plot the extrapolated points
figure
for tau=range
    K = kernal(x,x,tau,size(x,1),size(x,1));
    %K_X = zeros(size(K,1),1000);
    mu = zeros(1000,1);
    sig = zeros(1000,1);
    C = K+sigma*eye(size(x,1));
    C_inv = inv(C);
    for i=1:1000
        k = kernal(X(i),x,tau,size(x,2),size(x,1));
        c = kernal(X(i),X(i),tau,1,1) + sigma;
        mu(i)= k*C_inv*t;
        sig(i)= c - k*C_inv*k'; 
    end
    subplot(2,2,j);
    hold on
    plot(X,mu);
    plot(X,mu - 2*sqrt(sig),'--');
    plot(X,mu + 2*sqrt(sig),':');
    scatter(x,t,'.');
    xlabel('Year')
    ylabel('Temp Difference')
    title(['#6 Mean of the Data w/ SD +/- 2 tau^2=',num2str(tau)])
    legend('\mu' , '\mu - 2*\sigma', '\mu + 2*\sigma')
    j=j+1;
    hold off
end

% We see that on the evenly spaced points from the previous question that
% as tau increases or decreases the comnplexity of the GP prior changes in
% complexity. The higher the value of tau the less sensitive it is to the
% data making the curve less jagged and more smother. You can see the
% difference clearly when tau =10 and tau =10000 that the 10k value almost
% gives stright lines. Tau's effect on the mean of the distribution is that
% as tau decreases or approaches zero the kernal of the points approach
% 0 since tau is in the denom of the gaussian kernal. 
% This simulataneously occuring for small values of tau, as we
% extrapolate values that are much larger than the training data. The
% greater the value of tau also means the less close fit we have to the
% training data and extrapolated values are less likely to revert to zero
% as quickly since at large values of tau the mean is kept positive. Since 
% we are plotting the mean with 2 SDs we capture 95% of where the data
% should fall within this area. For values of low tau i.e. 10 the
% range of possible values for extrapolated points lies within in a larger
% range and therefore we are not providing very meaningful predictions this
% way and we are becoming  more uncertain as to what the value actually is. 
% The higher the value of tau i.e. 10k smaller this gap becomes and 
% the more and this results in us providing more meaningful predictions at 
% the cost of not modeling the original training points as well. The values
% of extrapolated points at small tau revert to zero quickly as discussed
% above due to the design of the kernal and the affect of tau on it. The
% higher values keep consistent with the curve and allow us to get good
% predictions assuming our assumptions on the data are reasonable.

%Note, tau models the complexity of the model and accounts for the
%bias-variance tradeoff in ML terms. Also as tau increases to infinity the
%extrapolation simply becomes a line since exp -> 1.

%%
%7 Tune Tau to get good recovery of the data 
% 50000 35000 20000 10000
range = [10000 15000 20000 25000];
j=1;
figure
for tau=range
    K = kernal(x,x,tau,size(x,1),size(x,1));
    %K_X = zeros(size(K,1),1000);
    mu = zeros(1000,1);
    sig = zeros(1000,1);
    C = K+sigma*eye(size(x,1));
    C_inv = inv(C);
    for i=1:1000
        k = kernal(X(i),x,tau,size(x,2),size(x,1));
        c = kernal(X(i),X(i),tau,1,1) + sigma;
        mu(i)= k*C_inv*t;
        sig(i)= c - k*C_inv*k'; 
    end
    %plot(X,mu - 2*sqrt(sig));
    %plot(X,mu + 2*sqrt(sig));
    
    subplot(2,2,j);
    hold on
    plot(X,mu);
    plot(X,mu - 2*sqrt(sig),'--');
    plot(X,mu + 2*sqrt(sig),':');
    scatter(x,t,'.');
    xlabel('Year')
    ylabel('Temp Difference')
    title(['#7 Mean of the Data with \tau^2 = ',num2str(tau)])
    legend('\mu' , '\mu - 2*\sigma', '\mu + 2*\sigma')
    j=j+1;
    hold off
end

% Comment on the value of tau and make sure to notice as t -> 0 the values
% approach infty which mean that the the mean approaches 0 for values of
% x-x' larger than tau. An other thing to notice is the fact that as tau
% approaches infty we get the fact that e^(kernal) goes to zero resulting
% in a line for the mean that is dependent entirely on C.

% The reason it is hard to extrapolate values that are greater than tau is
% because tau is in the denominator of an exponent of e. 
% e^(-(x-x')^2/tau^2) which means for tau less than the extrapolated
% values or in this case x' s.t. x>>x' or  x<<x' then e will more quickly
% approach zero since the value becomes a really large negative power
% quickly. Since it quickly goes to zero then this causes a problem for
% extrapolating points where tau is less than the difference from the
% test and training points since this value i.e. the kernal controls the
% value of the mean and the standard deviation. So to get reasonable values
% for 2100 we should ensure tau is greater than the difference between the
% values. Since we extrapolating to 2100 we want tau greater than the 
% difference of 2100 -2016 = 84. 
% Setting tau above this value should result in reasonable predictions for
% a range beyond this as well. I map out a few values of tau above the one
% mentioned and each do much better than 10, but it seems that 20k does the
% best at capturing the original data and having little uncertainty of
% where the data may lie. 10k also seems good but as it gets closer to 2100
% our uncertainty increases a lot which may be a problem. Values above 20k
% seem to be the same in terms of their uncertainty and their mean up to 
% 2100. 10k and 20k both share a similar mean as well. 20k is the best
% value i believe since the uncertainty is not too large and the values 
% seem to be reasonable considering the data. I already commented on the
% value above. But the farther the value of x gets from the data the more
% uncertainty there is which explains the the uncertainty interval around
% 2100 for 20k. 20k also does good in recovering the data while
% simultaneously extrapolating points making it a good choice overall

%%
%8 2500 1600 800 500 for 8 
% 
range = [30000 40000 50000 60000];
range_3 = [100 500 1000 2500];
%X = linspace(1880,2100,1000);
X_2 = linspace(1880,2200,1000);
X_3 = linspace(1880,2022,1000);
j=1;
figure
for tau=range
    K = kernal(x,x,tau,size(x,1),size(x,1));
    %K_X = zeros(size(K,1),1000);
%    mu = zeros(1000,1);
%    sig = zeros(1000,1);
    mu_2 = zeros(1000,1);
    sig_2 = zeros(1000,1);
    C = K+sigma*eye(size(x,1));
    C_inv = inv(C);
    for i=1:1000
%         k = kernal(X(i),x,tau,size(x,2),size(x,1));
%         c = kernal(X(i),X(i),tau,1,1) + sigma;
%         mu(i)= k*C_inv*t;
%         sig(i)= c - k*C_inv*k';
        k = kernal(X_2(i),x,tau,size(x,2),size(x,1));
        c = kernal(X_2(i),X_2(i),tau,1,1) + sigma;
        mu_2(i)= k*C_inv*t;
        sig_2(i)= c - k*C_inv*k';
     end
%     %plot(X,mu - 2*sqrt(sig));
%     %plot(X,mu + 2*sqrt(sig));
%     %subplot(5,2,j);
%     figure
%     hold on
%     plot(X,mu);
%     scatter(x,t,'.');
%     plot(X,mu - 2*sqrt(sig));
%     plot(X,mu + 2*sqrt(sig));
%     xlabel('Year')
%     ylabel('Temp Difference')
%     title(['Mean of the Data with \tau^2 = ',num2str(tau)])
%     hold off
    subplot(2,2,j);
    hold on
    plot(X_2,mu_2);
    plot(X_2,mu_2 - 2*sqrt(sig_2),'--');
    plot(X_2,mu_2 + 2*sqrt(sig_2),':');
    scatter(x,t,'.');
    xlabel('Year')
    ylabel('Temp Difference')
    title(['Mean of the Data to 2200 with \tau^2 = ',num2str(tau)])
    legend('\mu' , '\mu - 2*\sigma', '\mu + 2*\sigma')
    j=j+1;
    hold off
end
% for values closer to the range 2022
figure
j = 1;
for tau=range_3
    K = kernal(x,x,tau,size(x,1),size(x,1));
    mu = zeros(1000,1);
    sig = zeros(1000,1);
    C = K+sigma*eye(size(x,1));
    C_inv = inv(C);
    for i=1:1000
        k = kernal(X_3(i),x,tau,size(x,2),size(x,1));
        c = kernal(X_3(i),X_3(i),tau,1,1) + sigma;
        mu(i)= k*C_inv*t;
        sig(i)= c - k*C_inv*k'; 
    end
    subplot(2,2,j);
    hold on
    plot(X_3,mu);
    plot(X_3,mu - 2*sqrt(sig),'--');
    plot(X_3,mu + 2*sqrt(sig),':');
    scatter(x,t,'.');
    xlabel('Year')
    ylabel('Temp Difference')
    title(['Mean of the Data till 2022 with \tau^2 = ',num2str(tau)])
    legend('\mu' , '\mu - 2*\sigma', '\mu + 2*\sigma')
    xlim([1880 2022])
    hold off
    j = j+1;
end

% Conclusion: For extrapolating values up to 2200 using a similar idea from
% before we want values above 33k and the best one that seems to fit our
% mark is around 50k since the uncertainty is moderately low. The method 
% isnt well proven it just makes sense intuitively considering the
% parameters and the behaviour of or given kernel. 10k falls to quikly to
% zero and 30 seems to have high uncertainty compared to 50k.
% For 2100 as previously stated I would consider 20k as a good value that
% captures the data and predict the values well and with little
% uncertainty.
% For only a few years we want tau to be low that we capture the change
% within the training data as much as possible while still accounting for
% reasonable values when extrapolating. The method gives a value of 36
% which is too low and result in a lot of uncertainty for values
% immediately after the trianing range since we overfit the data so I checked
% values up to 2500 and it seems that best value occurs at around 1000 since
% we are still modeling the training data well and giving values that lie
% within a reasonable uncertainty range.

%%
%Gaussian Kernal for the entire Cov matrix
function K = kernal(x,x_hat,t_square,len1,len2)
K = zeros(len1,len2);
for i=1:len1
    for j=1:len2
        K(i,j) = exp(-((x(i) - x_hat(j))^2)/(2*t_square));
    end
end
end
