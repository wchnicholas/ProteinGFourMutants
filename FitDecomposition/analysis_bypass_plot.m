%% plot the frequency of success for two mechanism of bypass (Figure 2C,D)
%updated: 12/18/2015
clear;
clc;
load ./bypass/bypassfreq_1e5

%detour 
%if no missing genotype, count_all=19*(4-2)=38
freq_bypass=count_bypass./count_all;
%exclude entries if freq_bypass=NaN (because count_all=0)
freq_bypass=freq_bypass(~isnan(freq_bypass));

%conversion
%if no missing genotype, count_all_conversion=18*2=36, not 38!
% freq_bypass_conversion=count_bypass_conversion./count_all_conversion;
%fix bug, 12/18/2015 (substract 2 from count_all_conversion)
count_all_conversion_fixed=count_all_conversion-2;
freq_bypass_conversion=count_bypass_conversion./count_all_conversion_fixed;
%exclude entries if freq_bypass=NaN (because count_all=0)
freq_bypass_conversion=freq_bypass_conversion(~isnan(freq_bypass_conversion));

%% fit to poisson distribution
[lambdahat_detour,lambdaci_detour] = poissfit(freq_bypass)
[lambdahat_conversion,lambdaci_conversion] = poissfit(freq_bypass_conversion)


%% compare two mechanism: frequency
freq_max=max(max(freq_bypass),max(freq_bypass_conversion));
subplot(1,2,1);
hist(freq_bypass_conversion);
% set(gca,'xlim',[0 freq_max]);
set(gca,'ylim',[0 2.1e4],'fontsize',15);
xlabel('Bypass availability');
ylabel('Reciprocal epistasis');
title('Conversion');

subplot(1,2,2);
hist(freq_bypass);
% set(gca,'xlim',[0 freq_max]);
set(gca,'ylim',[0 2.1e4],'fontsize',15);
xlabel('Bypass availability');
ylabel('Reciprocal epistasis');
title('Detour');
% hold on;
% x = 0:0.001:freq_max;
% y = poisspdf(x,lambdahat_detour)*length(freq_bypass);
% plot(x,y,'--r')


%% histogram of successful bypass: frequency*possible routes
freq_success_bypass_conversion=freq_bypass_conversion(freq_bypass_conversion>0);
freq_success_bypass=freq_bypass(freq_bypass>0);

freq_max=max(max(freq_success_bypass),max(freq_success_bypass_conversion));
nbins=11;
subplot(2,1,1);
hist(freq_success_bypass_conversion*36); %fixed, 12/18/2015
% set(gca,'xlim',[0 freq_max]);
set(gca,'xlim',[1 22],'ylim',[0 7000]);
set(gca,'fontsize',15);
ylabel('Reciprocal epistasis');
title('Conversion');

subplot(2,1,2);
hist(freq_success_bypass*38); 
% set(gca,'xlim',[0 freq_max]);
set(gca,'xlim',[1 22],'ylim',[0 3000]);
set(gca,'fontsize',15);
xlabel('Available bypass');
ylabel('Reciprocal epistasis');
title('Detour');

%% bar plot: unavailable vs available
%plot the fraction of successful bypass
success_bypass_conversion=nnz(count_bypass_conversion)/length(count_bypass_conversion);

subplot(1,2,1);
b=bar([success_bypass_conversion,1-success_bypass_conversion]);
b(1).LineWidth = 2;
b(1).EdgeColor = 'b';
b(1).FaceColor = 'b';
box off;
set(gca,'XTickLabel',{'Success', 'Failure'})
set(gca','ylim',[0 1],'fontsize',20);
ylabel('Fraction');

success_bypass_detour=nnz(count_bypass)/length(count_bypass);
subplot(1,2,2);
b=bar([success_bypass_detour,1-success_bypass_detour]);
b(1).LineWidth = 2;
b(1).EdgeColor = 'b';
b(1).FaceColor = 'b';
box off;
set(gca,'XTickLabel',{'Success', 'Failure'})
set(gca','ylim',[0 1],'fontsize',20);
ylabel('Fraction');

