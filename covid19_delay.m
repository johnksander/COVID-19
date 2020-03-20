clear
clc; close all
format compact

%R0 "basic reproduction number" https://en.wikipedia.org/wiki/Basic_reproduction_number
options.R0 = 2.28; % value from https://www.ijidonline.com/action/showPdf?pii=S1201-9712%2820%2930091-6
options.sigR0 = .25; %variance in case-specific reproductive values
options.Vlife = 14; %14 day sickness
%https://en.wikipedia.org/wiki/Basic_reproduction_number#Reproductive_number_as_it_relates_to_contact_rate_and_infectious_period
%this R0 is new cases per 14 days with no effective quarantine

%symptom onset probabilities from: annals.org/aim/fullarticle/2762808/incubation-period-coronavirus-disease-2019-covid-19-from-publicly-reported
options.symp_pd = makedist('Lognormal','mu',1.621,'sigma',0.418); %values reported https://github.com/HopkinsIDD/ncov_incubation

options.I0 = 10442; %start with this many infections
%www.cdc.gov/coronavirus/2019-ncov/cases-updates/cases-in-us.html?CDC_AA_refVal=https%3A%2F%2Fwww.cdc.gov%2Fcoronavirus%2F2019-ncov%2Fcases-in-us.html
%CDC says 10,442 cases as of 3/19


fz = 18; %font

%look at some distributiopns  
figure(1);subplot(2,1,1)
histogram(normrnd(options.R0,options.sigR0,10e3,1),'Normalization','probability')
ylabel('p(x)')
xlabel('new infections per case ( R_0 )')
set(gca,'FontSize',fz)

x = 0:20;y = cdf(options.symp_pd,x);
subplot(2,1,2)
plot(x,y,'LineWidth',2);
ylabel({'proportion';'symptomatic'});xlabel('days after infection')
set(gca,'FontSize',fz)
%Compare to this fig https://github.com/HopkinsIDD/ncov_incubation/blob/master/README_files/figure-markdown_strict/dic-plots-1.png


iso_delays = [0,.5,1.5,options.Vlife];
options.Tmax = 360; %1 year 
num_delay = numel(iso_delays);
Inew = cell(num_delay,1);
Itotal = cell(num_delay,1);

Tday = 1:options.Vlife:options.Tmax;
Tn = numel(Tday);

cols = parula;
cols = cols(round(linspace(1,numel(cols(:,1)),num_delay)),:);

figure(2);
for idx = 1:num_delay
    
    options.delay = iso_delays(idx);
    fprintf('\n:::doing delay = %i hrs:::\n\n',options.delay*24)
   
    [Inew{idx},Itotal{idx}] = do_spread(options);
    
    pause(.1)
    ax(idx) = subplot(2,2,idx);
    if Itotal{idx}(end) > 1e6
        sc = 1e6;
        sclab = 'million';
    else
        sc = 10e3;
        sclab = '10k';
    end
    plot(Tday,Inew{idx}./sc,'LineWidth',2);hold on
    plot(Tday,Itotal{idx}./sc,'LineWidth',3,'LineStyle','--');
    xlabel('day','FontWeight','b')
    ylabel(sprintf('%s infections',sclab),'FontWeight','b')
    if options.delay < options.Vlife
        tit = sprintf('isolation %i hrs after S_x',options.delay*24);
    else
        tit = 'no isolation';
    end
    title(tit,'FontWeight','b')
    axis tight;drawnow

end

axes(ax(1))
legend({'new','total'},'FontWeight','b',...
    'Location','northwest','Box','off','Orientation','vertical')


function [Inew,Itotal] = do_spread(options)
USpop = 327e6; %US population
Vlife = options.Vlife;
Tday = 1:Vlife:options.Tmax;
Tn = numel(Tday);

Inew = zeros(Tn,1);
Itotal = zeros(Tn,1);

delay = options.delay; %delay from symptom onset to isolation

p_symp = cdf(options.symp_pd,1:Vlife); %per-day likelihood you show symptoms after onset
p_symp(end) = 1; %just to simplify
R0 = options.R0; sigR0 = options.sigR0;

Inew(1) = options.I0;
Itotal(1) = options.I0;
for idx = 2:Tn
    
    n = Inew(idx-1);
    if delay < Vlife
        Rx_time = rand(n,1);
        Rx_time = arrayfun(@(x) find(x <= p_symp,1),Rx_time);
        Rx_time = Rx_time + delay;
        Rx_time(Rx_time > Vlife) = Vlife; %can't go over the timeperiod
        R0_t = R0 * (Rx_time/Vlife);%adjusted R0 for time prior to treatment
    else
        R0_t = R0 + zeros(n,1);
    end
    Rt = normrnd(R0_t,sigR0,n,1); %draw r0 for each each
    Rt(Rt < 0) = 0;
    It = sum(round(Rt)); %new infections in this period
    
    Inew(idx) = It;
    Itotal(idx) = It + Itotal(idx-1);
    if Itotal(idx) >= USpop
        fprintf('----everyone infected in %i days\n',idx*Vlife)
        Itotal(idx:Tn) = USpop;
        break
    end
    if mod(idx,3) == 0,fprintf('week %i\n',(idx*Vlife)/7);end
end
end