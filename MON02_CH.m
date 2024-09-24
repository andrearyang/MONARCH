clearvars;% close all;

%% Constants
sRrange=[0.0001 0.001 0.01 0.1 1 10 100 1000];  %Cas13d mRNA synthesis rate

plot_hill = true; %boolean operator, plot hill is true
plot_hf = true; %boolean operator, plot hf is true

sG=100;       %crRNA synthesis rate for model 1
gG=0.1;     %crRNA degradation rate for model 1
%sR=1;
gR=0.1;     %Cas13d mRNA degradation
sC=1;       %Cas13d protein synthesis rate
gC=0.01;    %Cas13d protein decay rate
gA=gC;      %activated Cas13d protein decay rate
gH=gC;      %hyperactivated Cas13d protein decay rate
sT=1;       %target RNA synthesis rate
gT=0.1;     %target RNA decay rate
sN=1;       %nontarget RNA synthesis rate
gN=0.1;     %nontarget RNA decay rate
a=1;        %Cas13d-crRNA bindng rate, complex activation
h=10;       %actCas13d-targetRNA bindng rate, complex hyperactivation
cT=10;      %set this to 0 to eliminate target cutting by hyperactive Cas13d
cN=10;
%c=10;       %targetRNA and nontargetRNA binding/cutting rate
gU=10*gR;   %unstable, poly(A)-less mRNA decay rate

T0=sT/gT*ones(size(sRrange)); %target RNA synthesis / target RNA decay should give steady state target RNA level.
N0=sN/gN*ones(size(sRrange));

%T0 tests for the steady state at different Cas13d mRNA synthesis rates.
%N0 tests for steady state of nontarget RNA

%% Constitutive guide
%% binding affinity * synthesis rate of crRNA * gG degradation rate * degradation rate * synthesis rate.
%% formation of Cas13d-crRNA complex * total degradation of both crRNA/Cas13d * binding interaction between cas13d  & rna.
C1t1=(a*sG + gG*gC - a*sC*sRrange/gR);
% binding rate * crRNA synthesis rate + decay * decay
C1t2=4*a*gC*gG*sC*sRrange/gR;

C1=( -C1t1 + sqrt(C1t1.^2 + C1t2) )/(2*a*gC);
G1=sG./(gG+a*C1);

h = 1;

eta1 = a.*C1.*G1;
T1t1 = h*eta1+gA*gT-sT*h;
T1t2 = 4*(cT*h*eta1*sT+sT*gA*gT*h); % cT instead of c?
T1t3 = 2*((cT*h*eta1/gA)+gT*h);
T1 = ( -T1t1 + sqrt(T1t1.^2 + T1t2) )./(T1t3);

A1 = eta1./(h*T1+gA);
H1 = h.*A1.*T1./gA;

% adding hill term if needed, and extra logic for cT == 0
if cT == 0 && plot_hill
    % BELOW: H=1
    T1_hill = [9.99001, 9.9001, 9.00112, 0.311317, 0.00111099, 0.000101332, 0.000101009, 0.000101009];
    N1_hill = [9.99004, 5.03247, 0.0200048, 0.00104302, 0.00101001, 0.00100991, 0.00100991, 0.00100991];

else
    cN1 = cN;
    N1 = sN ./ (gN+cN1.*H1);
    % MM for cT -> copy and pasted w/ help from Mathematica
    % BELOW: CT = 1, h = 1
    T1_hill = [9.89227, 5.19493, 0.191304, 0.0109042, 0.00100901, 0.000100407, 0.00010009, 0.00010009];
    % BELOW: h = 1, cn = 10, L = 10
    N1_hill = [9.99004, 5.03928, 0.0221257, 0.00198716, 0.00110395, 0.0010193, 0.00101927, 0.00101927];
end




% hfCas13d

T1t1_hf = (h/10)*eta1+gA*gT-sT*(h/10);
T1t2_hf = 4*((h/10)*(h/10)*eta1*sT+sT*gA*gT*(h/10)); % cT instead of c?
T1t3_hf = 2*(((h/10)*(h/10)*eta1/gA)+gT*(h/10));
T1_hf = ( -T1t1_hf + sqrt(T1t2_hf.^2 + T1t2_hf) )./(T1t3_hf);

A1_hf = eta1./((h/10)*T1_hf+gA);
H1_hf = (h/10).*A1_hf.*T1_hf./gA;
N1_hf = sN./(gN+0*H1_hf);

% hfCas13d_hill, using Mathematica (just adjusting rates as needed), eta is
% the same
T1_hfhill = [9.98041, 9.09167, 1.62626, 0.107975, 0.0100809, 0.00100398, 0.00100081, 0.00100081];
N1_hfhill = N0; % cN = 0

%% Coexpressed guide
%R2t1=( a*sC*sRrange/gU - a*gU*sRrange/gU - gR*gC);
R2t1=( a*sC*sRrange/gU - a*gU*sRrange/gU + gR*gC);
R2t2= 4*a*sRrange*gC*(sC-(sC-gU)*gR/gU) ;

R2=( -R2t1 + sqrt( R2t1.^2 + R2t2 ) )./( 2*a*( sC - (sC-gU)*gR/gU ) ) ;
%C2=(sRrange-gR*R2)/a./R2;
C2=(sRrange-gR*R2)./(a.*R2);
U2=(sRrange-gR*R2)./gU;

h = 10;

T2t1 = (gT*gA-h*sT+sRrange*h-gR*h*R2)*gH;
T2t2= 4*h*sT*gH*gA*(gT*gH+sRrange*cT-gR*cT*R2);

T2=( -T2t1 +sqrt(T2t1.^2+T2t2) )./( 2*h*(gT*gH+sRrange*cT-gR*cT*R2) );
A2=a*C2.*R2./(gA+h*T2);
H2=h*T2.*(sRrange-gR*R2)./(gH*gA+h*gH*T2);


eta2 = a.*C2.*R2;
% adding hill term if needed, and extra logic for cT == 0
if cT == 0 && plot_hill
    T2_hill = [9.99962, 9.9927, 9.90952, 9.03123, 0.160721, 0.000111501, 0.0000101112, 1.00132*10^-6];
    N2_hill = [10., 10., 9.9926, 5.26111, 0.0206188, 0.0199605, 0.0199602, 0.01996];
else
    cN2 = cN;
    N2 = sN ./ (gN+cN2.*H2);
    % MM for cT -> copy and pasted w/ help from Mathematica
    % BELOW: CT=10,h=10
    T2_hill = [9.99816, 9.49092, 1.16635, 0.0204198, 0.00110339, 0.000101231, 0.0000100191, 1.00041*10^-6];
    % BELOW: h = 10, cn = 1, L = 100
    N2_hill = [10., 10., 9.99262, 5.61387, 0.0898242, 0.0238518, 0.0203252, 0.0199963];

end

% hfCas13d

T2t1_hf = (gT*gA-(h/10)*sT+sRrange*(h/10)-gR*(h/10)*R2)*gH;
T2t2_hf= 4*(h/10)*sT*gH*gA*(gT*gH+sRrange*(h/10)-gR*(h/10)*R2);

T2_hf=( -T2t1_hf +sqrt(T2t1_hf.^2+T2t2_hf) )./( 2*(h/10)*(gT*gH+sRrange*(h/10)-gR*(h/10)*R2) );
A2_hf=a*C2.*R2./(gA+(h/10)*T2_hf);
H2_hf=(h/10)*T2_hf.*(sRrange-gR*R2)./(gH*gA+(h/10)*gH*T2_hf);
N2_hf=sN./(gN+0*H2_hf);

% hfCas13d_hill, using Mathematica (just adjusting rates as needed), eta is
% the same
T2_hfhill = [9.99947,9.94025,5.6681,0.200425,0.0110229,0.00101222,0.00010019,0.0000100041];

%% Plotting
figure;hold on;

% NOTE: Plotting structure is a tad messy, comment out the plot calls
% accordingly to get desired plots

plot(sRrange,T0,'r.','LineWidth',2,'MarkerSize',20);
plot(sRrange,N0,'ko','LineWidth',2);
if plot_hill
    plot(sRrange,T1_hill,'r','LineWidth',2);
    plot(sRrange,N1_hill,'k','LineWidth',2);
    %plot(sRrange,T2_hill,'r--','LineWidth',2);
    %plot(sRrange,N2_hill,'k--','LineWidth',2);
    if plot_hf
    plot(sRrange,T1_hfhill,'g','LineWidth',2);
    plot(sRrange,N1_hf,'c','LineWidth',2);
    %plot(sRrange,T2_hfhill,'g--','LineWidth',2);
    %plot(sRrange,N2_hf,'c--','LineWidth',2);
    end
else
    plot(sRrange,T1,'r','LineWidth',2);
    plot(sRrange,N1,'k','LineWidth',2);
    plot(sRrange,T2,'r--','LineWidth',2);
    plot(sRrange,N2,'k--','LineWidth',2);
    if plot_hf
    plot(sRrange,T1_hf,'g','LineWidth',2);
    plot(sRrange,N1_hf,'c','LineWidth',2);
    plot(sRrange,T2_hf,'g--','LineWidth',2);
    plot(sRrange,N2_hf,'c--','LineWidth',2);
    end
end
set(gca,'FontSize',20,'XScale','log','YScale','Log','XTick',[0.001 0.1 10])
xlabel('sR, Cas13d mRNA synth. rate');ylabel('levels');
%legend({'T0','T1','T2'},'Location','SW');
legend({'T0','N0','T1','N1','T1\_hf','N1\_hf'},'Location','SW')
%legend({'T0','N0','T1','N1','T2','N2'},'Location','SW')





