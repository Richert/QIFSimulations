% initialization
%%%%%%%%%%%%%%%%

clear;                           % clear variables
format compact
close all;                       % close figures
addpath('../ddebiftool-git/ddebiftool',...
    '../ddebiftool-git/ddebiftool_extra_psol',...
    '../ddebiftool-git/ddebiftool_extra_nmfm',...
    '../ddebiftool-git/ddebiftool_utilities');

% definition of right-hand side evaluation of simple qif population
qif_rhs = @(ys, pars) qif_syn_adapt_rhs(ys, pars);

% definition of time-delay positions in parameter vector
qif_tau=@()[11,12];

% definition of continuation parameter range
iv=3;
iv_min=0.99;
iv_max=5.0;
iv_delta=0.05;

% initialization of user-defined input structure
funcs=set_funcs('sys_rhs',qif_rhs,'sys_tau',qif_tau);

% continuation in eta
%%%%%%%%%%%%%%%%%%%%%

% branch setup
branch1 = SetupStst(funcs,...
    'x',[20.17;-2.63;60.18;-0.378;285.95;683.69;962.12;458.00;35.98;0.06],...
    'parameter',[-4.39, 0.31, 1.0, 14.18, 11.36, 47.71, 7.61, 2.99, ...
    0.006, 0.014, 0.004, 0.006, 0.002, 0.004, 0.2, 2.0],...
    'contpar',iv,'max_step',[iv iv_delta],'max_bound',[iv iv_max],...
    'min_bound',[iv iv_min],'newheuristics_tests',0);

% continuation
branch1.method.continuation.plot = 0;
[branch1,s,f,r] = br_contn(funcs,branch1,1000);
branch1=br_rvers(branch1);
[branch1,s,f,r]=br_contn(funcs,branch1,1000); 

% stability calculation
branch1.method.stability.minimal_real_part = -10.0;
%branch1.method.stability.minimal_time_step = 0.000001;
[nunst_stst,dom,triv,branch1.point]=GetStability(branch1,'funcs',funcs);

% bifurcation detection
branch1.method.bifurcation.minimal_real_part = -10.0;
%branch1.method.bifurcation.imagthreshold = 1e-8;
branch1.method.point.newton_max_iterations = 10;
branch1.method.point.newton_nmon_iterations = 4;
branch1.method.point.halting_accuracy = 1e-8;
branch1.method.point.minimal_accuracy = 1e-6;
branch1 = LocateSpecialPoints(funcs, branch1);

% plotting I
%%%%%%%%%%%%

% bifurcation diagram
figure(1);
[xm,ym] = df_measr(0,branch1,1);
br_plot(branch1,xm,ym);
xlabel('k');
ylabel('r_e');

% eigenvalue plot
figure(2);
[xm,ym]=df_measr(1,branch1);
yyaxis left
br_plot(branch1,xm,ym,'b');
xlabel('k');ylabel('\Re\lambda');
yyaxis right
ym.func='imag';
br_plot(branch1,xm,ym,'r');
xlabel('k');ylabel('\Im\lambda');

% continuation of hopf in eta-k_ee
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get index of hopf and second continuation parameter
FPI = br_getflags(branch1);
idx_hopf = FPI(bif2num('hopf'),1);
iv2 = 2;

% setup hopf branch
[branch2, s] = SetupHopf(funcs, branch1, idx_hopf, 'contpar', [iv iv2], 'dir', iv, 'step', 0.0001);
branch2.parameter.min_bound(1:2,:)=[[iv iv_min]' [iv2 -100.0]']';
branch2.parameter.max_bound(1:2,:)=[[iv iv_max]' [iv2 1.0]']';
branch2.parameter.max_step(1:2,:)=[[iv iv_delta]' [iv2 0.05]']';
branch2.method.continuation.plot = 0;
branch2.method.stability.minimal_real_part = -10.0;

% continue hopf branch
[branch2,s,f,r]=br_contn(funcs,branch2,1000);            
branch2=br_rvers(branch2);
[branch2,s,f,r]=br_contn(funcs,branch2,1000);          

% bifurcation detection
branch2 = LocateSpecialPoints(funcs, branch2);

% plotting
figure(3); clf;

bifplot = df_bifplot();
bifplot.hoho = '*';

[xm,ym] = df_measr(0,branch2,1);
br_plot(branch2,xm,ym,'-');
br_bifplot(branch2,xm,ym,bifplot);

xlabel('k'); ylabel('\eta_i');

% continuation of hopf in eta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[branch3,suc]=SetupPsol(funcs,branch1,idx_hopf,'degree',4,'intervals',40);

% continuation
branch3.parameter.max_step(1,:)=[iv, iv_delta*10.0];
branch3.method.continuation.plot = 0;
[branch3,s,f,r]=br_contn(funcs,branch3,4000);        

% stability extraction
branch3.method.stability.minimal_real_part = -500.0;
%branch1.method.stability.minimal_time_step = 0.000001;
branch3 = br_stabl(funcs, branch3, 0, 0);

% bifurcation detection
branch3 = br_bifdet(funcs, branch3);

% look at the period along the branch:
etas=arrayfun(@(x)x.parameter(iv),branch3.point);
periods=[branch3.point.period];
figure(4); clf;
plot(etas,periods,'b.-');
xlabel('k');
ylabel('period');

% finish bifurcation diagram
figure(1);

bifplot = df_bifplot();
bifplot.hoho = '*';
[xm,ym] = df_measr(0,branch1,1);
br_bifplot(branch1,xm,ym,bifplot);

[xm,ym] = df_measr(0,branch3,1);
br_plot(branch3,xm,ym,'--');
br_bifplot(branch3,xm,ym,bifplot);

% eigenvalue plot
figure(5);clf;
%[xm,ym]=df_measr(1,branch3);
for i=1:length(branch3.point)
    p_splot(branch3.point(i))
end
%br_plot(branch3,xm,ym,'b');
