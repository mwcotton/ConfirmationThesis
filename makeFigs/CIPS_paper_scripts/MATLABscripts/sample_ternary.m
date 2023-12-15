function dataT = sample_ternary()

t_start = datestr(now,'yymmdd-HHMMSS');
data_file = ['../data-tables/',t_start,'-data-te.txt'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulation Parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=100;
Nx=1000;
T=10000;
Nt=200;
initnoise=0.01;

Ls=[];
Nxs=[];
Ts=[];
Nts=[];
initnoises=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% System Parameters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dse=1;
Dpe=4;
Dps=5;  %this isnt actually in the instability..?
Dew=10; %made this quite fast
Dsw=10; %made this quite fast
Dpw=10; %made this quite fast

kspo=1;
kcat=1;
e=-5;  %normally <0
mu=10;  %normally >0
vs=1;
ve=15;
vw=0.5;
phie0=0.15;
alpha=0.8;

s1 = [];
s2 = [];
s3 = [];
s4 = [];
Dses=[];
Dpes=[];
Dpss=[];
Dews=[];
Dsws=[];
Dpws=[];
kspos=[];
kcats=[];
es=[];
mus=[];
vss=[];
ves=[];
vws=[];
phie0s=[];
alphas=[];

%loop over some parameters
for alpha=0.8:0.007:0.9
    for phie0=0.01:0.0063:0.1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% Steady State %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Rs0=kspo.*exp(e)+kcat.*phie0;
        Rp0=kspo+kcat.*phie0.*exp(-e-mu);
        phis0=alpha.*Rp0./(Rs0+Rp0);
        phip0=alpha-phis0;
        phiw0=1-phie0-alpha;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% Do Numerics %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xmesh=linspace(0,L,Nx);
        tspan=linspace(0,T,Nt);
        sol = pdepe(0,@pdefun,@icfun,@bcfun,xmesh,tspan);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Plot Final State %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure('visible','off');
        hold all
        plot(xmesh,sol(end,:,1), color=[0, 0, 1, 0.5], linewidth=2)
        plot(xmesh,sol(end,:,2), color=[1, 0, 0, 0.5], linewidth=2)
        plot(xmesh,sol(end,:,3), color=[0, 0.5, 0, 0.5], linewidth=2)
        plot(xmesh,sol(end,:,4), color=[0, 0, 0, 0.5], linewidth=2)
        
        Rs=kspo*exp(e)+kcat.*sol(end,:,1);
        Rp=kspo+kcat.*sol(end,:,1)*exp(-e-mu);

        phis=(sol(end,:,2)+sol(end,:,3)).*Rp./(Rs+Rp);
        phip=(sol(end,:,2)+sol(end,:,3)).*Rs./(Rs+Rp);
        plot(xmesh,phis,'r--', linewidth=2);
        plot(xmesh,phip,color=[0, 0.5, 0],linewidth=2,linestyle='--');

        plot(xmesh,sol(1,:,1),color=[0, 0.5, 1], linestyle=':')
        plot(xmesh,sol(1,:,2),color=[1,  0,  0], linestyle=':')
        plot(xmesh,sol(1,:,3),color=[0, 0.5, 0], linestyle=':')
        plot(xmesh,sol(1,:,4),color=[0, 0, 0], linestyle=':')

        titl1 = sprintf('$k_{spo}: %.2g, k_{cat}: %.2g, D_{PS}: %.2g, D_{SE}: %.2g, D_{PE}: %.2g$', kspo, kcat, Dps, Dse, Dpe);
        titl2 = sprintf('$\\Delta e: %.2g, \\Delta\\mu: %.2g, v_S: %.2g, v_E: %.2g, \\phi_E: %.2g$', e, mu, vs, ve, phie0);
        titl3 = sprintf('$D_{EW}: %.2g, D_{SW}: %.2g, D_{PW}: %.2g, \\alpha: %.2g$', Dew, Dsw, Dpw, alpha);
        title({titl1, titl2, titl3},'Interpreter','latex');

        xlabel('x')
        ylabel('\phi')

        legend('E', 'S', 'P', 'W', '$S_{R=0}$', '$P_{R=0}$','Interpreter','latex');
        fig = gcf;

        FileName=['../Figures/sampling/n',datestr(now,'yymmdd-HHMMSS'),'.pdf'];

        saveas(fig,FileName)
        
        Dses=[Dses;Dse];
        Dpes=[Dpes;Dpe];
        Dpss=[Dpss;Dps];
        Dews=[Dews;Dew];
        Dsws=[Dsws;Dsw];
        Dpws=[Dpws;Dpw];

        kspos=[kspos;kspo];
        kcats=[kcats;kcat];
        es=[es;e];
        mus=[mus;mu];
        vss=[vss;vs];
        ves=[ves;ve];
        vws=[vws;vw];
        phie0s=[phie0s;phie0];
        alphas=[alphas;alpha];
        
        Ls=[Ls;L];
        Nxs=[Nxs;Nx];
        Ts=[Ts;T];
        Nts=[Nts;Nt];
        initnoises=[initnoises;initnoise];
        
        diff = zeros(1, 4);
        for compon=1:4
            diff(compon) = max(sol(end, :, compon))-min(sol(end, :, compon));
        end
        
        start_diff = zeros(1, 4);
        for compon=1:4
            start_diff(compon) = max(sol(1, :, compon))-min(sol(1, :, compon));
        end
        
        s1 = [s1;diff(1) > start_diff(1)];
        s2 = [s2;diff(2) > start_diff(2)];
        s3 = [s3;diff(3) > start_diff(3)];
        s4 = [s4;diff(4) > start_diff(4)];
        dataT = table(Dses,Dpes,Dpss,Dews,Dsws,Dpws,kspos,kcats,es,mus,vss,ves,phie0s,alphas,s1,s2,s3,s4,Ls,Nxs,Ts,Nts,initnoises)
        
        writetable(dataT,data_file,'Delimiter',',')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decalare Numeric Functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [c,f,s] = pdefun(x,t,u,dudx)
        
        Mpe = -Dpe*u(1)*u(3);
        Mse = -Dse*u(1)*u(2);
        Mps = -Dps*u(2)*u(3);
        Mps = -Dps*u(2)*u(3);
        Mps = -Dps*u(2)*u(3);
        Mps = -Dps*u(2)*u(3);

        Mes = Mse*ve/vs;
        Mep = Mpe*ve/vs;
        Msp = Mps;
        
        Mew = -Dew*u(1)*u(4);
        Msw = -Dsw*u(2)*u(4);
        Mpw = -Dpw*u(3)*u(4);
        
        Mwe = Mew*vw/ve;
        Mws = Msw*vw/vs;
        Mwp = Mpw*vw/vs;

        Mss = -(Mes+Mps+Mws);
        Mpp = -(Mep+Msp+Mwp);
        Mee = -(Mse+Mpe+Mwe);
        Mww = -(Mew+Msw+Mpw);

        dmue=dudx(1)/u(1);
        dmus=dudx(2)/u(2);
        dmup=dudx(3)/u(3);
        dmuw=dudx(4)/u(4);
        
        rspo=kspo*(exp(e)*u(2)-u(3));
        rcat=kcat*u(1)*(u(2)-(u(3)*exp(-e-mu)));
        R = rspo+rcat;
        
        c=ones(4,1);

        f=zeros(4,1);
        f(1) = Mee*dmue+Mes*dmus+Mep*dmup+Mew*dmuw;
        f(2) = Mse*dmue+Mss*dmus+Msp*dmup+Msw*dmuw;
        f(3) = Mpe*dmue+Mps*dmus+Mpp*dmup+Mpw*dmuw;
        f(4) = Mwe*dmue+Mws*dmus+Mwp*dmup+Mww*dmuw;
        
        s=zeros(4,1);
        s(1) = 0;
        s(2) = -R;
        s(3) = +R;
        s(4) = 0;
    end

    function u0 = icfun(x)
        u0=zeros(4,1);
        r_num = rand-0.5;
        u0(1) = phie0+initnoise*r_num;
        u0(2) = (1-initnoise*r_num)*phis0/(1-phie0);
        u0(3) = (1-initnoise*r_num)*phip0/(1-phie0);
        u0(4) = 1-u0(1)-u0(2)-u0(3);
        
        % in the limit of \phi near either 1 or 0, and larger
        % initnoise might be some funny business going on
    end

    function [pL,qL,pR,qR] = bcfun(xL,uL,xR,uR,t)
        pL=zeros(4,1);    
        pR=zeros(4,1);
        qL = [1;1;1;1];
        qR = [1;1;1;1];
    end


end