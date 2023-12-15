function dataT = sample_binary()

t_start = datestr(now,'yymmdd-HHMMSS');
data_file = ['./data-tables/',t_start,'-data-bi.txt'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulation Parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=100;
Nx=1000;
T=5000000;
Nt=20000;
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
Dpe=5;
Dps=50;  %this isnt actually in the instability..?
kspo=1;
kcat=1.25;
e=-5;  %normally <0
mu=5;  %normally >0
vs=1;
ve=20;
phie0=0.1;





s1 = [];
s2 = [];
s3 = [];
Dses=[];
Dpes=[];
Dpss=[];
kspos=[];
kcats=[];
es=[];
mus=[];
vss=[];
ves=[];
phie0s=[];

%loop over some parameters
for mu=3.5:0.5:8.5
    for phie0=0.08:0.08
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% Steady State %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Rs0=kspo.*exp(e)+kcat.*phie0;
        Rp0=kspo+kcat.*phie0.*exp(-e-mu);
        phis0=(1-phie0).*Rp0./(Rs0+Rp0);
        phip0=1-phie0-phis0;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% Do Numerics %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xmesh=linspace(0,L,Nx);
        tspan=linspace(0,T,Nt);
        sol = pdepe(0,@pdefun,@icfun,@bcfun,xmesh,tspan);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Plot Final State %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure('visible','on');
        hold all
        plot(xmesh,sol(end,:,1), color=[0, 0, 1, 0.5], linewidth=2)
        plot(xmesh,sol(end,:,2), color=[1, 0, 0, 0.5], linewidth=2)
        plot(xmesh,sol(end,:,3), color=[0, 0.5, 0, 0.5], linewidth=2)

        Rs=kspo*exp(e)+kcat.*sol(end,:,1);
        Rp=kspo+kcat.*sol(end,:,1)*exp(-e-mu);

        phis=(1-sol(end,:,1)).*Rp./(Rs+Rp);
        phip=(1-sol(end,:,1)).*Rs./(Rs+Rp);
        plot(xmesh,phis,'r--', linewidth=2);
        plot(xmesh,phip,color=[0, 0.5, 0],linewidth=2,linestyle='--');

        plot(xmesh,sol(1,:,1),color=[0, 0.5, 1], linestyle=':')
        plot(xmesh,sol(1,:,2),color=[1,  0,  0], linestyle=':')
        plot(xmesh,sol(1,:,3),color=[0, 0.5, 0], linestyle=':')


        titl1 = sprintf('$k_{spo}: %.2g, k_{cat}: %.2g, D_{PS}: %.2g, D_{SE}: %.2g, D_{PE}: %.2g$', kspo, kcat, Dps, Dse, Dpe);
        titl2 = sprintf('$\\Delta e: %.2g, \\Delta\\mu: %.2g, v_S: %.2g, v_E: %.2g, \\phi_E: %.2g$', e, mu, vs, ve, phie0);
        title({titl1, titl2},'Interpreter','latex');

        xlabel('x')
        ylabel('\phi')

        legend('E', 'S', 'P', '$S_{R=0}$', '$P_{R=0}$','Interpreter','latex');
        fig = gcf;

        FileName=['Figures/sampling/n',datestr(now,'yymmdd-HHMMSS'),'.pdf'];

        saveas(fig,FileName)
        
        Dses=[Dses;Dse];
        Dpes=[Dpes;Dpe];
        Dpss=[Dpss;Dps];
        kspos=[kspos;kspo];
        kcats=[kcats;kcat];
        es=[es;e];
        mus=[mus;mu];
        vss=[vss;vs];
        ves=[ves;ve];
        phie0s=[phie0s;phie0];
        
        Ls=[Ls;L];
        Nxs=[Nxs;Nx];
        Ts=[Ts;T];
        Nts=[Nts;Nt];
        initnoises=[initnoises;initnoise];
        
        diff = zeros(1, 3);
        for compon=1:3
            diff(compon) = max(sol(end, :, compon))-min(sol(end, :, compon));
        end
        
        s1 = [s1;diff(1) > initnoise];
        s2 = [s2;diff(2) > initnoise];
        s3 = [s3;diff(3) > initnoise];
        dataT = table(Dses,Dpes,Dpss,kspos,kcats,es,mus,vss,ves,phie0s,s1,s2,s3,Ls,Nxs,Ts,Nts,initnoises)
        
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

        Mes = Mse*ve/vs;
        Mep = Mpe*ve/vs;
        Msp = Mps;

        Mss = -(Mes+Mps);
        Mpp = -(Mep+Msp);
        Mee = -(Mse+Mpe);

        dmue=dudx(1)/u(1);
        dmus=dudx(2)/u(2);
        dmup=dudx(3)/u(3);
        
        rspo=kspo*(exp(e)*u(2)-u(3));
        rcat=kcat*u(1)*(u(2)-(u(3)*exp(-e-mu)));
        R = rspo+rcat;
        
        c=ones(3,1);

        f=zeros(3,1);
        f(1) = Mee*dmue+Mes*dmus+Mep*dmup;
        f(2) = Mse*dmue+Mss*dmus+Msp*dmup;
        f(3) = Mpe*dmue+Mps*dmus+Mpp*dmup;
        
        s=zeros(3,1);
        s(1) = 0;
        s(2) = -R;
        s(3) = +R;
        
    end

    function u0 = icfun(x)
        u0=zeros(3,1);
        r_num = rand-0.5;
        u0(1) = phie0+initnoise*r_num;
        u0(2) = (1-initnoise*r_num)*phis0;
        u0(3) = 1-u0(1)-u0(2);
        % in the limit of \phi near either 1 or 0, and larger
        % initnoise might be some funny business going on
    end

    function [pL,qL,pR,qR] = bcfun(xL,uL,xR,uR,t)
        pL=zeros(3,1);    
        pR=zeros(3,1);
        qL = [1;1;1];
        qR = [1;1;1];
    end


end