function [xmesh,tspan,sol] = specify_binary(L, Nx, T, Nt, initnoise, Dse, Dpe, Dps, kspo, kcat, e, mu, vs, ve, phie0, figsave, figshow)

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
if figshow
    figure('visible','on');
else
    figure('visible','off');
end

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

legend('E', 'S', 'P', '$S(\phi_E)$', '$P(\phi_E)$','Interpreter','latex');
fig = gcf;

FileName=['../Figures/numsol',datestr(now,'yymmdd-HHMMSS'),'.pdf'];

if figsave
    saveas(fig,FileName)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save end solution to file %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = table(sol(end, :, 1)',sol(end, :, 2)',sol(end, :, 3)', xmesh', 'VariableNames', {'Enzymes', 'Substrate', 'Product', 'x'});
TableFile=['../numerics-out/',datestr(now,'yymmdd-HHMMSS'),'.txt'];
writetable(T, TableFile)

sys1 = sprintf('k_{spo}: %.2g, k_{cat}: %.2g, D_{PS}: %.2g, D_{SE}: %.2g, D_{PE}: %.2g', kspo, kcat, Dps, Dse, Dpe);
sys2 = sprintf('Delta e: %.2g, Delta mu: %.2g, v_S: %.2g, v_E: %.2g, phi_E: %.2g', e, mu, vs, ve, phie0);

sys_pars = strcat('\n', sys1, ', ', sys2);
fid = fopen(TableFile, 'a');
fprintf(fid, sys_pars);
fclose(fid);

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