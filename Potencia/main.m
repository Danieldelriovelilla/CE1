%{
    Programa que resuelve el Trabajo 1 de la asignatura
    Generacion y Gestion de Potencia Eelectrica
%}

clc;
clear all;
close all;
fig = 1;


%% DATOS

% Tierra
mu = 398600;                % km^3/s^2
rT = 6378;                  % km
J2 = 1.0827*10^-3;          % -

% Sol
beta = deg2rad(0);                  % rad (angulo beta -> sol/tierra)
beta_v = [cos(beta) sin(beta) 0];   % versor solar
G = 1361 ;                          % W/m2

% Orbita
h = [420];        % km
r = rT + h;                 % km
RAAN = deg2rad(0);         % rad

% Satelite
%w = [0.05, 0.1, 0.5];       % rad/s
w = [0];       % rad/s
%   A = (30.18 + 11*2.17)/100^2;
  A = (30.18*2)/100^2;
% m^2 Area 1U AZURSPACE
% A = 4.9290*12/100^2;      % m^2 Area 1U EspectroLab TriSol X   
fc = 1;                   % factor de ocupacion REF : c.pdf ( pindado) pag 10.
rend = 0.295;              % aprox valor DataSheet Azure triple joint



%% CALCULO INCLINACION

cte = 2*pi/(365.25*24*3600);
inc = acos(((-3*rT^2*J2*mu^0.5)./(2*cte*r.^(7/2))).^(-1));

%% PERIODO ORBITAL Y ANGULOS EN FUNCION DEL TIEMPO

% Periodo de orbitas
T = 2*pi*sqrt(r.^3/mu);                 % s
anom_ver_punto = 1./sqrt(r.^3/mu);      % Velocidad angular anomalia verdadera

N=1e3+1; % Mallado temporal
for orb = 1:1
    time(orb,:) = linspace(0,T(orb),N);            % Vector de tiempos 1 periodo
    anom_ver(orb,:) = time(orb,:)*anom_ver_punto(orb);   % Anomalia verdadera
    for vel = 1:1
        roll(:,orb,vel) = wrapTo2Pi(time(orb,:)*w(vel));              % Rotacion sat sobre su eje Z
    end
end



%% ECLIPSE
disp('*** Eclipses ***')
eclipse = ones(size(anom_ver)); % Vector señal eclipse booleano (0-1)

% Angulos de Sol y eclipse para cada orbita
for orb=1:length(h)
    
    inclinacion = inc(orb);
    
    Reo = Rx(inclinacion)*Rz(RAAN);         % Tierra - Orbita -> angulo con sol + inclinacion
    rho(orb) = asin(rT./(rT + h(orb)));
    beta_s(orb) = pi/2 - acos((Reo*beta_v')'*[0,0,1]');
    phi(orb) = real(2*acos(cos(rho(orb))/cos(beta_s(orb))));
    
    if phi(orb) ~= 0
        eclipse(orb, anom_ver(orb,:) >= (pi - phi(orb)/2) & anom_ver(orb,:) < (pi + phi(orb)/2) ) = 0;
        t_eclipse(orb) = max(time(eclipse(orb,:) == 0)) - min(time(eclipse(orb,:) == 0));
        anom_ver_eclipse_ini(orb) = min(anom_ver(eclipse(orb,:) == 0));
        anom_ver_eclipse_fin(orb) = max(anom_ver(eclipse(orb,:) == 0));

        disp(['Eclipse para h = ',num2str(h(orb)), ' km'])
        disp(['  ','t_eclipse = ',num2str(round(t_eclipse(orb))),' s'])
        disp(['  ','anom_ver_eclipse_ini = ',num2str(rad2deg(anom_ver_eclipse_ini(orb))),' deg'])
        disp(['  ','anom_ver_eclipse_fin = ',num2str(rad2deg(anom_ver_eclipse_fin(orb))+180),' deg'])
        
    else
        disp(['No hay eclipse para h = ',num2str(h(orb)), 'km'])
    end
    
end

% Plot eclipses
h_plot(fig) = figure(fig);
    hold on
    plot(eclipse(1,:),'DisplayName','450 km')
    legend()
    title('Eclipse')
    fig = fig+1;


    
%% ANGULO PANELES

for orb = 1:length(h)       % Bucle en alturas
    for vel = 1:length(w)   % Bucle en velocidades angulares  
        for p = 1:4         % Bucle en paneles  
            [angulo_panel(:,orb,vel,p), senal_panel(:,orb,vel,p)] = ...
                panel(w(vel), time(orb,:), (p-1)*pi/2);  
        end
    end
end

% h_plot(fig)=figure(fig);
%     hold on
%     for p = 1:1
%         %plot(time(1,:), cos(angulo_panel(:,1,1,p)).*senal_panel(:,1,1,p),...
%         %     'DisplayName',['Panel ' num2str(p)])
%         plot(time(1,:), cos(angulo_panel(:,1,1,p)),'DisplayName',['Panel ' num2str(p)])
%          plot(time(1,:), senal_panel(:,1,1,p),'DisplayName',['Señal Panel ' num2str(p)])
%         % plot(t(1,:), cos_panel(i,:).*( cos(anom_ver(1,:)*2 + pi) + 1 )/2)
%     end
%     legend()
%     title('cos(angulo panel)')
%     fig = fig+1;


%% SIMULACION

for orb = 1:length(h)                   % Bucle en alturas
    
    inclinacion = inc(orb);             % Inclinacion para cada orbita
    
    for vel = 1:length(w)               % Bucle en velocidades angulares        
        for p = 1:4                      % Bucle en paneles                 
            for t = 1:length(time(orb,:))   % Bucle en tiempo
                
                C_tierra_sol = Rz(beta);                        % Sol -> Tierra -> beta
                C_plano_tierra = Rx(inclinacion)*Rz(RAAN);      % Tierra -> plano orbital  
                C_orbita_plano = Rz(anom_ver(orb,t));
                % plano orbital -> orbita 
                C_sat_orbita = eye(3);  
                C_sat_tierra = C_sat_orbita*C_orbita_plano*C_plano_tierra*C_tierra_sol ;
                r_tierra = beta_v;
                r_orbita = C_sat_tierra*r_tierra';
%               potencia_panel(t,p,orb,vel) = G*rend*A*fc*(r_orbita'*[0 0 1]')*eclipse(orb,t);  % Caras con panel: Y Z
                potencia_panel(t,1,orb,vel) = G*rend*A*fc*(r_orbita'*[1/sqrt(2) -1/sqrt(2) 0]'); %Modificar en los plot los ejes
                potencia_panel(t,2,orb,vel) = G*rend*A*fc*(r_orbita'*[-1/sqrt(2) 1/sqrt(2) 0]');
                potencia_panel(t,3,orb,vel) = G*rend*A*fc*(r_orbita'*[0 1 0]');
                potencia_panel(t,4,orb,vel) = G*rend*(30.18 + 6*2.17)/100^2*fc*(r_orbita'*[0 -1 0]');
                potencia_panel(t,5,orb,vel) = G*rend*A*fc*(r_orbita'*[1/sqrt(2) 1/sqrt(2) 0]');
                potencia_panel(t,6,orb,vel) = G*rend*A*fc*(r_orbita'*[-1/sqrt(2) -1/sqrt(2) 0]');
                
            end
            
            % potencia_panel(:,p,orb,vel) = potencia_panel(:,p,orb,vel).*eclipse(orb,:)';  % Caras con panel: Y Z
            
        end
    end
end
% Suma de la contribucion de los paneles
potencia_panel = max(0,potencia_panel); %Hacer 0 cuando las caras están de espaldas al Sol
P_m = sum(potencia_panel,2);

%% Plot potencias
for orb = 1:length(h)
    
    for vel = 1:length(w)
        
        [sup, inf] = envelope(P_m(:,1,orb,vel),300,'peak');
        
        sup = sup'.*eclipse(orb,:); 
        inf = inf'.*eclipse(orb,:);
        P_m(:,1,orb,vel) = P_m(:,1,orb,vel).*eclipse(orb,:)';
        
%         sup(sup == 0) = NaN;
%         inf(inf == 0) = NaN;
        
        h_ = figure(5);
            hold on
            plot(time(1,:),P_m(:,1,orb,vel),'-','LineWidth',1.5,'Color',[150, 150, 150]/255,'DisplayName','Potencia total')
            %plot(rad2deg(anom_ver(1,:)),sup,'-.','LineWidth',2,'Color','k','DisplayName','Envolvente superior')
            %plot(rad2deg(anom_ver(1,:)),inf,'--','LineWidth',2,'Color','k','DisplayName','Envolvente inferior') 
            box on; grid on
%             axis([0,361, 0, 16])
            legend('Interpreter', 'Latex', 'location', 'SouthEast')
            xlh = xlabel('$t$ [s]','Interpreter','latex');
            xlh.Position(1) = xlh.Position(1) + abs(xlh.Position(1) * 0.75);
            ylh = ylabel({'$P$ [W]'},'Interpreter','latex');
            ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1) * 0.7);
            ylh.Position(2) = ylh.Position(2) + abs(ylh.Position(2) * 0.45);
            Save_as_PDF(h_,['Figures/Simulacion_Potencia'],0);
          

%             for p = 1:4
%                 plot(rad2deg(anom_ver(1,:)),potencia_panel(:,p,orb,vel),'DisplayName',['Panel ' num2str(p)])
%                 %plot(time(1,:),potencia_panel(:,p,orb,vel),'DisplayName',['Panel ' num2str(p)])
%             end
%             plot(rad2deg(anom_ver(1,:)),P_m(:,1,orb,vel),'LineWidth',0.1,'DisplayName','Potencia total')
%             % plot(time(1,:),P_m(:,1,orb,vel),'DisplayName','Potencia total')
%             box on
%             legend()
%             title(['Simuacion ',num2str(h(orb)),' ',num2str(w(vel))])
%             xlabel('\nu[deg]')
%             %xlabel('tiempo [s]')
%             ylabel('Potencia')
%             hold off
%             fig = fig+1;
    end
end

% Potencias medias
disp(' *** Potencias medias ***')

for orb = 1:length(h)                   % Bucle en alturas    
    disp(['Potencias medias generadas para h = ',num2str(h(orb)), ' km'])    
    for vel = 1:length(w)               % Bucle en velocidades angulares     
        Potencia_media_generada(orb,vel) = trapz(time(orb,:), P_m(:,1,orb,vel))/(T(orb)-t_eclipse);
        c = 0.025; % Degradacion por año
        tf = 1; % Tiempo mision CubeSat
        I_d = 0.85;
        Potencia_media_generada_BOL(orb,vel) = Potencia_media_generada(orb,vel)*I_d;
        Potencia_media_generada_EOL(orb,vel) = Potencia_media_generada(orb,vel)*(1-c)^tf*I_d;
        disp(['  ','w = ', num2str(w(vel)),' rad/s -> ','Pm = ',num2str(Potencia_media_generada(orb,vel)), ' W'])
        disp(['  ','w = ', num2str(w(vel)),' rad/s -> ','Pm_BOL = ',num2str(Potencia_media_generada_BOL(orb,vel)), ' W'])
        disp(['  ','w = ', num2str(w(vel)),' rad/s -> ','Pm_EOL = ',num2str(Potencia_media_generada_EOL(orb,vel)), ' W'])
    end
end



%% PLOTS PARA INFOMRE

h_ = figure();        % Paneles x
    hold on
    plot(time(1,:),potencia_panel(:,5,1,1), '-', 'LineWidth', 2, 'Color', 'k','DisplayName','Panel Y$^+$')
    plot(time(1,:),potencia_panel(:,6,1,1), '--', 'LineWidth', 2, 'Color', 'k','DisplayName','Panel Y$^-$')
    axis([0,2*pi/w(1)*1.55, 0, 4.3])
    box on; grid on
    legend('Interpreter', 'Latex', 'location', 'SouthEast')
    xlh = xlabel('$t$ [s]','Interpreter','latex');
    xlh.Position(1) = xlh.Position(1) + abs(xlh.Position(1) * 0.75);
    ylh = ylabel({'$P$ [W]'},'Interpreter','latex');
    ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1) * 0.7);
    ylh.Position(2) = ylh.Position(2) + abs(ylh.Position(2) * 0.45);
    Save_as_PDF(h_,"Figures/Simulacion_Y",0);
    hold off

h_ = figure();        % Paneles y
    hold on
    plot(time(1,:),potencia_panel(:,3),'-','LineWidth',2,'Color','k','DisplayName','Panel Z$^+$')
%     plot(time(1,:),potencia_panel(:,1),'-.','LineWidth',2,'Color','k','DisplayName','Panel Zenith')
    plot(time(1,:),potencia_panel(:,4),'--','LineWidth',2,'Color','k','DisplayName','Panel Z$^-$')
    axis([0,2*pi/w(1)*1.75, 0, 4.3])
    box on; grid on
    legend('Interpreter', 'Latex', 'location', 'SouthEast')
    xlh = xlabel('$t$ [s]','Interpreter','latex');
    xlh.Position(1) = xlh.Position(1) + abs(xlh.Position(1) * 0.75);
    ylh = ylabel({'$P$ [W]'},'Interpreter','latex');
    ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1) * 0.7);
    ylh.Position(2) = ylh.Position(2) + abs(ylh.Position(2) * 0.45);
    Save_as_PDF(h_,"Figures/Simulacion_Z",0);
    hold off
    
    h_ = figure();        % Paneles z
    hold on
    plot(time(1,:),potencia_panel(:,1),'-','LineWidth',2,'Color','k','DisplayName','Panel X$^+$')    
    plot(time(1,:),potencia_panel(:,2),'--','LineWidth',2,'Color','k','DisplayName','Panel X$^-$')
    axis([0,2*pi/w(1)*1.75, 0, 4.3])
    box on; grid on
    legend('Interpreter', 'Latex', 'location', 'SouthEast')
    xlh = xlabel('$t$ [s]','Interpreter','latex');
    xlh.Position(1) = xlh.Position(1) + abs(xlh.Position(1) * 0.75);
    ylh = ylabel({'$P$ [W]'},'Interpreter','latex');
    ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1) * 0.7);
    ylh.Position(2) = ylh.Position(2) + abs(ylh.Position(2) * 0.45);
    Save_as_PDF(h_,"Figures/Simulacion_X",0);
    hold off

%% FUNCIONES

% Calculo de angulo sol-panel
function [angulo, senal] = panel(w, t, desfase)

    angulo = acos(cos(w*t(:) + desfase));      %rad
    
    senal = ones(size(angulo));
    senal(angulo>pi/2) = 0;
    %angulo(angulo>pi/2) = pi/2;

    %cos_angulo = cos(angulo);
end


% Matrices de cambio de base
function [Rx] = Rx(angle)

    Rx = [1 0 0;... 
          0 cos(angle) sin(angle);...
          0 -sin(angle) cos(angle)];

end

function [Ry] = Ry(angle)

    Ry = [cos(angle) 0 -sin(angle);...
          0 1 0;...
          sin(angle) 0 cos(angle)];

end

function [Rz] = Rz(angle)

    Rz = [cos(angle) sin(angle) 0;...
        -sin(angle) cos(angle) 0;...
        0 0 1];

end