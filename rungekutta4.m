function [vals] = rungekutta4(y0, t0, t_end, dt, func)
    %rungekutta4: Implementierung des Runge-Kutta-Verfahrens 4ter Ordnung.
    %Input:
    %   initial_v  Startwert
    %   t_end     Zeit bis zu der integriert wird
    %   dt           Zeitschritt
    %   func        zu integrierende Funktion

    steps = ceil((t_end-t0)/dt)+1;      % Anzahl Schritte
    vals = zeros(1, steps);             % berechnete Werte
    vals(1) = y0;                       % übernehme Startwert
    
    for t=1:steps-1                     % iteriere bis Zielzeit
        k1 = dt * func(vals(t), t0+dt*(t-1));
        k2 = dt * func(vals(t)+0.5*k1, t0+dt*(t-0.5));
        k3 = dt * func(vals(t)+0.5*k2, t0+dt*(t-0.5));
        k4 = dt * func(vals(t)+k3, t0+dt*(t));
        vals(t+1) = vals(t) + (k1+2*k2+2*k3+k4)/6;
    end
end