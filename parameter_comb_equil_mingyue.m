function parameter_combs
    tspan = [0 200];
    rH = 0.1;  fH = 0.01;
    rC_vals = 0.1:0.1:0.5;                 % 0.1 to 0.5, step 0.1
    fC_vals = 0.01:0.01:0.05;              % 0.01 to 0.05, step 0.01

    H0_0 = 0.9;
    C0_0 = 0.1;
    H1_0 = 0;
    C1_0 = 0;
    y0 = [H0_0; C0_0; H1_0; C1_0];


    % UI
    fig = uifigure('Name','Logistic Growth System','Position',[100 100 850 480]); % [left bottom width height]
    ax  = uiaxes(fig,'Position',[50 140 600 300]); hold(ax,'on'); grid(ax,'on');
    xlabel(ax,'Time'); ylabel(ax,'Population'); title(ax,'H_0, C_0, H_1, C_1');
    lg = legend(ax,'Location','best');

    % Lines
    hH0 = plot(ax, NaN, NaN, 'b-', 'LineWidth', 2, 'DisplayName','H_0');
    hC0 = plot(ax, NaN, NaN, 'r-', 'LineWidth', 2, 'DisplayName','C_0');
    hH1 = plot(ax, NaN, NaN, 'b-', 'LineWidth', 4, 'DisplayName','H_1');
    hC1 = plot(ax, NaN, NaN, 'r-', 'LineWidth', 4, 'DisplayName','C_1');

    % rC slider 0.1 to 0.5
    uilabel(fig,'Position',[50 95 60 22],'Text','r_C');
    rCslider = uislider(fig,'Position',[100 105 500 3], ...
        'Limits',[min(rC_vals) max(rC_vals)], ...
        'MajorTicks', rC_vals, ...
        'MajorTickLabels', compose('%.2f', rC_vals), ...
        'Value', rC_vals(1));
    rClbl = uilabel(fig,'Position',[620 95 120 22],'Text',sprintf('r_C = %.2f', rC_vals(1))); % shows current value

    % fC slider 0.01 to 0.05
    uilabel(fig,'Position',[50 55 60 22],'Text','f_C');
    fCslider = uislider(fig,'Position',[100 65 500 3], ...
        'Limits',[min(fC_vals) max(fC_vals)], ...
        'MajorTicks', fC_vals, ...
        'MajorTickLabels', compose('%.2f', fC_vals), ...
        'Value', fC_vals(1));
    fClbl = uilabel(fig,'Position',[620 55 120 22],'Text',sprintf('f_C = %.2f', fC_vals(1)));

    % Equilibrium
    eqLbl = uilabel(fig,'Position',[680 300 150 60],'FontSize',11,'Text','Final cancer: --');

    % Callbacks (snap to nearest allowed value)
    rCslider.ValueChangingFcn = @(~,evt) onChange('rC',evt.Value); % evt: current
    fCslider.ValueChangingFcn = @(~,evt) onChange('fC',evt.Value);

    % Initial draw
    recompute();

    % Helper function 1
    function onChange(which,val) % calls when moving a slider
        if strcmp(which,'rC')
            rCslider.Value = snap(val, rC_vals); % snap the slider to the nearest discrete tick
        else
            fCslider.Value = snap(val, fC_vals);
        end
        recompute();
    end

    % Helper 2
    function v = snap(val, allowed) % decide the closest tick
        [~,k] = min(abs(allowed - val));
        v = allowed(k);
    end

    function recompute()
        rC = rCslider.Value; fC = fCslider.Value;
        rClbl.Text = sprintf('r_C = %.2f', rC);
        fClbl.Text = sprintf('f_C = %.2f', fC);

        [t,Y] = ode45(@(t,y) odesys(t,y,rH,rC,fH,fC), tspan, y0);

        set(hH0,'XData',t,'YData',Y(:,1));
        set(hC0,'XData',t,'YData',Y(:,2));
        set(hH1,'XData',t,'YData',Y(:,3));
        set(hC1,'XData',t,'YData',Y(:,4));

        % Display the final populations    
        C0_end = Y(end,2);
        C1_end = Y(end,4);
        Cavg   = 0.5*(C0_end + C1_end);
        Ctot   = C0_end + C1_end;
    
        eqLbl.Text = sprintf([ ...
            'Final cancer:\n' ...
            'C_0:   %.4f\n' ...
            'C_1:   %.4f\n' ...
            'Avg:   %.4f\n' ...
            'Total: %.4f'], C0_end, C1_end, Cavg, Ctot);
    
        drawnow;

    end
end

function dydt = odesys(~, y, rH, rC, fH, fC)
    H0 = y(1); C0 = y(2); H1 = y(3); C1 = y(4);
    dH0dt = rH*H0*(1 - (H0 + C0)) + fH*(H1 - H0);
    dC0dt = rC*C0*(1 - (H0 + C0)) + fC*(C1 - C0);
    dH1dt = rH*H1*(1 - (H1 + C1)) + fH*(H0 - H1);
    dC1dt = rC*C1*(1 - (H1 + C1)) + fC*(C0 - C1);
    dydt = [dH0dt; dC0dt; dH1dt; dC1dt];
end