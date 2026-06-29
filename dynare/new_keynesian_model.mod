/*
 * Baseline New Keynesian Model -- Dynare Implementation
 * Source: Galí (2008), Monetary Policy, Inflation, and the Business Cycle,
 *         Princeton University Press, Chapter 3
 * Author: Aryan Pandey
 *
 * All variables expressed in deviations from steady state.
 *
 * Solves the canonical 3-equation New Keynesian model (NK Phillips curve,
 * dynamic IS curve, smoothed Taylor rule) via Dynare's rational-expectations
 * solver, then runs three comparative experiments:
 *   Experiment 1 -- effect of price stickiness (Calvo theta) on the response
 *                   to a monetary policy shock
 *   Experiment 2 -- effect of interest-rate smoothing (rho_i) on the response
 *                   to a monetary policy shock
 *   Experiment 3 -- effect of price stickiness on the response to a
 *                   technology shock
 */

%----------------------------------------------------------------
%  0. Preprocessor Directives
%----------------------------------------------------------------

@#define money_growth_rule = 0

%----------------------------------------------------------------
%  1. Variable Declaration
%----------------------------------------------------------------

var
    pi          ${\pi}$             (long_name='Inflation')
    y_gap       ${\tilde{y}}$       (long_name='Output Gap')
    y_nat       ${y^{nat}}$         (long_name='Natural Output')
    y           ${y}$               (long_name='Output')
    r_nat       ${r^{nat}}$         (long_name='Natural Interest Rate')
    r_real      ${r^{r}}$           (long_name='Real Interest Rate')
    i           ${i}$               (long_name='Nominal Interest Rate')
    n           ${n}$               (long_name='Hours Worked')
    m_real      ${m-p}$             (long_name='Real Money Stock')
    m_growth_ann ${\Delta m}$       (long_name='Money Growth Annualized')
    a           ${a}$               (long_name='AR(1) Technology Shock')
    r_real_ann  ${r^{r,ann}}$       (long_name='Annualized Real Interest Rate')
    i_ann       ${i^{ann}}$         (long_name='Annualized Nominal Interest Rate')
    r_nat_ann   ${r^{nat,ann}}$     (long_name='Annualized Natural Interest Rate')
    pi_ann      ${\pi^{ann}}$       (long_name='Annualized Inflation Rate')
    @#if money_growth_rule == 0
        nu      ${\nu}$             (long_name='AR(1) Monetary Policy Shock')
    @#else
        money_growth ${\Delta m_q}$ (long_name='Money Growth')
    @#endif
    ;

%----------------------------------------------------------------
%  2. Shock Declaration
%----------------------------------------------------------------

varexo
    eps_a   ${\varepsilon_a}$       (long_name='Technology Shock')
    @#if money_growth_rule == 0
        eps_nu  ${\varepsilon_\nu}$ (long_name='Monetary Policy Shock')
    @#else
        eps_m   ${\varepsilon_m}$   (long_name='Money Growth Rate Shock')
    @#endif
    ;

%----------------------------------------------------------------
%  3. Parameter Declaration
%----------------------------------------------------------------

parameters
    alppha      ${\alpha}$          (long_name='Capital Share')
    betta       ${\beta}$           (long_name='Discount Factor')
    siggma      ${\sigma}$          (long_name='Inverse EIS (Log Utility)')
    phi         ${\phi}$            (long_name='Inverse Frisch Elasticity')
    phi_pi      ${\phi_{\pi}}$      (long_name='Inflation Feedback  -  Taylor Rule')
    phi_y       ${\phi_{y}}$        (long_name='Output Feedback  -  Taylor Rule')
    eta         ${\eta}$            (long_name='Semi-Elasticity of Money Demand')
    epsilon     ${\epsilon}$        (long_name='Demand Elasticity')
    theta       ${\theta}$          (long_name='Calvo Parameter')
    rho_a       ${\rho_a}$          (long_name='Autocorrelation  -  Technology Shock')
    rho_i       ${\rho_i}$          (long_name='Interest Rate Smoothing')
    @#if money_growth_rule == 0
        rho_nu  ${\rho_{\nu}}$      (long_name='Autocorrelation  -  Monetary Policy Shock')
    @#else
        rho_m   ${\rho_{m}}$        (long_name='Autocorrelation  -  Money Growth Shock')
    @#endif
    ;

%----------------------------------------------------------------
%  4. Calibration (Galí 2008, p. 52)
%----------------------------------------------------------------

siggma  = 1;
phi     = 1;
phi_pi  = 1.5;
phi_y   = 0.5/4;
eta     = 4;
epsilon = 6;
alppha  = 1/3;
betta   = 0.99;
rho_a   = 0.9;
theta   = 0.66;     % Baseline (overwritten in Experiment 1 loop)
rho_i   = 0;        % Baseline (overwritten in Experiment 2 loop)

@#if money_growth_rule == 0
    rho_nu = 0.5;
@#else
    rho_m  = 0.5;
@#endif

%----------------------------------------------------------------
%  5. Model Equations
%----------------------------------------------------------------

model(linear);

    // Composite Parameters
    #Omega      = (1 - alppha) / (1 - alppha + alppha*epsilon);             // p. 47
    #psi_n_ya   = (1 + phi) / (siggma*(1 - alppha) + phi + alppha);        // p. 48
    #lambda     = (1 - theta)*(1 - betta*theta) / theta * Omega;           // p. 47
    #kappa      = lambda * (siggma + (phi + alppha)/(1 - alppha));          // p. 49

    // 1. New Keynesian Phillips Curve [eq. 21]
    pi = betta*pi(+1) + kappa*y_gap;

    // 2. Dynamic IS Curve [eq. 22]
    y_gap = -1/siggma * (i - pi(+1) - r_nat) + y_gap(+1);

    // 3. Monetary Policy Rule [eq. 25]  -  Smoothed Taylor Rule
    @#if money_growth_rule == 0
        i = rho_i*i(-1) + (1 - rho_i)*(phi_pi*pi + phi_y*y_gap) + nu;
    @#endif

    // 4. Natural Rate of Interest [eq. 23]
    r_nat = siggma * psi_n_ya * (a(+1) - a);

    // 5. Real Interest Rate Definition
    r_real = i - pi(+1);

    // 6. Natural Output [eq. 19]
    y_nat = psi_n_ya * a;

    // 7. Output Gap Definition
    y_gap = y - y_nat;

    // 8. Shock Processes
    @#if money_growth_rule == 0
        nu = rho_nu*nu(-1) + eps_nu;
    @#endif
    a = rho_a*a(-1) + eps_a;

    // 9. Production Function [eq. 13]
    y = a + (1 - alppha)*n;

    // 10. Money Growth [derived from eq. 4]
    m_growth_ann = 4*(y - y(-1) - eta*(i - i(-1)) + pi);

    // 11. Real Money Demand [eq. 4]
    m_real = y - eta*i;

    @#if money_growth_rule == 1
        money_growth = m_real - m_real(-1) + pi;
        money_growth = rho_m*money_growth(-1) + eps_m;
    @#endif

    // 12. Annualized Variables
    i_ann       = 4*i;
    r_real_ann  = 4*r_real;
    r_nat_ann   = 4*r_nat;
    pi_ann      = 4*pi;

end;

%----------------------------------------------------------------
%  6. Shock Variances  -  Monetary Policy Shock ON, Technology OFF
%----------------------------------------------------------------

shocks;
    @#if money_growth_rule == 0
        var eps_nu = 0.25^2;    % Monetary policy shock ON
    @#else
        var eps_m  = 0.25^2;
    @#endif
    var eps_a = 0;              % Technology shock OFF
end;

%----------------------------------------------------------------
%  7. Steady State & Stability Check
%----------------------------------------------------------------

resid;
steady;
check;

%----------------------------------------------------------------
%  8. Initialize stoch_simul
%----------------------------------------------------------------

stoch_simul(order=1, irf=15, nograph) y_gap pi_ann i_ann r_real_ann m_growth_ann nu;

%================================================================
%  EXPERIMENT 1: IRFs Across Different Price Stickiness Levels (theta)
%================================================================

theta_vec    = [0.2, 0.5, 0.66, 0.75];
theta_labels = {'\theta = 0.2', '\theta = 0.5', '\theta = 0.66 (Baseline)', '\theta = 0.75'};

% Reset rho_i to baseline for Experiment 1
set_param_value('rho_i', 0);

% Pre-allocate storage
IRF_y_gap_th        = zeros(15, length(theta_vec));
IRF_pi_ann_th       = zeros(15, length(theta_vec));
IRF_i_ann_th        = zeros(15, length(theta_vec));
IRF_r_real_ann_th   = zeros(15, length(theta_vec));
IRF_m_growth_ann_th = zeros(15, length(theta_vec));
IRF_nu_th           = zeros(15, length(theta_vec));

% Loop over theta values
for j = 1:length(theta_vec)
    set_param_value('theta', theta_vec(j));
    [info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
    IRF_y_gap_th(:, j)        = oo_.irfs.y_gap_eps_nu;
    IRF_pi_ann_th(:, j)       = oo_.irfs.pi_ann_eps_nu;
    IRF_i_ann_th(:, j)        = oo_.irfs.i_ann_eps_nu;
    IRF_r_real_ann_th(:, j)   = oo_.irfs.r_real_ann_eps_nu;
    IRF_m_growth_ann_th(:, j) = oo_.irfs.m_growth_ann_eps_nu;
    IRF_nu_th(:, j)           = oo_.irfs.nu_eps_nu;
end

% Reset theta to baseline after Experiment 1 loop
set_param_value('theta', 0.66);

% Plot Experiment 1
newcolors = [0.00 0.45 0.74;
             0.47 0.67 0.19;
             0.85 0.33 0.10;
             0.93 0.69 0.13];
t = 1:15;

figure('Name', 'Experiment 1: Monetary Policy Shock  -  Different Theta', 'Position', [100, 100, 1000, 800]);
sgtitle('Monetary Policy Shock: Different \theta', 'FontSize', 14, 'FontWeight', 'bold');

subplot(3,2,1); hold on;
for k = 1:4; plot(t, IRF_y_gap_th(:,k), '-', 'Color', newcolors(k,:), 'LineWidth', 2); end
title('Output Gap (y_t^{gap})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(3,2,2); hold on;
for k = 1:4; plot(t, IRF_pi_ann_th(:,k), '-', 'Color', newcolors(k,:), 'LineWidth', 2); end
title('Inflation (\pi_{ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(3,2,3); hold on;
for k = 1:4; plot(t, IRF_i_ann_th(:,k), '-', 'Color', newcolors(k,:), 'LineWidth', 2); end
title('Nominal Interest Rate (i_{ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(3,2,4); hold on;
for k = 1:4; plot(t, IRF_r_real_ann_th(:,k), '-', 'Color', newcolors(k,:), 'LineWidth', 2); end
title('Real Interest Rate (r_{real,ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(3,2,5); hold on;
for k = 1:4; plot(t, IRF_m_growth_ann_th(:,k), '-', 'Color', newcolors(k,:), 'LineWidth', 2); end
title('Money Growth (m\_growth_{ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

s6 = subplot(3,2,6); hold on;
handles = gobjects(4,1);
for k = 1:4; handles(k) = plot(t, IRF_nu_th(:,k), '-', 'Color', newcolors(k,:), 'LineWidth', 2); end
title('Monetary Shock (\nu)'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);
legend(s6, handles, theta_labels, 'Location', 'southoutside', 'Orientation', 'horizontal');

%================================================================
%  EXPERIMENT 1 (continued): Monetary Policy Shock IRF  -  Baseline Only (theta = 2/3, rho_i = 0)
%================================================================

% Reset parameters to baseline
set_param_value('theta', 2/3);
set_param_value('rho_i', 0);

% Shocks: Monetary Policy ON, Technology OFF
shocks;
    @#if money_growth_rule == 0
        var eps_nu = 0.25^2;
    @#else
        var eps_m  = 0.25^2;
    @#endif
    var eps_a = 0;
end;

resid;
steady;
check;

stoch_simul(order=1, irf=15, nograph) y_gap pi_ann i_ann r_real_ann m_growth_ann nu;

base_color = [0.85 0.33 0.10];   % Red

figure('Name', 'Experiment 1: Monetary Policy Shock  -  Baseline (theta = 2/3)', 'Position', [100, 100, 1000, 800]);
sgtitle('Monetary Policy Shock: Baseline \theta = 2/3', 'FontSize', 14, 'FontWeight', 'bold');

subplot(3,2,1);
plot(t, oo_.irfs.y_gap_eps_nu, '-', 'Color', base_color, 'LineWidth', 2);
title('Output Gap (y_{gap})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(3,2,2);
plot(t, oo_.irfs.pi_ann_eps_nu, '-', 'Color', base_color, 'LineWidth', 2);
title('Inflation (\pi_{ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(3,2,3);
plot(t, oo_.irfs.i_ann_eps_nu, '-', 'Color', base_color, 'LineWidth', 2);
title('Nominal Interest Rate (i_{ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(3,2,4);
plot(t, oo_.irfs.r_real_ann_eps_nu, '-', 'Color', base_color, 'LineWidth', 2);
title('Real Interest Rate (r_{real,ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(3,2,5);
plot(t, oo_.irfs.m_growth_ann_eps_nu, '-', 'Color', base_color, 'LineWidth', 2);
title('Money Growth (m\_growth_{ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(3,2,6);
plot(t, oo_.irfs.nu_eps_nu, '-', 'Color', base_color, 'LineWidth', 2);
title('Monetary Shock (\nu)'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

%================================================================
%  EXPERIMENT 2: IRFs Across Different Interest Rate Smoothing Levels (rho_i)
%================================================================

rho_vec    = [0, 0.5, 0.8, 0.95];
rho_labels = {'\rho_i = 0 (Baseline)', '\rho_i = 0.5', '\rho_i = 0.8', '\rho_i = 0.95'};

% Pre-allocate storage
IRF_y_gap_rh          = zeros(15, length(rho_vec));
IRF_pi_ann_rh         = zeros(15, length(rho_vec));
IRF_i_ann_rh          = zeros(15, length(rho_vec));
IRF_r_real_ann_rh     = zeros(15, length(rho_vec));
IRF_m_growth_ann_rh   = zeros(15, length(rho_vec));
IRF_nu_rh             = zeros(15, length(rho_vec));

% Loop over rho_i values
for j = 1:length(rho_vec)
    set_param_value('rho_i', rho_vec(j));
    [info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
    IRF_y_gap_rh(:, j)        = oo_.irfs.y_gap_eps_nu;
    IRF_pi_ann_rh(:, j)       = oo_.irfs.pi_ann_eps_nu;
    IRF_i_ann_rh(:, j)        = oo_.irfs.i_ann_eps_nu;
    IRF_r_real_ann_rh(:, j)   = oo_.irfs.r_real_ann_eps_nu;
    IRF_m_growth_ann_rh(:, j) = oo_.irfs.m_growth_ann_eps_nu;
    IRF_nu_rh(:, j)           = oo_.irfs.nu_eps_nu;
end

% Plot Experiment 2
figure('Name', 'Experiment 2: Monetary Policy Shock  -  Different Rho_i', 'Position', [100, 100, 1000, 800]);
sgtitle('Monetary Policy Shock: Different \rho_i', 'FontSize', 14, 'FontWeight', 'bold');

subplot(3,2,1); hold on;
for k = 1:4; plot(t, IRF_y_gap_rh(:,k), '-', 'Color', newcolors(k,:), 'LineWidth', 2); end
title('Output Gap (y_t^{gap})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(3,2,2); hold on;
for k = 1:4; plot(t, IRF_pi_ann_rh(:,k), '-', 'Color', newcolors(k,:), 'LineWidth', 2); end
title('Inflation (\pi_{ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(3,2,3); hold on;
for k = 1:4; plot(t, IRF_i_ann_rh(:,k), '-', 'Color', newcolors(k,:), 'LineWidth', 2); end
title('Nominal Interest Rate (i_{ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(3,2,4); hold on;
for k = 1:4; plot(t, IRF_r_real_ann_rh(:,k), '-', 'Color', newcolors(k,:), 'LineWidth', 2); end
title('Real Interest Rate (r_{real,ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(3,2,5); hold on;
for k = 1:4; plot(t, IRF_m_growth_ann_rh(:,k), '-', 'Color', newcolors(k,:), 'LineWidth', 2); end
title('Money Growth (m\_growth_{ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

s6 = subplot(3,2,6); hold on;
handles = gobjects(4,1);
for k = 1:4; handles(k) = plot(t, IRF_nu_rh(:,k), '-', 'Color', newcolors(k,:), 'LineWidth', 2); end
title('Monetary Shock (\nu)'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);
legend(s6, handles, rho_labels, 'Location', 'southoutside', 'Orientation', 'horizontal');

%================================================================
%  EXPERIMENT 3: IRFs for Technology Shock  -  theta = 2/3 (Baseline) and theta = 0.2
%================================================================

% Reset rho_i to baseline
set_param_value('rho_i', 0);

% Update shocks: Technology ON, Monetary Policy OFF
shocks;
    @#if money_growth_rule == 0
        var eps_nu = 0;
    @#else
        var eps_m  = 0;
    @#endif
    var eps_a = 1^2;
end;

resid;
steady;
check;

stoch_simul(order=1, irf=15, nograph) y_gap pi_ann y n i_ann r_real_ann m_growth_ann a;

theta_tech        = [2/3, 0.2];
theta_labels_tech = {'\theta = 2/3 (Baseline)', '\theta = 0.2 (Flexible)'};

% Pre-allocate storage
IRF_y_gap_tech        = zeros(15, length(theta_tech));
IRF_pi_ann_tech       = zeros(15, length(theta_tech));
IRF_y_tech            = zeros(15, length(theta_tech));
IRF_n_tech            = zeros(15, length(theta_tech));
IRF_i_ann_tech        = zeros(15, length(theta_tech));
IRF_r_real_ann_tech   = zeros(15, length(theta_tech));
IRF_m_growth_ann_tech = zeros(15, length(theta_tech));
IRF_a_tech            = zeros(15, length(theta_tech));

% Loop over theta values
for j = 1:length(theta_tech)
    set_param_value('theta', theta_tech(j));
    [info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
    IRF_y_gap_tech(:, j)        = oo_.irfs.y_gap_eps_a;
    IRF_pi_ann_tech(:, j)       = oo_.irfs.pi_ann_eps_a;
    IRF_y_tech(:, j)            = oo_.irfs.y_eps_a;
    IRF_n_tech(:, j)            = oo_.irfs.n_eps_a;
    IRF_i_ann_tech(:, j)        = oo_.irfs.i_ann_eps_a;
    IRF_r_real_ann_tech(:, j)   = oo_.irfs.r_real_ann_eps_a;
    IRF_m_growth_ann_tech(:, j) = oo_.irfs.m_growth_ann_eps_a;
    IRF_a_tech(:, j)            = oo_.irfs.a_eps_a;
end

tech_colors = [0.00 0.45 0.74;
               0.85 0.33 0.10];

figure('Name', 'Experiment 3: Technology Shock  -  Baseline vs Flexible Prices', 'Position', [100, 100, 1000, 800]);
sgtitle('Technology Shock: \theta = 2/3 (Baseline) vs \theta = 0.2 (Flexible)', ...
        'FontSize', 14, 'FontWeight', 'bold');

subplot(4,2,1); hold on;
for k = 1:2; plot(t, IRF_y_gap_tech(:,k), '-', 'Color', tech_colors(k,:), 'LineWidth', 2); end
title('Output Gap (y_{gap})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(4,2,2); hold on;
for k = 1:2; plot(t, IRF_pi_ann_tech(:,k), '-', 'Color', tech_colors(k,:), 'LineWidth', 2); end
title('Inflation (\pi_{ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(4,2,3); hold on;
for k = 1:2; plot(t, IRF_y_tech(:,k), '-', 'Color', tech_colors(k,:), 'LineWidth', 2); end
title('Output (y)'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(4,2,4); hold on;
for k = 1:2; plot(t, IRF_n_tech(:,k), '-', 'Color', tech_colors(k,:), 'LineWidth', 2); end
title('Hours Worked (n)'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(4,2,5); hold on;
for k = 1:2; plot(t, IRF_i_ann_tech(:,k), '-', 'Color', tech_colors(k,:), 'LineWidth', 2); end
title('Nominal Interest Rate (i_{ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(4,2,6); hold on;
for k = 1:2; plot(t, IRF_r_real_ann_tech(:,k), '-', 'Color', tech_colors(k,:), 'LineWidth', 2); end
title('Real Interest Rate (r_{real,ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(4,2,7); hold on;
for k = 1:2; plot(t, IRF_m_growth_ann_tech(:,k), '-', 'Color', tech_colors(k,:), 'LineWidth', 2); end
title('Money Growth (m\_growth_{ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

s8 = subplot(4,2,8); hold on;
handles = gobjects(2,1);
for k = 1:2; handles(k) = plot(t, IRF_a_tech(:,k), '-', 'Color', tech_colors(k,:), 'LineWidth', 2); end
title('Technology (a)'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);
legend(s8, handles, theta_labels_tech, 'Location', 'southoutside', 'Orientation', 'horizontal');

%================================================================
%  EXPERIMENT 3 (continued): Technology Shock IRF  -  Baseline Only (theta = 2/3)
%================================================================

% Reset parameters to baseline
set_param_value('theta', 2/3);
set_param_value('rho_i', 0);

% Shocks: Technology ON, Monetary Policy OFF
shocks;
    @#if money_growth_rule == 0
        var eps_nu = 0;
    @#else
        var eps_m  = 0;
    @#endif
    var eps_a = 1^2;
end;

resid;
steady;
check;

stoch_simul(order=1, irf=15, nograph) y_gap pi_ann y n i_ann r_real_ann m_growth_ann a;

base_color = [0.85 0.33 0.10];   % Red

figure('Name', 'Experiment 3: Technology Shock  -  Baseline Only', 'Position', [100, 100, 1000, 800]);
sgtitle('Technology Shock: Baseline \theta = 2/3', 'FontSize', 14, 'FontWeight', 'bold');

subplot(4,2,1);
plot(t, oo_.irfs.y_gap_eps_a, '-', 'Color', base_color, 'LineWidth', 2);
title('Output Gap (y_{gap})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(4,2,2);
plot(t, oo_.irfs.pi_ann_eps_a, '-', 'Color', base_color, 'LineWidth', 2);
title('Inflation (\pi_{ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(4,2,3);
plot(t, oo_.irfs.y_eps_a, '-', 'Color', base_color, 'LineWidth', 2);
title('Output (y)'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(4,2,4);
plot(t, oo_.irfs.n_eps_a, '-', 'Color', base_color, 'LineWidth', 2);
title('Hours Worked (n)'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(4,2,5);
plot(t, oo_.irfs.i_ann_eps_a, '-', 'Color', base_color, 'LineWidth', 2);
title('Nominal Interest Rate (i_{ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(4,2,6);
plot(t, oo_.irfs.r_real_ann_eps_a, '-', 'Color', base_color, 'LineWidth', 2);
title('Real Interest Rate (r_{real,ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(4,2,7);
plot(t, oo_.irfs.m_growth_ann_eps_a, '-', 'Color', base_color, 'LineWidth', 2);
title('Money Growth (m\_growth_{ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(4,2,8);
plot(t, oo_.irfs.a_eps_a, '-', 'Color', base_color, 'LineWidth', 2);
title('Technology (a)'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

%================================================================
%  EXPERIMENT 2 (continued): Monetary Policy Shock  -  rho_i = 0 vs rho_i = 0.8 (Baseline theta = 2/3)
%================================================================

% Reset theta to baseline
set_param_value('theta', 2/3);

% Shocks: Monetary Policy ON, Technology OFF
shocks;
    @#if money_growth_rule == 0
        var eps_nu = 0.25^2;
    @#else
        var eps_m  = 0.25^2;
    @#endif
    var eps_a = 0;
end;

resid;
steady;
check;

stoch_simul(order=1, irf=15, nograph) y_gap pi_ann i_ann r_real_ann m_growth_ann nu;

rho_compare     = [0, 0.8];
rho_labels_comp = {'\rho_i = 0 (Baseline)', '\rho_i = 0.8'};

% Pre-allocate storage
IRF_y_gap_comp      = zeros(15, length(rho_compare));
IRF_pi_ann_comp     = zeros(15, length(rho_compare));
IRF_i_ann_comp      = zeros(15, length(rho_compare));
IRF_r_real_ann_comp = zeros(15, length(rho_compare));
IRF_m_growth_comp   = zeros(15, length(rho_compare));
IRF_nu_comp         = zeros(15, length(rho_compare));

% Loop over rho_i values
for j = 1:length(rho_compare)
    set_param_value('rho_i', rho_compare(j));
    [info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
    IRF_y_gap_comp(:, j)      = oo_.irfs.y_gap_eps_nu;
    IRF_pi_ann_comp(:, j)     = oo_.irfs.pi_ann_eps_nu;
    IRF_i_ann_comp(:, j)      = oo_.irfs.i_ann_eps_nu;
    IRF_r_real_ann_comp(:, j) = oo_.irfs.r_real_ann_eps_nu;
    IRF_m_growth_comp(:, j)   = oo_.irfs.m_growth_ann_eps_nu;
    IRF_nu_comp(:, j)         = oo_.irfs.nu_eps_nu;
end

comp_colors = [0.00 0.45 0.74;
               0.85 0.33 0.10];

figure('Name', 'Experiment 2: Monetary Policy Shock  -  rho_i = 0 vs rho_i = 0.8', 'Position', [100, 100, 1000, 800]);
sgtitle('Monetary Policy Shock: \rho_i = 0 vs \rho_i = 0.8', 'FontSize', 14, 'FontWeight', 'bold');

subplot(3,2,1); hold on;
for k = 1:2; plot(t, IRF_y_gap_comp(:,k), '-', 'Color', comp_colors(k,:), 'LineWidth', 2); end
title('Output Gap (y_{gap})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(3,2,2); hold on;
for k = 1:2; plot(t, IRF_pi_ann_comp(:,k), '-', 'Color', comp_colors(k,:), 'LineWidth', 2); end
title('Inflation (\pi_{ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(3,2,3); hold on;
for k = 1:2; plot(t, IRF_i_ann_comp(:,k), '-', 'Color', comp_colors(k,:), 'LineWidth', 2); end
title('Nominal Interest Rate (i_{ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(3,2,4); hold on;
for k = 1:2; plot(t, IRF_r_real_ann_comp(:,k), '-', 'Color', comp_colors(k,:), 'LineWidth', 2); end
title('Real Interest Rate (r_{real,ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

subplot(3,2,5); hold on;
for k = 1:2; plot(t, IRF_m_growth_comp(:,k), '-', 'Color', comp_colors(k,:), 'LineWidth', 2); end
title('Money Growth (m\_growth_{ann})'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);

s6 = subplot(3,2,6); hold on;
handles = gobjects(2,1);
for k = 1:2; handles(k) = plot(t, IRF_nu_comp(:,k), '-', 'Color', comp_colors(k,:), 'LineWidth', 2); end
title('Monetary Shock (\nu)'); xlabel('Quarters'); ylabel('Deviation');
xlim([1 15]); yline(0,'k-','HandleVisibility','off'); grid on; set(gca, 'GridAlpha', 0.4, 'GridColor', [0.3 0.3 0.3]);
legend(s6, handles, rho_labels_comp, 'Location', 'southoutside', 'Orientation', 'horizontal');

write_latex_dynamic_model;

