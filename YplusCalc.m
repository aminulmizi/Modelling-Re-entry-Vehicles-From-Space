function starCCM_PrismCalculator
% STAR-CCM+ Prism Layer Configuration Tool
% Calculates first cell height, prism layer parameters, and boundary layer estimates

    % Input parameters
    p = input('Enter fluid density (kg/m³): ');
    mu = input('Enter dynamic viscosity (Pa·s): ');
    U_inf = input('Enter freestream velocity (m/s): ');
    L = input('Enter characteristic length (m): ');
    targetYplus = input('Enter target Y+ value: ');
    Gr = input('Enter growth ratio (recommended intial value of 1.2): ');
    Nl = input('Enter number of layers (recommended intial value of 25): ');

    % Reynolds number calculation
    Re = p * U_inf * L / mu;
    
    % Skin friction coefficient (Schlichting for turbulent flow)
    Cf = (2*log10(Re) - 0.65).^(-2.3);  % Valid for Re < 1e9 [4]
    
    % Friction velocity calculation
    u_tau = U_inf * sqrt(Cf/2);
    
    % First cell height calculation [6][8]
    y1 = (targetYplus * mu) / (p * u_tau);
    
    % Boundary layer thickness estimation (flat plate turbulent) [7]
    delta = 0.37 * L / (Re^(1/5));
    
    % Prism layer recommendations
    num_layers = Nl;  % Recommended starting point [2][6]
    growth_ratio = Gr;  % Conservative growth rate [2][6]
    
    % Total prism thickness calculation
    total_thickness = y1 * (growth_ratio^num_layers - 1)/(growth_ratio - 1);
    
    % STAR-CCM+ specific validation
    min_wall_distance = 1e-9;  % Default reference value [6]
    if y1 < min_wall_distance
        warning(['First cell height below minimum wall distance! '...
            'Adjust Physics > Reference Values > Minimum Allowable Wall Distance']);
    end
    
    % Output display
    fprintf('\n--- STAR-CCM+ Prism Layer Parameters ---\n');
    fprintf('Reynolds Number: %.4f m\n\n', Re);
    fprintf('First Cell Height (Near Wall Thickness): %.10e m\n', y1);
    fprintf('Recommended Layers: %d (25-30 optimal)\n', num_layers);
    fprintf('Growth Ratio: %.2f (≤1.3 recommended)\n', growth_ratio);
    fprintf('Estimated BL Thickness: %.4f m\n', delta);
    fprintf('Prism Layer Total Thickness: %.10f m\n\n', total_thickness);
    
    
    % Validation checks
    if total_thickness < delta
        warning('Total prism thickness < boundary layer - increase layers/growth ratio');
    end
    
    if growth_ratio > 1.3
        warning('Growth ratio exceeds recommended maximum - reduce for mesh quality');
    end
end