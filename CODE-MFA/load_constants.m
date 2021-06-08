% constants for thermodynamics
RT = 2.5;

% constants for cytosol/mitochondria compartments
CY_WC_VOLUME=0.8;
MT_WC_VOLUME=1-CY_WC_VOLUME;

% epsilon value
EPSILON_VALUE = 1e-5;    

% chi square distribution  - one degree of freedom
CONSTANT_VALUE_FOR_CONFIDENCE_INTERVAL = 3.84;

% minimum Gibbs free energy distance from zero
dG_MIN_DISTANCE_FROM_ZERO = 1;

% concentration STD factor - standard deviation of the lowest concentration
% among the measured or the simulated. This way, the simulated WC
% concentration will try to be close to the measured WC concentration
WC_CONCENTRATION_STD_FACTOR = 0.2;

% constant for figures
RED_COLOR           = [0.8500, 0.3250, 0.0980];
GREEN_COLOR         = [0.4660, 0.6740, 0.1880];
GREEN_DARK_COLOR    = [0.31, 0.4, 0.1880];
BLUE_COLOR          = [0, 0.4470, 0.7410];
BLUE_DARK_COLOR     = [0, 0.3, 0.5];
GRAY_COLOR          = [0.8 0.8 0.8];
GRAY_DARK_COLOR     = [0.4 0.4 0.4];
YELLOW_COLOR        = [0.9373 0.6902 0.1294];

