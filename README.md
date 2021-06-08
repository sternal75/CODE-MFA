# CODE-MFA
Inferring mitochondrial and cytosolic metabolism by coupling isotope tracing and deconvolution

The inability to inspect metabolic activities within distinct subcellular compartments has been a major barrier to our understanding of eukaryotic cell metabolism. Previous work addressed this challenge by analyzing metabolism in isolated organelles, which grossly bias metabolic activity. Here, we developed a method for inferring physiological metabolic fluxes and metabolite concentrations in mitochondria and cytosol based on isotope tracing experiments performed with intact cells. This is made possible by computational deconvolution of metabolite isotopic labeling patterns and concentrations into cytosolic and mitochondrial counterparts, coupled with metabolic and thermodynamic modelling. Our approach lowers the uncertainty regarding compartmentalized fluxes and concentrations by one and three orders of magnitude compared to existing modelling approaches, respectively. We derive a first quantitative view of mitochondrial and cytosolic metabolic activities in central carbon metabolism across cultured cell lines, finding major variability in compartmentalized malate-aspartate shuttle fluxes. We expect our approach for directly probing metabolism at a subcellular resolution to be instrumental for a variety of studies of metabolic dysfunction in human disease and for bioengineering.

The code is written in Matlab.
To run the code use Matlab 2014 or later version.

Pre-running installation requirements:
Install COBRA toolbox, available at the following link https://opencobra.github.io/cobratoolbox/stable/index.html

Inputs to run CODE-MFA method:
Model: 
xls_input_files/input.xlsx: reactions with lower and upper bounds for fluxes, carbon mapping, list of metabolites including extracellular metabolites with their labeled atoms.
xls_input_files/input_for_thermodynamics.xlsx: reactions with their Gibbs free energies, and list of metabolites with lower and upper bounds measured concentrations.
Isotope tracing of metabolites:
processIsotopicLabel/data: Mass Isotopomer Distribution (MID) of metabolites while feeding cells with Glucose, Glutamine, and Lactate.

Running CODE-MFA:
run calc.m.
Note that running CODE-MFA algorithm is computationaly hard and may take a few days, thus using parallel computing is mandatory.

Plotting resulting figures:
After running the CODE-MFA, run plot/output_all_figures_for_paper.m. All figures will be saved in plot/output_images directory in pdf format.

