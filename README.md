# simulating-variant-effect-screens
Simulation of variant effect screens based on yeast surface display for experimental optimization.

## Objectives
- Model error-prone PCR (epPCR) for library generation.
- Model mutant distribution of synthetic phenotypes and / or integrate gold-standard datasets for ground-truth phenotypes.
- Model experimental processes translating variant in yeast to read count from sequencer.
- Simulate experiments with different manipulable experimental variables to optimize conditions.

## Modeling epPCR
`library_analysis.ipynb` and `simulate_epPCR.ipynb`
1. Analyze Pikh-1 library to obtain mutation parameters.
2. Create computational model of epPCR.
3. Validate by comparing model with experimental Pikh-1 library.
4. Support multiple WT templates.

## Model synthetic phenotypes
`assign_phenotypes.ipynb`
1. Define interface residues.
2. Use ESM to:
3. Simulate ddG fold as the sum of the pseudo-log-likelihood of each mutation not at the interface.
4. Simulate ddG bind as the sume of the pseudo-loglikelihood of each mutation at the interface.

## Import gold-standard datasets for ground-truth phenotypes
`setting_up_gb1_data.ipynb`
1. Download GB1 dataset from Otowinowsky et al.
2. Assign counts to mutants based on nt distance.

## Modeling experimental processes
`simulate_experiment.ipynb`
### Biophysics
Use biophysical models of protein stability and binding to connect true ddG values to simulated expression and binding probabilities.

### Sorting
Use statistical models for signal variability and FACS gating to connect biophysics to probability of being sorted.

### Sequencing
Use sampling probabilities to translate into discrete read counts.

## Simulate experiments
`experiment_optimization.ipynb`
1. Create synthetic library with assigned true latent phenotypes and/or use gold-standard datasets.
2. Run through experimental model using baseline parameters.
4. Plot correlation as a function of each parameter.