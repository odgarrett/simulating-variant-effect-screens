# simulating-variant-effect-screens
Simulation of variant effect screens based on yeast surface display and yeast two hybrid for experimental optimization.

## Objectives
- Model error-prone PCR (epPCR) for library generation.
- Model mutant distribution of synthetic phenotypes
- Model experimental processes translating variant in yeast to read count from sequencer.
- Simulate experiments with different manipulable experimental variables to optimize conditions.

## Modeling epPCR
1. Analyze Pikh-1 library to obtain mutation parameters.
2. Create computational model of epPCR.
3. Validate by comparing model with experimental Pikh-1 library.
4. Support multiple WT templates.

## Model synthetic phenotypes
1. Define interface residues.
2. Use ESM to:
3. Simulate ddG fold as the sum of the pseudo-log-likelihood of each mutation not at the interface.
4. Simulate ddG bind as the sume of the pseudo-loglikelihood of each mutation at the interface.

## Modeling experimental processes
### Biophysics
Connect the biophysical phenotypes to expression and binding.

### Sorting
Connect expression and binding to probability of being sorted.

### Sequencing
Connect probability of being sorted to read counts.

## Simulate experiments
1. Create synthetic library with assigned true latent phenotypes.
2. Run through experimental model using baseline parameters.
3. Quantify correlation of different scoring methods with true data.
4. Plot correlation as a function of each parameter.
5. Use Bayesian optimization algorithm to find best experimental conditions (best correlation).