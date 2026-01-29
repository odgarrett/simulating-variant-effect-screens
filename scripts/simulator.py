import numpy as np
from scipy.stats import norm



def scale_library_to_physical_count(df, num_cells_sorted):
        current_total = df['count'].sum()
        scaling_factor = num_cells_sorted / current_total
        
        # Scale and convert to integer
        df['count'] = (df['count'] * scaling_factor).astype(int)
        
        return df


def simulate_biophysics(
        df,
        RT,
        wt_dG_fold,
        wt_dG_bind,
        ligand_conc,
        hill_coeff
):
    # Calculate aboslute energies from wild-type and ddG values
    df['dG_fold'] = wt_dG_fold + df['ddG_fold']
    df['p_fold'] = 1.0 / (1.0 + np.exp(df['dG_fold'] / RT))

    df['dG_bind'] = wt_dG_bind + df['ddG_bind']
    df['Kd_variant'] = np.exp(df['dG_bind'] / RT)
    df['occupancy'] = (ligand_conc)**hill_coeff / ((ligand_conc)**hill_coeff + (df['Kd_variant'])**hill_coeff)
    return df


def simulate_facs(
        df,
        t_exp,
        t_bind,
        facs_noise,
):
    # Determine effective threshold
    with np.errstate(divide='ignore'): # handle division by 0
        req_expression_for_binding = t_bind / df['occupancy']
    
    t_effective = np.maximum(t_exp, req_expression_for_binding)

    # Log transform to normal distribution
    mu_log = np.log(df['p_fold'] + 1e-12)
    thresh_log = np.log(t_effective + 1e-12)

    # Calculate negative z-scores
    z_scores = -(thresh_log - mu_log) / facs_noise

    # Use CDF to calculate probability of sorting
    df['p_sort'] = norm.cdf(z_scores)

    # Handle when occupancy is 0, which produces infinite t_effective
    df.loc[np.isinf(t_effective), 'p_sort'] = 0.0

    # Execute the sort
    df['n_cells_sorted'] = np.random.binomial(
        n=df['count'].values.astype(int),
        p=df['p_sort'].values
    )

    return df


def simulate_template_extraction(
        df,
        miniprep_input_cells,
        plasmid_copy_number,
        miniprep_efficiency,
        elution_vol,
        template_vol,
        source_col
):
    # Calculate variant frequencies in outgrowth culture
    total_sorted = df[source_col].sum()
    sorted_freqs = df[source_col].values / total_sorted

    # Sample from outgrowth culture
    df['n_cells_pelleted'] = np.random.multinomial(n=int(miniprep_input_cells), pvals=sorted_freqs)

    # Determine plasmids availble for isolation
    potential_plasmids = df['n_cells_pelleted'].values * plasmid_copy_number

    # Sample from those plasmids according to kit efficiency
    df['n_plasmids_eluted'] = np.random.binomial(
        n=potential_plasmids.astype(np.int64),
        p=miniprep_efficiency
    )

    # Aliquot to PCR tube
    vol_fraction = template_vol / elution_vol
    df['Tv'] = np.random.binomial(
        n=df['n_plasmids_eluted'].values.astype(np.int64),
        p=vol_fraction
    )

    return df


def simulate_PCR(
        df,
        pcr_cycles,
        pcr_efficiency,
        pcr_noise
):
    # Determine the Gamma function parameters
    mean_gain = (1 + pcr_efficiency) ** pcr_cycles
    shape = df['Tv'].values * pcr_noise
    scale = mean_gain / pcr_noise

    # Sample from distribution for amplicon counts
    amplicons = np.zeros(len(df))
    mask = shape > 0
    if np.any(mask):
        amplicons[mask] = np.random.gamma(shape[mask], scale)
    df['n_amplicons'] = amplicons

    return df


def simulate_sequencing(
        df,
        seq_depth,
        source
):
    # Calculate amplicon frequencies
    total_molecules = df['n_amplicons'].sum()
    amp_freqs = df['n_amplicons'].values / total_molecules

    # Sample
    df[f'{source}_read_counts'] = np.random.multinomial(n=seq_depth, pvals=amp_freqs)

    return df


def simulate_experiment(
        df,
        num_cells_sorted,
        RT,
        wt_dG_fold,
        wt_dG_bind,
        ligand_conc,
        hill_coeff,
        t_exp,
        t_bind,
        facs_noise,
        miniprep_input_cells,
        plasmid_copy_number,
        miniprep_efficiency,
        elution_vol,
        template_vol,
        pcr_cycles,
        pcr_efficiency,
        pcr_noise,
        seq_depth
):
    # Outgrowth post-transformation
    df = scale_library_to_physical_count(
        df,
        num_cells_sorted
    )

    # Sequence library
    df = simulate_template_extraction(
        df,
        miniprep_input_cells,
        plasmid_copy_number,
        miniprep_efficiency,
        elution_vol,
        template_vol,
        source_col="count"
    )
    df = simulate_PCR(
        df,
        pcr_cycles,
        pcr_efficiency,
        pcr_noise
        )
    df = simulate_sequencing(
        df,
        seq_depth,
        source='lib'
        )

    # Expose to ligand
    df = simulate_biophysics(
        df,
        RT,
        wt_dG_fold,
        wt_dG_bind,
        ligand_conc,
        hill_coeff
        )
    
    # Sort
    df = simulate_facs(
        df,
        t_exp,
        t_bind,
        facs_noise
        )
    
    # Sequence sorted population
    df = simulate_template_extraction(
        df,
        miniprep_input_cells,
        plasmid_copy_number,
        miniprep_efficiency,
        elution_vol,
        template_vol,
        source_col="n_cells_sorted"
        )
    df = simulate_PCR(
        df,
        pcr_cycles,
        pcr_efficiency,
        pcr_noise
        )
    df = simulate_sequencing(
        df,
        seq_depth,
        source='sort'
        )
    
    return df