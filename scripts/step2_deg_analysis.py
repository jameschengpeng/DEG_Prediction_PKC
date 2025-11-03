"""
Step 2: DEG Analysis on GSE43217 dataset
Identify differentially expressed genes between PKC-inhibited and control samples
"""

import os
import sys
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import logging

# Add parent directory to path to import config
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import config

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def load_processed_data(data_file, metadata_file):
    """
    Load processed expression data and metadata
    
    Parameters:
    -----------
    data_file : str
        Path to processed expression data
    metadata_file : str
        Path to metadata file
        
    Returns:
    --------
    expr_df : pd.DataFrame
        Expression matrix
    metadata_df : pd.DataFrame
        Sample metadata
    """
    logger.info("Loading processed data...")
    
    expr_df = pd.read_csv(data_file, index_col=0)
    metadata_df = pd.read_csv(metadata_file)
    
    logger.info(f"Loaded expression data: {expr_df.shape}")
    logger.info(f"Loaded metadata: {metadata_df.shape}")
    
    return expr_df, metadata_df


def assign_groups(metadata_df):
    """
    Assign samples to treatment groups (control vs PKC-inhibited)
    
    NOTE: This function needs to be customized based on the actual
    sample characteristics in GSE43217. Check the metadata to identify
    which samples are controls and which are treated with PKC inhibitor.
    
    Parameters:
    -----------
    metadata_df : pd.DataFrame
        Sample metadata
        
    Returns:
    --------
    groups : dict
        Dictionary with 'control' and 'treated' sample lists
    """
    logger.info("Assigning samples to groups...")
    
    # TODO: Customize this based on actual metadata
    # This is a placeholder - you'll need to inspect the metadata
    # and modify this logic to correctly identify control vs treated samples
    
    # Example logic (modify as needed):
    control_samples = []
    treated_samples = []
    
    for idx, row in metadata_df.iterrows():
        sample_info = str(row['characteristics']).lower()
        
        # Look for keywords indicating PKC inhibition
        if 'inhibitor' in sample_info or 'pkc' in sample_info or 'treated' in sample_info:
            treated_samples.append(row['sample_id'])
        else:
            control_samples.append(row['sample_id'])
    
    # If automatic assignment fails, prompt user to manually specify
    if len(control_samples) == 0 or len(treated_samples) == 0:
        logger.warning("Could not automatically identify sample groups!")
        logger.warning("Please manually inspect metadata and update this function.")
        logger.info("Metadata sample:")
        print(metadata_df.head())
        
        # Fallback: split samples evenly (not recommended for real analysis)
        all_samples = metadata_df['sample_id'].tolist()
        mid = len(all_samples) // 2
        control_samples = all_samples[:mid]
        treated_samples = all_samples[mid:]
    
    groups = {
        'control': control_samples,
        'treated': treated_samples
    }
    
    logger.info(f"Control samples: {len(control_samples)}")
    logger.info(f"Treated samples: {len(treated_samples)}")
    
    return groups


def perform_deg_analysis(expr_df, groups, p_threshold, log2fc_threshold):
    """
    Perform differential expression analysis using t-test
    
    Parameters:
    -----------
    expr_df : pd.DataFrame
        Expression matrix
    groups : dict
        Dictionary with 'control' and 'treated' sample lists
    p_threshold : float
        P-value threshold for significance
    log2fc_threshold : float
        Log2 fold change threshold
        
    Returns:
    --------
    deg_results : pd.DataFrame
        DEG analysis results
    """
    logger.info("Performing DEG analysis...")
    
    # Separate gene annotation from expression values
    annotation_cols = [col for col in expr_df.columns if col not in groups['control'] + groups['treated']]
    
    # Get expression values for each group
    control_expr = expr_df[groups['control']]
    treated_expr = expr_df[groups['treated']]
    
    # Calculate statistics for each gene
    results = []
    
    for idx in expr_df.index:
        control_values = control_expr.loc[idx].values
        treated_values = treated_expr.loc[idx].values
        
        # Calculate means
        control_mean = np.mean(control_values)
        treated_mean = np.mean(treated_values)
        
        # Calculate log2 fold change
        log2fc = treated_mean - control_mean
        
        # Perform t-test
        try:
            t_stat, p_value = stats.ttest_ind(treated_values, control_values)
        except:
            t_stat, p_value = np.nan, 1.0
        
        # Get gene symbol if available
        gene_symbol = expr_df.loc[idx, 'Gene_Symbol'] if 'Gene_Symbol' in expr_df.columns else idx
        
        results.append({
            'probe_id': idx,
            'gene_symbol': gene_symbol,
            'control_mean': control_mean,
            'treated_mean': treated_mean,
            'log2_fold_change': log2fc,
            't_statistic': t_stat,
            'p_value': p_value
        })
    
    # Create results dataframe
    deg_results = pd.DataFrame(results)
    
    # Apply multiple testing correction (Benjamini-Hochberg FDR)
    deg_results['adj_p_value'] = multipletests(deg_results['p_value'], method='fdr_bh')[1]
    
    # Classify genes as up/down/not significant
    deg_results['regulation'] = 'not_significant'
    deg_results.loc[
        (deg_results['adj_p_value'] < p_threshold) & (deg_results['log2_fold_change'] > log2fc_threshold),
        'regulation'
    ] = 'upregulated'
    deg_results.loc[
        (deg_results['adj_p_value'] < p_threshold) & (deg_results['log2_fold_change'] < -log2fc_threshold),
        'regulation'
    ] = 'downregulated'
    
    # Sort by adjusted p-value
    deg_results = deg_results.sort_values('adj_p_value')
    
    # Print summary
    n_up = (deg_results['regulation'] == 'upregulated').sum()
    n_down = (deg_results['regulation'] == 'downregulated').sum()
    n_total = len(deg_results)
    
    logger.info(f"Total genes analyzed: {n_total}")
    logger.info(f"Upregulated genes: {n_up}")
    logger.info(f"Downregulated genes: {n_down}")
    logger.info(f"Not significant: {n_total - n_up - n_down}")
    
    return deg_results


def save_deg_results(deg_results, output_file):
    """
    Save DEG analysis results
    
    Parameters:
    -----------
    deg_results : pd.DataFrame
        DEG analysis results
    output_file : str
        Output file path
    """
    logger.info(f"Saving DEG results to {output_file}")
    
    deg_results.to_csv(output_file, index=False)
    
    # Also save significant DEGs separately
    sig_degs = deg_results[deg_results['regulation'] != 'not_significant']
    sig_output = output_file.replace('.csv', '_significant.csv')
    sig_degs.to_csv(sig_output, index=False)
    
    logger.info(f"Saved {len(sig_degs)} significant DEGs to {sig_output}")


def main():
    """
    Main function to perform DEG analysis
    """
    logger.info("Starting Step 2: DEG Analysis")
    
    try:
        # Load processed data
        data_file = config.OUTPUT_FILES['processed_gse43217']
        metadata_file = data_file.replace('.csv', '_metadata.csv')
        
        expr_df, metadata_df = load_processed_data(data_file, metadata_file)
        
        # Assign samples to groups
        groups = assign_groups(metadata_df)
        
        # Perform DEG analysis
        deg_results = perform_deg_analysis(
            expr_df,
            groups,
            config.DEG_PARAMS['adj_p_value_threshold'],
            config.DEG_PARAMS['log2fc_threshold']
        )
        
        # Save results
        save_deg_results(deg_results, config.OUTPUT_FILES['deg_results'])
        
        logger.info("Step 2 completed successfully!")
        logger.info(f"DEG results saved to: {config.OUTPUT_FILES['deg_results']}")
        
    except Exception as e:
        logger.error(f"Error in Step 2: {str(e)}")
        raise


if __name__ == "__main__":
    main()
