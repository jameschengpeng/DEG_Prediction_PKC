"""
Step 1: Download and preprocess GSE43217 dataset
This script downloads the proxy PKC inhibition dataset from GEO
"""

import os
import sys
import pandas as pd
import numpy as np
import GEOparse
import logging

# Add parent directory to path to import config
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import config

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def download_geo_dataset(gse_id, cache_dir):
    """
    Download dataset from GEO database
    
    Parameters:
    -----------
    gse_id : str
        GEO series ID (e.g., 'GSE43217')
    cache_dir : str
        Directory to cache downloaded data
        
    Returns:
    --------
    gse : GEOparse.GEOTypes.GSE
        GEO series object
    """
    logger.info(f"Downloading dataset {gse_id} from GEO...")
    
    try:
        gse = GEOparse.get_GEO(geo=gse_id, destdir=cache_dir)
        logger.info(f"Successfully downloaded {gse_id}")
        return gse
    except Exception as e:
        logger.error(f"Error downloading {gse_id}: {str(e)}")
        raise


def preprocess_expression_data(gse):
    """
    Preprocess expression data from GEO dataset
    
    Parameters:
    -----------
    gse : GEOparse.GEOTypes.GSE
        GEO series object
        
    Returns:
    --------
    expr_df : pd.DataFrame
        Expression matrix with genes as rows and samples as columns
    metadata_df : pd.DataFrame
        Sample metadata
    """
    logger.info("Preprocessing expression data...")
    
    # Get the platform (microarray annotation)
    platform = list(gse.gpls.values())[0]
    
    # Get expression data from all samples
    gsm_names = list(gse.gsms.keys())
    
    # Extract expression values and create dataframe
    expr_data = []
    sample_names = []
    
    for gsm_name in gsm_names:
        gsm = gse.gsms[gsm_name]
        sample_names.append(gsm_name)
        expr_data.append(gsm.table['VALUE'].values)
    
    # Get probe IDs from the first sample
    first_sample = list(gse.gsms.values())[0]
    probe_ids = first_sample.table.index
    
    # Create expression dataframe
    expr_df = pd.DataFrame(
        np.array(expr_data).T,
        index=probe_ids,
        columns=sample_names
    )
    
    # Map probe IDs to gene symbols using platform annotation
    gene_symbols = []
    
    # Check various possible column names for gene symbols in the platform
    gene_symbol_col = None
    possible_names = ['Symbol', 'Gene Symbol', 'GENE_SYMBOL', 'Gene_Symbol', 
                      'ILMN_Gene', 'ORF', 'gene_assignment', 'Gene']
    
    for col_name in possible_names:
        if col_name in platform.table.columns:
            gene_symbol_col = col_name
            logger.info(f"Found gene symbol column: {col_name}")
            break
    
    if gene_symbol_col:
        # Create a mapping from probe ID to gene symbol
        platform.table.index = platform.table.index.astype(str)
        probe_to_gene = platform.table[gene_symbol_col].to_dict()
        
        # Map each probe ID to its gene symbol
        for probe_id in probe_ids:
            probe_id_str = str(probe_id)
            if probe_id_str in probe_to_gene:
                gene_sym = str(probe_to_gene[probe_id_str])
                # Handle multiple genes or empty values
                if gene_sym and gene_sym != 'nan' and gene_sym != '':
                    # Take first gene if multiple are listed (separated by ///)
                    gene_symbols.append(gene_sym.split('///')[0].strip())
                else:
                    gene_symbols.append(probe_id_str)
            else:
                gene_symbols.append(probe_id_str)
        
        expr_df.insert(0, 'Gene_Symbol', gene_symbols)
        logger.info(f"Mapped {sum([1 for g in gene_symbols if not g.isdigit()])} probes to gene symbols")
    else:
        logger.warning("No gene symbol column found in platform annotation")
        logger.warning(f"Available platform columns: {list(platform.table.columns)}")
        expr_df.insert(0, 'Gene_Symbol', probe_ids)
    
    # Create metadata dataframe
    metadata = []
    for gsm_name in gsm_names:
        gsm = gse.gsms[gsm_name]
        metadata.append({
            'sample_id': gsm_name,
            'title': gsm.metadata['title'][0] if 'title' in gsm.metadata else '',
            'source': gsm.metadata['source_name_ch1'][0] if 'source_name_ch1' in gsm.metadata else '',
            'characteristics': str(gsm.metadata.get('characteristics_ch1', [])),
        })
    
    metadata_df = pd.DataFrame(metadata)
    
    logger.info(f"Expression matrix shape: {expr_df.shape}")
    logger.info(f"Number of samples: {len(sample_names)}")
    
    return expr_df, metadata_df


def normalize_expression_data(expr_df):
    """
    Normalize expression data (log2 transform and quantile normalization)
    
    Parameters:
    -----------
    expr_df : pd.DataFrame
        Raw expression matrix
        
    Returns:
    --------
    norm_df : pd.DataFrame
        Normalized expression matrix
    """
    logger.info("Normalizing expression data...")
    
    # Separate gene annotation columns from expression values
    annotation_cols = [col for col in expr_df.columns if not expr_df[col].dtype in [np.float64, np.int64]]
    expr_values = expr_df.drop(columns=annotation_cols)
    
    # Convert to numeric
    expr_values = expr_values.apply(pd.to_numeric, errors='coerce')
    
    # Log2 transform (add small constant to avoid log(0))
    expr_log = np.log2(expr_values + 1)
    
    # Simple quantile normalization (rank-based)
    # For each column, replace values with sorted mean across all columns
    sorted_values = np.sort(expr_log.values, axis=0)
    ranks = expr_log.rank(method='min')
    
    # Calculate mean of sorted values across samples
    sorted_mean = np.mean(sorted_values, axis=1)
    
    # Replace values with corresponding sorted mean
    norm_values = pd.DataFrame(
        np.zeros_like(expr_log.values),
        index=expr_log.index,
        columns=expr_log.columns
    )
    
    for col in expr_log.columns:
        for idx in expr_log.index:
            rank = int(ranks.loc[idx, col]) - 1
            if rank < len(sorted_mean):
                norm_values.loc[idx, col] = sorted_mean[rank]
    
    # Combine annotation columns with normalized values
    norm_df = pd.concat([expr_df[annotation_cols], norm_values], axis=1)
    
    logger.info("Normalization complete")
    
    return norm_df


def save_processed_data(expr_df, metadata_df, output_file):
    """
    Save processed expression data and metadata
    
    Parameters:
    -----------
    expr_df : pd.DataFrame
        Processed expression matrix
    metadata_df : pd.DataFrame
        Sample metadata
    output_file : str
        Output file path
    """
    logger.info(f"Saving processed data to {output_file}")
    
    # Save expression data
    expr_df.to_csv(output_file, index=True)
    
    # Save metadata
    metadata_file = output_file.replace('.csv', '_metadata.csv')
    metadata_df.to_csv(metadata_file, index=False)
    
    logger.info("Data saved successfully")


def main():
    """
    Main function to download and preprocess GSE43217 dataset
    """
    logger.info("Starting Step 1: Download and preprocess GSE43217")
    
    # Check if data directory exists
    os.makedirs(config.PROCESSED_DATA_DIR, exist_ok=True)
    
    try:
        # Download dataset from GEO
        gse = download_geo_dataset(
            config.GEO_PARAMS['gse_id'],
            config.GEO_PARAMS['cache_dir']
        )
        
        # Preprocess expression data
        expr_df, metadata_df = preprocess_expression_data(gse)
        
        # Normalize expression data
        norm_df = normalize_expression_data(expr_df)
        
        # Save processed data
        save_processed_data(
            norm_df,
            metadata_df,
            config.OUTPUT_FILES['processed_gse43217']
        )
        
        logger.info("Step 1 completed successfully!")
        logger.info(f"Processed data saved to: {config.OUTPUT_FILES['processed_gse43217']}")
        
    except Exception as e:
        logger.error(f"Error in Step 1: {str(e)}")
        raise


if __name__ == "__main__":
    main()
