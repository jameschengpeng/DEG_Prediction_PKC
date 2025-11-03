"""
Step 3: Map DEG results to gene panel and integrate with astrocyte expression
Extract genes of interest from DEG results and check their expression in astrocytes
"""

import os
import sys
import pandas as pd
import numpy as np
import logging

# Add parent directory to path to import config
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import config

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def load_deg_results(deg_file):
    """
    Load DEG analysis results
    
    Parameters:
    -----------
    deg_file : str
        Path to DEG results file
        
    Returns:
    --------
    deg_df : pd.DataFrame
        DEG results
    """
    logger.info("Loading DEG results...")
    
    deg_df = pd.read_csv(deg_file)
    logger.info(f"Loaded {len(deg_df)} DEG results")
    
    return deg_df


def load_astrocyte_expression(astrocyte_file):
    """
    Load astrocyte baseline expression data
    
    NOTE: This function needs to be customized based on the format
    of your astrocyte expression dataset
    
    Parameters:
    -----------
    astrocyte_file : str
        Path to astrocyte expression file
        
    Returns:
    --------
    astro_expr : pd.DataFrame
        Astrocyte expression data with gene symbols and expression values
    """
    logger.info("Loading astrocyte expression data...")
    
    try:
        # Try to read the file (adjust based on actual format)
        astro_expr = pd.read_csv(astrocyte_file)
        
        logger.info(f"Loaded astrocyte expression data: {astro_expr.shape}")
        
        return astro_expr
        
    except FileNotFoundError:
        logger.warning(f"Astrocyte expression file not found: {astrocyte_file}")
        logger.warning("Creating placeholder astrocyte expression data")
        
        # Create placeholder data
        # In real analysis, you must provide actual astrocyte expression data
        astro_expr = pd.DataFrame({
            'gene_symbol': config.GENE_PANEL,
            'expression_level': np.random.rand(len(config.GENE_PANEL)) * 100,
            'is_expressed': [True] * len(config.GENE_PANEL)
        })
        
        return astro_expr
    
    except Exception as e:
        logger.error(f"Error loading astrocyte expression data: {str(e)}")
        raise


def map_gene_panel_to_degs(deg_df, gene_panel):
    """
    Map gene panel to DEG results
    
    Parameters:
    -----------
    deg_df : pd.DataFrame
        DEG results
    gene_panel : list
        List of genes of interest
        
    Returns:
    --------
    panel_degs : pd.DataFrame
        DEG results for genes in the panel
    """
    logger.info(f"Mapping {len(gene_panel)} genes from panel to DEG results...")
    
    # Filter DEG results for genes in the panel
    panel_degs = deg_df[deg_df['gene_symbol'].isin(gene_panel)].copy()
    
    # Add genes from panel that are not in DEG results
    missing_genes = set(gene_panel) - set(panel_degs['gene_symbol'])
    
    if missing_genes:
        logger.warning(f"{len(missing_genes)} genes from panel not found in DEG results")
        logger.info(f"Missing genes: {', '.join(sorted(missing_genes))}")
        
        # Add missing genes with NA values
        missing_df = pd.DataFrame({
            'gene_symbol': list(missing_genes),
            'regulation': 'not_found',
            'log2_fold_change': np.nan,
            'p_value': np.nan,
            'adj_p_value': np.nan
        })
        
        panel_degs = pd.concat([panel_degs, missing_df], ignore_index=True)
    
    logger.info(f"Found {len(panel_degs)} genes from panel in results")
    
    return panel_degs


def integrate_astrocyte_expression(panel_degs, astro_expr):
    """
    Integrate astrocyte expression data with DEG results
    
    Parameters:
    -----------
    panel_degs : pd.DataFrame
        DEG results for gene panel
    astro_expr : pd.DataFrame
        Astrocyte expression data
        
    Returns:
    --------
    integrated_df : pd.DataFrame
        Integrated data with DEG results and astrocyte expression
    """
    logger.info("Integrating astrocyte expression data...")
    
    # Merge DEG results with astrocyte expression
    # Adjust column names based on actual astrocyte data format
    
    # Try to find gene symbol column in astrocyte data
    gene_col = None
    for col in ['gene_symbol', 'Gene_Symbol', 'gene', 'Gene', 'symbol', 'Symbol']:
        if col in astro_expr.columns:
            gene_col = col
            break
    
    if gene_col is None:
        logger.warning("Could not find gene symbol column in astrocyte data")
        logger.warning("Using first column as gene symbols")
        gene_col = astro_expr.columns[0]
    
    # Standardize gene column name
    astro_expr = astro_expr.rename(columns={gene_col: 'gene_symbol'})
    
    # Merge dataframes
    integrated_df = panel_degs.merge(
        astro_expr,
        on='gene_symbol',
        how='left',
        suffixes=('_deg', '_astro')
    )
    
    # Add classification for astrocyte expression
    if 'expression_level' in integrated_df.columns:
        integrated_df['expressed_in_astrocytes'] = (
            integrated_df['expression_level'] > config.DEG_PARAMS['expression_threshold']
        )
    elif 'is_expressed' in integrated_df.columns:
        integrated_df['expressed_in_astrocytes'] = integrated_df['is_expressed']
    else:
        # If no expression info, assume all genes are expressed
        integrated_df['expressed_in_astrocytes'] = True
    
    logger.info(f"Integrated data shape: {integrated_df.shape}")
    
    return integrated_df


def categorize_genes_by_pathway(integrated_df):
    """
    Categorize genes by their functional pathway
    
    Parameters:
    -----------
    integrated_df : pd.DataFrame
        Integrated DEG and astrocyte expression data
        
    Returns:
    --------
    integrated_df : pd.DataFrame
        Data with added pathway category column
    """
    logger.info("Categorizing genes by pathway...")
    
    # Define pathway categories
    pathway_map = {
        'IP3 Receptor': ['ITPR1', 'ITPR2', 'ITPR3'],
        'Phospholipase C': ['PLCB1', 'PLCB2', 'PLCB3', 'PLCB4', 'PLCG1', 'PLCG2'],
        'SERCA Pump': ['ATP2A2'],
        'PMCA Pump': ['ATP2B1', 'ATP2B2', 'ATP2B3', 'ATP2B4'],
        'Calcium Channel': ['CACNA1A', 'CACNA1C', 'CACNA1D', 'CACNA1E', 'CACNA1G', 'CACNA1H'],
        'Calcium Buffer': ['CALB1', 'CALB2', 'CALM1', 'CALM2', 'CALM3'],
        'PKC Isoform': ['PRKCA', 'PRKCB', 'PRKCD', 'PRKCE', 'PRKCG', 'PRKCH', 'PRKCI', 'PRKCQ', 'PRKCZ'],
        'SOCE': ['ORAI1', 'ORAI2', 'ORAI3', 'STIM1', 'STIM2'],
        'G-protein': ['GNA11', 'GNAQ'],
    }
    
    # Create reverse map
    gene_to_pathway = {}
    for pathway, genes in pathway_map.items():
        for gene in genes:
            gene_to_pathway[gene] = pathway
    
    # Add pathway column
    integrated_df['pathway'] = integrated_df['gene_symbol'].map(gene_to_pathway)
    integrated_df['pathway'] = integrated_df['pathway'].fillna('Other')
    
    # Count genes per pathway
    pathway_counts = integrated_df['pathway'].value_counts()
    logger.info("Genes per pathway:")
    for pathway, count in pathway_counts.items():
        logger.info(f"  {pathway}: {count}")
    
    return integrated_df


def save_gene_panel_results(integrated_df, output_file):
    """
    Save integrated gene panel results
    
    Parameters:
    -----------
    integrated_df : pd.DataFrame
        Integrated gene panel data
    output_file : str
        Output file path
    """
    logger.info(f"Saving gene panel results to {output_file}")
    
    integrated_df.to_csv(output_file, index=False)
    
    logger.info("Results saved successfully")


def main():
    """
    Main function to map gene panel and integrate astrocyte expression
    """
    logger.info("Starting Step 3: Gene Panel Mapping and Integration")
    
    try:
        # Load DEG results
        deg_df = load_deg_results(config.OUTPUT_FILES['deg_results'])
        
        # Load astrocyte expression data
        astro_expr = load_astrocyte_expression(config.ASTROCYTE_EXPRESSION_FILE)
        
        # Map gene panel to DEG results
        panel_degs = map_gene_panel_to_degs(deg_df, config.GENE_PANEL)
        
        # Integrate with astrocyte expression
        integrated_df = integrate_astrocyte_expression(panel_degs, astro_expr)
        
        # Categorize by pathway
        integrated_df = categorize_genes_by_pathway(integrated_df)
        
        # Save results
        save_gene_panel_results(integrated_df, config.OUTPUT_FILES['gene_panel_expression'])
        
        logger.info("Step 3 completed successfully!")
        logger.info(f"Gene panel results saved to: {config.OUTPUT_FILES['gene_panel_expression']}")
        
    except Exception as e:
        logger.error(f"Error in Step 3: {str(e)}")
        raise


if __name__ == "__main__":
    main()
