"""
Configuration file for DEG Prediction PKC project
Update file paths after downloading datasets
"""

import os

# Project root directory
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Data directories - STORED ON D: DRIVE TO SAVE SSD SPACE
DATA_DIR = r'D:\DEG_Prediction_PKC\data'
RAW_DATA_DIR = os.path.join(DATA_DIR, 'raw')
PROCESSED_DATA_DIR = os.path.join(DATA_DIR, 'processed')

# Results directories
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')
FIGURES_DIR = os.path.join(RESULTS_DIR, 'figures')
TABLES_DIR = os.path.join(RESULTS_DIR, 'tables')

# Predictions directory
PREDICTIONS_DIR = os.path.join(PROJECT_ROOT, 'predictions')

# ============================================================================
# DATA FILE PATHS - UPDATE THESE AFTER DOWNLOADING YOUR DATA
# ============================================================================

# GSE43217 proxy dataset (PKC inhibition in HT1080 cells)
# Download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43217
GSE43217_FILE = os.path.join(RAW_DATA_DIR, 'GSE43217_series_matrix.txt')
# Alternative: if you download raw CEL files, specify the directory
GSE43217_CEL_DIR = os.path.join(RAW_DATA_DIR, 'GSE43217_RAW')

# Astrocyte baseline expression data
# Update with your chosen astrocyte dataset
ASTROCYTE_EXPRESSION_FILE = os.path.join(RAW_DATA_DIR, 'astrocyte_expression.csv')

# Optional: Additional knockout/knockdown datasets
OPTIONAL_KO_FILES = {
    # 'dataset_name': os.path.join(RAW_DATA_DIR, 'filename.csv'),
}

# ============================================================================
# ANALYSIS PARAMETERS
# ============================================================================

# DEG analysis parameters
DEG_PARAMS = {
    'p_value_threshold': 0.05,
    'adj_p_value_threshold': 0.05,  # FDR/adjusted p-value
    'log2fc_threshold': 0.5,  # Log2 fold change threshold
    'expression_threshold': 1.0,  # Minimum expression level in astrocytes
}

# Gene panel of interest (Ca2+ signaling related genes)
# This is a starting list - can be expanded based on literature
GENE_PANEL = [
    # IP3 receptors
    'ITPR1', 'ITPR2', 'ITPR3',
    
    # Phospholipase C isoforms
    'PLCB1', 'PLCB2', 'PLCB3', 'PLCB4',
    'PLCG1', 'PLCG2',
    
    # Calcium pumps
    'ATP2A2',  # SERCA2
    'ATP2B1', 'ATP2B2', 'ATP2B3', 'ATP2B4',  # PMCA isoforms
    
    # Calcium channels
    'CACNA1A', 'CACNA1C', 'CACNA1D', 'CACNA1E', 'CACNA1G', 'CACNA1H',
    
    # Calcium buffers
    'CALB1', 'CALB2', 'CALM1', 'CALM2', 'CALM3',
    
    # PKC isoforms
    'PRKCA', 'PRKCB', 'PRKCD', 'PRKCE', 'PRKCG', 'PRKCH', 'PRKCI', 'PRKCQ', 'PRKCZ',
    
    # Store-operated calcium entry
    'ORAI1', 'ORAI2', 'ORAI3',
    'STIM1', 'STIM2',
    
    # Other relevant genes
    'GNA11', 'GNAQ',  # Gq proteins
]

# Output file names
OUTPUT_FILES = {
    'processed_gse43217': os.path.join(PROCESSED_DATA_DIR, 'gse43217_processed.csv'),
    'deg_results': os.path.join(PROCESSED_DATA_DIR, 'deg_analysis_results.csv'),
    'gene_panel_expression': os.path.join(PROCESSED_DATA_DIR, 'gene_panel_astrocyte_expression.csv'),
    'final_predictions': os.path.join(PREDICTIONS_DIR, 'pkc_ko_predictions.csv'),
    'final_predictions_excel': os.path.join(PREDICTIONS_DIR, 'pkc_ko_predictions.xlsx'),
}

# Visualization parameters
VIZ_PARAMS = {
    'figure_format': 'png',
    'figure_dpi': 300,
    'figure_size': (10, 8),
}

# GEO database parameters
GEO_PARAMS = {
    'gse_id': 'GSE43217',
    'cache_dir': RAW_DATA_DIR,
}
