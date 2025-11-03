# Scripts Directory

This directory contains all analysis scripts for the DEG Prediction PKC project.

## Main Scripts

### `run_analysis.py`
Main runner script that executes the entire analysis pipeline or individual steps.

**Usage:**
```bash
# Run all steps
python run_analysis.py

# Run individual step
python run_analysis.py [1-5]

# Get help
python run_analysis.py --help
```

### `config.py`
Configuration file containing:
- File paths for input/output data
- Analysis parameters (p-value thresholds, fold change cutoffs)
- Gene panel definition
- Visualization settings

**⚠️ IMPORTANT: Update file paths in this file after downloading your datasets!**

## Analysis Pipeline Steps

### Step 1: `step1_download_preprocess.py`
Downloads and preprocesses the GSE43217 proxy dataset from GEO.

**What it does:**
- Downloads GSE43217 data using GEOparse
- Extracts expression values and sample metadata
- Normalizes expression data (log2 transform + quantile normalization)
- Saves processed data to `data/processed/`

**Output:**
- `data/processed/gse43217_processed.csv`
- `data/processed/gse43217_processed_metadata.csv`

**Run individually:**
```bash
python scripts/step1_download_preprocess.py
```

### Step 2: `step2_deg_analysis.py`
Performs differential expression analysis on the proxy dataset.

**What it does:**
- Loads processed expression data
- Assigns samples to control vs PKC-inhibited groups
- Performs t-test for each gene
- Applies FDR correction (Benjamini-Hochberg)
- Classifies genes as up/down/not significant

**Output:**
- `data/processed/deg_analysis_results.csv` (all genes)
- `data/processed/deg_analysis_results_significant.csv` (significant only)

**⚠️ Note:** You may need to manually edit the `assign_groups()` function to correctly identify control vs treated samples based on the actual metadata.

**Run individually:**
```bash
python scripts/step2_deg_analysis.py
```

### Step 3: `step3_map_gene_panel.py`
Maps the gene panel to DEG results and integrates astrocyte expression data.

**What it does:**
- Loads DEG results from Step 2
- Loads astrocyte baseline expression data
- Filters DEG results for genes in the predefined gene panel
- Integrates proxy DEG data with astrocyte expression levels
- Categorizes genes by functional pathway

**Output:**
- `data/processed/gene_panel_astrocyte_expression.csv`

**Run individually:**
```bash
python scripts/step3_map_gene_panel.py
```

### Step 4: `step4_generate_predictions.py`
Generates final predictions using mechanistic reasoning.

**What it does:**
- Loads integrated gene panel data
- Applies biological/mechanistic logic to predict direction of change
- Combines proxy data direction with known PKC biology
- Assigns confidence levels based on:
  - Evidence strength from proxy data
  - Known mechanistic relationships
  - Expression in astrocytes
- Generates rationale for each prediction

**Output:**
- `predictions/pkc_ko_predictions.csv`
- `predictions/pkc_ko_predictions.xlsx`
- `results/tables/prediction_summary_by_pathway.csv`

**Customization:** Edit the `apply_mechanistic_logic()` function to modify prediction rules based on your biological knowledge.

**Run individually:**
```bash
python scripts/step4_generate_predictions.py
```

### Step 5: `step5_visualize.py`
Creates visualizations of the prediction results.

**What it does:**
- Loads final predictions
- Creates multiple plots:
  - Summary bar chart by pathway
  - Heatmap of predicted changes
  - Confidence distribution pie chart
  - Volcano plot of proxy data

**Output:**
- `results/figures/prediction_summary_by_pathway.png`
- `results/figures/predictions_heatmap.png`
- `results/figures/confidence_distribution.png`
- `results/figures/volcano_plot.png`

**Run individually:**
```bash
python scripts/step5_visualize.py
```

## Workflow Diagram

```
Step 1: Download & Preprocess
    ↓
    GSE43217 raw data → Normalized expression matrix
    ↓
Step 2: DEG Analysis
    ↓
    Control vs PKC-inhibited → DEG list with statistics
    ↓
Step 3: Map Gene Panel
    ↓
    Gene panel + Astrocyte expression → Integrated data
    ↓
Step 4: Generate Predictions
    ↓
    Mechanistic reasoning → Final predictions with confidence
    ↓
Step 5: Visualize
    ↓
    Figures and plots
```

## Customization Guide

### Adding New Genes to the Panel

Edit `config.py`:
```python
GENE_PANEL = [
    'ITPR1', 'ITPR2', 'ITPR3',
    'YOUR_GENE_1', 'YOUR_GENE_2',  # Add here
]
```

### Modifying DEG Thresholds

Edit `config.py`:
```python
DEG_PARAMS = {
    'p_value_threshold': 0.05,        # Change as needed
    'adj_p_value_threshold': 0.05,    # FDR threshold
    'log2fc_threshold': 0.5,          # Fold change cutoff
}
```

### Customizing Mechanistic Predictions

Edit `step4_generate_predictions.py`, function `apply_mechanistic_logic()`:
```python
def apply_mechanistic_logic(gene_panel_df):
    # Add your custom logic here
    # Example:
    if gene in ['YOUR_GENE']:
        if deg_regulation == 'upregulated':
            predicted_change = 'up'
            confidence = 'high'
            rationale = 'Your reasoning here'
```

### Adding New Visualizations

Edit `step5_visualize.py` and add new plotting functions. Follow the pattern of existing functions.

## Troubleshooting

### Script Fails at Step 1
- Check internet connection (for GEO download)
- Verify GEOparse is installed: `pip install GEOparse`
- Try downloading data manually and updating file paths

### Script Fails at Step 2
- Check that Step 1 completed successfully
- Verify processed data files exist in `data/processed/`
- May need to manually edit `assign_groups()` function

### Script Fails at Step 3
- Ensure astrocyte expression data is downloaded
- Check file path in `config.py`
- Verify data format matches expected structure

### Script Fails at Step 4 or 5
- Check that previous steps completed successfully
- Verify intermediate files exist

## Best Practices

1. **Always run steps in order** (unless you know what you're doing)
2. **Check intermediate outputs** after each step
3. **Keep backups** of your modified scripts
4. **Document changes** you make to the code
5. **Version control** your modifications with git

## Dependencies

See `requirements.txt` for full list. Key libraries:
- `pandas`: Data manipulation
- `numpy`: Numerical computations
- `scipy`: Statistical tests
- `GEOparse`: Download GEO datasets
- `matplotlib`, `seaborn`: Visualization
- `statsmodels`: Advanced statistics

## Performance Notes

- Step 1 (download): May take 5-15 minutes depending on internet speed
- Step 2 (DEG): Usually < 1 minute
- Step 3 (mapping): Usually < 1 minute
- Step 4 (predictions): Usually < 1 minute
- Step 5 (visualization): Usually < 2 minutes

Total runtime: ~10-20 minutes for complete analysis

## Output File Sizes (Approximate)

- Processed data: ~5-50 MB (depends on platform)
- DEG results: ~1-10 MB
- Gene panel: < 1 MB
- Predictions: < 1 MB
- Figures: ~1-5 MB total
