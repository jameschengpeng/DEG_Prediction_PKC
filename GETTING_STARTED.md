# DEG Prediction PKC - Getting Started Guide

This guide will help you get started with the DEG prediction analysis for PKC knockout in astrocytes.

## Prerequisites

- Python 3.8 or higher
- Basic familiarity with command line
- (Optional) Jupyter Notebook for interactive analysis

## Step-by-Step Setup Instructions

### 1. Install Python Dependencies

Open a terminal (PowerShell on Windows) and navigate to the project directory:

```powershell
cd c:\Users\james\CBIL\Astrocyte\DEG_Prediction_PKC
```

Install the required Python packages:

```powershell
pip install -r requirements.txt
```

**Note:** If you encounter any errors, you may need to install packages individually or use a virtual environment (see Optional Setup below).

### 2. Download Required Datasets

**IMPORTANT:** Data files are stored on D: drive at `D:\DEG_Prediction_PKC\data` to save SSD space.

You need to download two main datasets:

#### A. GSE43217 Proxy Dataset

1. Visit the GEO database: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43217
2. Click on "Series Matrix File(s)" and download the file
3. Save it to: `D:\DEG_Prediction_PKC\data\raw\` directory
4. The file path is already configured in `scripts/config.py` to point to D: drive

**Alternative automated download:** 
The script `step1_download_preprocess.py` can automatically download from GEO using the GEOparse library.

#### B. Astrocyte Baseline Expression Data

You need to select and download astrocyte-specific expression data. Recommended sources:

- **Zhang et al. Brain Cell Atlas**: https://www.brainrnaseq.org/
- **GTEx Portal** (brain tissue): https://gtexportal.org/
- **Allen Brain Atlas**: https://portal.brain-map.org/
- **Single-cell RNA-seq databases**: Search for astrocyte populations

After downloading:
1. Save the file to: `D:\DEG_Prediction_PKC\data\raw\` directory
2. Update the file path in `scripts/config.py` (line ~31):
   ```python
   ASTROCYTE_EXPRESSION_FILE = os.path.join(RAW_DATA_DIR, 'your_astrocyte_file.csv')
   ```

### 3. Verify Configuration

Open `scripts/config.py` and verify:

- Data directory is correctly set to `D:\DEG_Prediction_PKC\data`
- File paths are correct
- Gene panel includes your genes of interest (default includes Ca²⁺ signaling genes)
- Analysis parameters are appropriate (p-value thresholds, fold change cutoffs)

### 4. Run the Analysis

You have two options:

#### Option A: Run All Steps at Once

```powershell
python scripts/run_analysis.py
```

This will execute all 5 steps in sequence:
1. Download and preprocess GSE43217
2. Perform DEG analysis
3. Map gene panel to results
4. Generate predictions
5. Create visualizations

#### Option B: Run Steps Individually

```powershell
# Step 1: Download and preprocess
python scripts/run_analysis.py 1

# Step 2: DEG analysis
python scripts/run_analysis.py 2

# Step 3: Map gene panel
python scripts/run_analysis.py 3

# Step 4: Generate predictions
python scripts/run_analysis.py 4

# Step 5: Visualize results
python scripts/run_analysis.py 5
```

Running steps individually is useful for:
- Debugging errors
- Checking intermediate results
- Customizing parameters between steps

### 5. Check Results

After successful execution, check the following directories:

- **`data/processed/`**: Processed datasets
- **`results/tables/`**: Summary statistics
- **`results/figures/`**: Plots and visualizations
- **`predictions/`**: Final prediction tables (CSV and Excel)

The main output file is:
```
predictions/pkc_ko_predictions.xlsx
```

## Understanding the Output

### Predictions File Columns:

- **gene**: Gene symbol
- **pathway**: Functional pathway/category
- **proxy_regulation**: Direction of change in proxy data (up/down/not significant)
- **proxy_log2fc**: Log2 fold change from proxy data
- **expressed_in_astrocytes**: Whether gene is expressed in astrocytes
- **predicted_change**: Final prediction for PKC KO in astrocytes (up/down/no_change)
- **confidence**: Confidence level (high/medium/low/very_low)
- **rationale**: Biological reasoning for the prediction

### Figures:

1. **prediction_summary_by_pathway.png**: Bar chart showing predicted changes by pathway
2. **predictions_heatmap.png**: Heatmap of all predicted changes
3. **confidence_distribution.png**: Pie chart of prediction confidence levels
4. **volcano_plot.png**: Volcano plot highlighting gene panel in proxy data

## Troubleshooting Common Issues

### Issue 1: "Module not found" error

**Solution:** Install the missing package:
```powershell
pip install <package_name>
```

### Issue 2: "File not found" error

**Solution:** 
- Check that you downloaded the datasets
- Verify file paths in `scripts/config.py`
- Make sure file names match exactly

### Issue 3: "Could not automatically identify sample groups"

**Solution:** 
- Open the metadata file: `data/processed/gse43217_processed_metadata.csv`
- Manually identify which samples are control vs PKC-inhibited
- Edit `step2_deg_analysis.py` function `assign_groups()` to correctly assign samples

### Issue 4: GEOparse download fails

**Solution:**
- Download the series matrix file manually from GEO website
- Place it in `data/raw/` directory
- Update the file path in config.py

## Optional: Using a Virtual Environment

To avoid conflicts with other Python projects, use a virtual environment:

```powershell
# Create virtual environment
python -m venv venv

# Activate it
.\venv\Scripts\Activate.ps1

# Install requirements
pip install -r requirements.txt

# Run analysis
python scripts/run_analysis.py
```

## Customizing the Analysis

### Modifying the Gene Panel

Edit `scripts/config.py` and update the `GENE_PANEL` list (around line 60):

```python
GENE_PANEL = [
    'ITPR1', 'ITPR2', 'ITPR3',  # Your genes here
    # Add more genes as needed
]
```

### Adjusting Analysis Parameters

In `scripts/config.py`, modify `DEG_PARAMS`:

```python
DEG_PARAMS = {
    'p_value_threshold': 0.05,           # Raw p-value cutoff
    'adj_p_value_threshold': 0.05,       # FDR-adjusted p-value
    'log2fc_threshold': 0.5,             # Fold change cutoff (log2)
    'expression_threshold': 1.0,         # Min expression in astrocytes
}
```

### Customizing Mechanistic Logic

Edit `scripts/step4_generate_predictions.py`, function `apply_mechanistic_logic()` to modify how predictions are made based on your biological knowledge.

## Next Steps

After completing the basic analysis:

1. **Review predictions**: Check the Excel file and figures
2. **Validate findings**: Compare with literature
3. **Refine gene panel**: Add/remove genes based on results
4. **Integrate with Ca²⁺ model**: Use predictions in your calcium signaling model
5. **Plan experiments**: Design validation experiments for high-confidence predictions

## Getting Help

If you encounter issues not covered in this guide:

1. Check the log files in `results/` directory
2. Review the error messages carefully
3. Ensure all input data is in the correct format
4. Try running steps individually to isolate the problem

## Tips for Beginners in Omics Analysis

- **Start small**: Run with a subset of genes first to test the pipeline
- **Check intermediate outputs**: Review files in `data/processed/` after each step
- **Understand the data**: Look at the raw expression values and metadata
- **Read the code comments**: Each script has detailed comments explaining the logic
- **Keep notes**: Document any manual changes or decisions you make

Good luck with your analysis!
