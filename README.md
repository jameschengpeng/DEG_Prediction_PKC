# Astrocyte cPKC KO Omics Prediction  
This repository supports the analysis of predicted gene expression changes (“DEGs”) when conventional Protein Kinase C (cPKC: PRKCA/PRKCB) is knocked out in astrocytes.

## Purpose  
- We want to *predict* which genes will go up or down when cPKC is removed in astrocytes, based on mechanistic logic + available omics signatures.  
- We do *not* necessarily need new experimental KO data, but we use public omics to infer direction of change.  
- The predictions will be paired with our calcium-signaling model (astrocyte Ca²⁺ dynamics) in the companion manuscript.

## Dataset sources  
1. Proxy PKC inhibition dataset: GSE43217 (human HT1080 cells, PKC inhibitor + stimulation)  
2. Baseline astrocyte expression data: [insert link once selected]  
3. Optional: other knockdown/knockout datasets for PRKCA/PRKCB from ENCODE/KnockTF if found.

## Workflow  
1. Download and preprocess the proxy dataset (GSE43217).  
2. Identify differentially expressed genes (DEGs) in proxy (PKC-inhibited vs control).  
3. Map the directional changes (↑/↓) from proxy to a small mechanistic gene panel relevant to astrocyte Ca²⁺ signaling (e.g., IP₃R isoforms, PLCβ, SERCA2, PMCA, buffers, TF targets).  
4. Retrieve baseline expression levels of that gene panel in astrocytes (from astrocyte atlas).  
5. Combine mechanistic reasoning + proxy direction + baseline expression to assign predictions (↑, ↓, ↔) for each gene under cPKC KO in astrocytes.  
6. Output a table (Gene | Pathway | Predicted Change | Confidence | Rationale).  
7. (Optional) Plot summary heatmap or bar-plot of predicted changes.  
8. Save all scripts, results, and predictions in this repo for transparency.

## Code structure  
- `data/` : **NOTE: Data files are stored on D: drive at `D:\DEG_Prediction_PKC\data` to save SSD space**  
  - Raw and processed dataset files are stored separately on disk D  
- `scripts/` : analysis scripts (e.g., download, preprocess, DEG analysis, mapping to gene panel).  
- `results/` : output tables, plots, predictions.  
- `README.md` : this file.  
- `requirements.txt` : list of Python (or R) libraries needed (e.g., pandas, numpy, scipy, statsmodels, seaborn).  
- `predictions/` : final prediction table in CSV or Excel format.

## Prerequisites

- Python 3.8 or higher
- Basic familiarity with command line (PowerShell on Windows)
- (Optional) Jupyter Notebook for interactive analysis

## Quick Start

### 1. Clone and Install Dependencies

Clone the repository:  
```bash  
git clone https://github.com/jameschengpeng/DEG_Prediction_PKC.git  
cd DEG_Prediction_PKC
```

Install required Python packages:
```powershell
pip install -r requirements.txt
```

**Optional:** Use a virtual environment to avoid conflicts:
```powershell
python -m venv venv
.\venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### 2. Download Required Datasets

**IMPORTANT:** Data files are stored on **D: drive** at `D:\DEG_Prediction_PKC\data` to save SSD space.

#### Dataset 1: GSE43217 (Proxy PKC inhibition data)
1. Visit GEO database: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43217
2. Download the Series Matrix File
3. Save to: `D:\DEG_Prediction_PKC\data\raw\`

**Alternative:** The script can auto-download using GEOparse library.

#### Dataset 2: Astrocyte Baseline Expression Data
Choose from recommended sources:
- **Zhang et al. Brain Cell Atlas**: https://www.brainrnaseq.org/
- **GTEx Portal**: https://gtexportal.org/
- **Allen Brain Atlas**: https://portal.brain-map.org/
- **Single-cell RNA-seq databases** with astrocyte populations

Save to: `D:\DEG_Prediction_PKC\data\raw\`

### 3. Configure File Paths

⚠️ **IMPORTANT:** Edit `scripts/config.py` and update the astrocyte data file path (line ~31):
```python
ASTROCYTE_EXPRESSION_FILE = os.path.join(RAW_DATA_DIR, 'your_astrocyte_file.csv')
```

The GSE43217 path is already configured to point to D: drive.

### 4. Run the Analysis

Run all steps at once:
```powershell
python scripts/run_analysis.py
```

Or run steps individually:
```powershell
python scripts/run_analysis.py 1  # Download & preprocess GSE43217
python scripts/run_analysis.py 2  # Perform DEG analysis
python scripts/run_analysis.py 3  # Map gene panel to results
python scripts/run_analysis.py 4  # Generate predictions with mechanistic logic
python scripts/run_analysis.py 5  # Create visualizations
```

### 5. Check Results

After successful execution:
- **`D:\DEG_Prediction_PKC\data\processed\`**: Processed datasets
- **`results/tables/`**: Summary statistics
- **`results/figures/`**: Plots and visualizations (heatmaps, bar charts, volcano plots)
- **`predictions/`**: Final prediction tables (`pkc_ko_predictions.xlsx`)

## Analysis Pipeline Details

### Step 1: Download and Preprocess (`step1_download_preprocess.py`)
- Downloads GSE43217 from GEO (or loads local file)
- Extracts expression values and sample metadata
- Normalizes expression data (log2 transform + quantile normalization)
- Saves to `D:\DEG_Prediction_PKC\data\processed\`

### Step 2: DEG Analysis (`step2_deg_analysis.py`)
- Assigns samples to control vs PKC-inhibited groups
- Performs t-test for each gene
- Applies FDR correction (Benjamini-Hochberg)
- Classifies genes as upregulated/downregulated/not significant

### Step 3: Map Gene Panel (`step3_map_gene_panel.py`)
- Filters DEG results for Ca²⁺ signaling gene panel
- Integrates with astrocyte baseline expression
- Categorizes genes by functional pathway

### Step 4: Generate Predictions (`step4_generate_predictions.py`)
- Applies mechanistic reasoning based on PKC biology
- Combines proxy data with biological knowledge
- Assigns predictions (↑/↓/↔) with confidence levels
- Generates biological rationale for each prediction

### Step 5: Visualize (`step5_visualize.py`)
- Summary bar chart by pathway
- Heatmap of predicted changes
- Confidence distribution pie chart
- Volcano plot of proxy data

## Output Files Explained

### Predictions File (`predictions/pkc_ko_predictions.xlsx`)
Contains predictions for each gene with:
- **gene**: Gene symbol
- **pathway**: Functional category (IP3 Receptor, PLC, SERCA, PMCA, PKC, SOCE, etc.)
- **proxy_regulation**: Direction in proxy data (up/down/not significant)
- **proxy_log2fc**: Log2 fold change from GSE43217
- **expressed_in_astrocytes**: Expression status in astrocytes
- **predicted_change**: Final prediction (up/down/no_change/unknown)
- **confidence**: Confidence level (high/medium/low/very_low)
- **rationale**: Biological reasoning for the prediction

### Gene Panel (Default)
The analysis focuses on Ca²⁺ signaling genes:
- IP₃ receptors (ITPR1/2/3)
- Phospholipase C isoforms (PLCB1-4, PLCG1-2)
- SERCA pump (ATP2A2)
- PMCA pumps (ATP2B1-4)
- Calcium channels (CACNA1 family)
- Calcium buffers (CALB1/2, CALM1-3)
- PKC isoforms (PRKCA-Z)
- SOCE components (ORAI1-3, STIM1-2)
- G-proteins (GNA11, GNAQ)

**To customize:** Edit `GENE_PANEL` list in `scripts/config.py`

## Customization

### Adjust Analysis Parameters
Edit `scripts/config.py`:
```python
DEG_PARAMS = {
    'p_value_threshold': 0.05,           # Raw p-value cutoff
    'adj_p_value_threshold': 0.05,       # FDR threshold
    'log2fc_threshold': 0.5,             # Fold change cutoff
    'expression_threshold': 1.0,         # Min expression in astrocytes
}
```

### Modify Mechanistic Logic
Edit `apply_mechanistic_logic()` function in `step4_generate_predictions.py` to incorporate your biological knowledge.

### Add Custom Genes
Update `GENE_PANEL` list in `scripts/config.py` with additional genes of interest.

## Troubleshooting

### "Module not found" error
```powershell
pip install <package_name>
```

### "File not found" error
- Verify datasets are in `D:\DEG_Prediction_PKC\data\raw\`
- Check file paths in `scripts/config.py`
- Ensure file names match exactly

### "Could not identify sample groups"
- Check metadata: `D:\DEG_Prediction_PKC\data\processed\gse43217_processed_metadata.csv`
- Manually edit `assign_groups()` in `step2_deg_analysis.py`

### GEOparse download fails
- Download series matrix manually from GEO
- Place in `D:\DEG_Prediction_PKC\data\raw\`

## Additional Documentation

- **`scripts/README.md`**: Detailed documentation for each script
- **`D:\DEG_Prediction_PKC\data\README.md`**: Data directory structure

## Tips for Beginners

1. **Start small**: Test pipeline with default gene panel first
2. **Check intermediate outputs**: Review files after each step
3. **Understand the data**: Examine expression values and metadata
4. **Read code comments**: Scripts contain detailed explanations
5. **Keep notes**: Document any manual changes or decisions

## Citation and Contact

If you use this analysis pipeline, please cite the companion manuscript on astrocyte Ca²⁺ signaling and PKC regulation.

For questions or issues, please open an issue on GitHub
