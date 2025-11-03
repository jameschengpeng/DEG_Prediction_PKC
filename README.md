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
- `data/` : raw and processed dataset files.  
- `scripts/` : analysis scripts (e.g., download, preprocess, DEG analysis, mapping to gene panel).  
- `results/` : output tables, plots, predictions.  
- `README.md` : this file.  
- `requirements.txt` : list of Python (or R) libraries needed (e.g., pandas, numpy, scipy, statsmodels, seaborn).  
- `predictions/` : final prediction table in CSV or Excel format.

## Usage instructions  
1. Clone the repository:  
   ```bash  
   git clone https://github.com/jameschengpeng/DEG_Prediction_PKC.git  
   cd DEG_Prediction_PKC
