"""
Step 5: Visualize Results
Create plots and figures to visualize predictions
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging

# Add parent directory to path to import config
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import config

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = config.VIZ_PARAMS['figure_dpi']


def load_predictions(predictions_file):
    """
    Load prediction results
    
    Parameters:
    -----------
    predictions_file : str
        Path to predictions file
        
    Returns:
    --------
    predictions_df : pd.DataFrame
        Predictions data
    """
    logger.info("Loading predictions...")
    
    predictions_df = pd.read_csv(predictions_file)
    logger.info(f"Loaded {len(predictions_df)} predictions")
    
    return predictions_df


def plot_prediction_summary(predictions_df, output_file):
    """
    Create bar plot summarizing signaling activity predictions by pathway
    
    Parameters:
    -----------
    predictions_df : pd.DataFrame
        Predictions data
    output_file : str
        Output file path
    """
    logger.info("Creating signaling prediction summary plot...")
    
    # Count predictions by pathway and signaling change type
    summary_data = []
    
    for pathway in predictions_df['pathway'].unique():
        pathway_data = predictions_df[predictions_df['pathway'] == pathway]
        
        # Categorize signaling changes
        increased = pathway_data['predicted_signaling_change'].str.contains(
            'increased|enhanced', case=False, na=False).sum()
        decreased = pathway_data['predicted_signaling_change'].str.contains(
            'decreased|reduced|loss', case=False, na=False).sum()
        altered = pathway_data['predicted_signaling_change'].str.contains(
            'altered|changed', case=False, na=False).sum()
        no_change = (pathway_data['predicted_signaling_change'] == 'no_change').sum() + \
                    (pathway_data['predicted_signaling_change'] == 'no_direct_change').sum()
        
        summary_data.append({
            'Pathway': pathway,
            'Increased Activity': increased,
            'Decreased Activity': decreased,
            'Altered/Complex': altered,
            'No Change': no_change,
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.set_index('Pathway')
    
    # Create stacked bar plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    summary_df.plot(
        kind='barh',
        stacked=True,
        ax=ax,
        color=['#2ca02c', '#d62728', '#ff7f0e', '#7f7f7f'],
        alpha=0.8
    )
    
    ax.set_xlabel('Number of Genes', fontsize=12, fontweight='bold')
    ax.set_ylabel('Pathway', fontsize=12, fontweight='bold')
    ax.set_title('Predicted Protein-Level Signaling Changes by Pathway\n(PKC KO in Astrocytes)',
                 fontsize=14, fontweight='bold', pad=20)
    ax.legend(title='Signaling Activity', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=config.VIZ_PARAMS['figure_dpi'], bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved summary plot to: {output_file}")


def plot_heatmap(predictions_df, output_file):
    """
    Create heatmap of predicted signaling activity changes
    
    Parameters:
    -----------
    predictions_df : pd.DataFrame
        Predictions data
    output_file : str
        Output file path
    """
    logger.info("Creating heatmap...")
    
    # Prepare data for heatmap
    # Convert signaling predictions to numeric values
    def map_signaling_change(change):
        if pd.isna(change):
            return 0
        change_lower = str(change).lower()
        if 'increased' in change_lower or 'enhanced' in change_lower:
            return 1
        elif 'decreased' in change_lower or 'reduced' in change_lower or 'loss' in change_lower:
            return -1
        elif 'altered' in change_lower or 'changed' in change_lower:
            return 0.5  # Intermediate value for complex changes
        else:
            return 0
    
    predictions_df['signaling_numeric'] = predictions_df['predicted_signaling_change'].apply(map_signaling_change)
    
    # Create pivot table
    heatmap_data = predictions_df.pivot_table(
        values='signaling_numeric',
        index='gene',
        columns='pathway',
        aggfunc='first'
    )
    
    # If pivot doesn't work well, create alternative visualization
    if heatmap_data.empty or heatmap_data.shape[1] == 0:
        # Alternative: simple heatmap with genes as rows
        heatmap_data = predictions_df[['gene', 'signaling_numeric']].set_index('gene')
        heatmap_data.columns = ['Signaling Change']
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, max(10, len(predictions_df) * 0.3)))
    
    # Create heatmap
    sns.heatmap(
        heatmap_data,
        cmap='RdBu_r',
        center=0,
        vmin=-1,
        vmax=1,
        cbar_kws={'label': 'Signaling Activity', 'ticks': [-1, 0, 0.5, 1]},
        linewidths=0.5,
        linecolor='white',
        ax=ax,
        square=False
    )
    
    # Customize colorbar
    cbar = ax.collections[0].colorbar
    cbar.ax.set_yticklabels(['Decreased', 'No Change', 'Altered', 'Increased'])
    
    ax.set_title('Predicted Protein-Level Signaling Activity Changes\n(PKC KO in Astrocytes)',
                 fontsize=14, fontweight='bold', pad=20)
    ax.set_xlabel('Pathway', fontsize=12, fontweight='bold')
    ax.set_ylabel('Gene', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=config.VIZ_PARAMS['figure_dpi'], bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved heatmap to: {output_file}")


def plot_confidence_distribution(predictions_df, output_file):
    """
    Create pie chart showing confidence distribution
    
    Parameters:
    -----------
    predictions_df : pd.DataFrame
        Predictions data
    output_file : str
        Output file path
    """
    logger.info("Creating confidence distribution plot...")
    
    # Count by confidence level
    confidence_counts = predictions_df['confidence'].value_counts()
    
    # Create pie chart
    fig, ax = plt.subplots(figsize=(10, 8))
    
    colors = ['#2ca02c', '#ff7f0e', '#d62728', '#9467bd']
    explode = [0.05 if conf == 'high' else 0 for conf in confidence_counts.index]
    
    ax.pie(
        confidence_counts.values,
        labels=confidence_counts.index,
        autopct='%1.1f%%',
        startangle=90,
        colors=colors[:len(confidence_counts)],
        explode=explode,
        textprops={'fontsize': 12, 'fontweight': 'bold'}
    )
    
    ax.set_title('Confidence Distribution of Predictions',
                 fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=config.VIZ_PARAMS['figure_dpi'], bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved confidence plot to: {output_file}")


def plot_volcano(predictions_df, output_file):
    """
    Create volcano plot of proxy data colored by predicted signaling changes
    
    Parameters:
    -----------
    predictions_df : pd.DataFrame
        Predictions data
    output_file : str
        Output file path
    """
    logger.info("Creating volcano plot...")
    
    # Filter for genes with valid p-values and fold changes
    plot_data = predictions_df.dropna(subset=['proxy_log2fc'])
    
    if len(plot_data) == 0:
        logger.warning("No valid data for volcano plot")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Color by signaling change prediction
    def get_color(signaling_change):
        if pd.isna(signaling_change):
            return '#7f7f7f'
        change_lower = str(signaling_change).lower()
        if 'increased' in change_lower or 'enhanced' in change_lower:
            return '#2ca02c'  # Green for increased
        elif 'decreased' in change_lower or 'reduced' in change_lower or 'loss' in change_lower:
            return '#d62728'  # Red for decreased
        elif 'altered' in change_lower or 'changed' in change_lower:
            return '#ff7f0e'  # Orange for altered
        else:
            return '#7f7f7f'  # Gray for no change
    
    colors = plot_data['predicted_signaling_change'].apply(get_color)
    
    # Create scatter plot
    scatter = ax.scatter(
        plot_data['proxy_log2fc'],
        np.abs(plot_data['proxy_log2fc']),  # Use absolute fold change as proxy for significance
        c=colors,
        s=100,
        alpha=0.6,
        edgecolors='black',
        linewidths=0.5
    )
    
    # Add gene labels for notable changes or PKC genes
    for idx, row in plot_data.iterrows():
        if (abs(row['proxy_log2fc']) > 0.5) or ('PKC' in str(row['gene']).upper()):
            ax.annotate(
                row['gene'],
                (row['proxy_log2fc'], abs(row['proxy_log2fc'])),
                fontsize=8,
                alpha=0.7
            )
    
    # Add threshold lines
    ax.axvline(x=0, color='black', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.axvline(x=0.5, color='gray', linestyle=':', linewidth=0.5, alpha=0.5)
    ax.axvline(x=-0.5, color='gray', linestyle=':', linewidth=0.5, alpha=0.5)
    
    ax.set_xlabel('Log2 Fold Change (Proxy mRNA Data)', fontsize=12, fontweight='bold')
    ax.set_ylabel('|Log2 Fold Change|', fontsize=12, fontweight='bold')
    ax.set_title('Proxy Data vs Predicted Signaling Changes\n(Gene Panel - PKC Inhibition)',
                 fontsize=14, fontweight='bold', pad=20)
    
    # Create custom legend for signaling predictions
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#2ca02c', label='Signaling Increased'),
        Patch(facecolor='#d62728', label='Signaling Decreased'),
        Patch(facecolor='#ff7f0e', label='Signaling Altered'),
        Patch(facecolor='#7f7f7f', label='No Change'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', title='Predicted Signaling')
    
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=config.VIZ_PARAMS['figure_dpi'], bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved volcano plot to: {output_file}")


def main():
    """
    Main function to create visualizations
    """
    logger.info("Starting Step 5: Visualize Results")
    
    try:
        # Load predictions
        predictions_df = load_predictions(config.OUTPUT_FILES['final_predictions'])
        
        # Create output directory
        os.makedirs(config.FIGURES_DIR, exist_ok=True)
        
        # Generate plots
        plot_prediction_summary(
            predictions_df,
            os.path.join(config.FIGURES_DIR, 'prediction_summary_by_pathway.png')
        )
        
        plot_heatmap(
            predictions_df,
            os.path.join(config.FIGURES_DIR, 'predictions_heatmap.png')
        )
        
        plot_confidence_distribution(
            predictions_df,
            os.path.join(config.FIGURES_DIR, 'confidence_distribution.png')
        )
        
        plot_volcano(
            predictions_df,
            os.path.join(config.FIGURES_DIR, 'volcano_plot.png')
        )
        
        logger.info("\n" + "="*60)
        logger.info("Step 5 completed successfully!")
        logger.info(f"Figures saved to: {config.FIGURES_DIR}")
        logger.info("="*60)
        
    except Exception as e:
        logger.error(f"Error in Step 5: {str(e)}")
        raise


if __name__ == "__main__":
    main()
