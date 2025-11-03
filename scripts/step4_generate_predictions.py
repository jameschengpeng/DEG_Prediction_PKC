"""
Step 4: Generate Final Predictions
Combine mechanistic reasoning with proxy data to predict PKC KO effects in astrocytes
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


def load_gene_panel_data(gene_panel_file):
    """
    Load integrated gene panel data
    
    Parameters:
    -----------
    gene_panel_file : str
        Path to gene panel file
        
    Returns:
    --------
    gene_panel_df : pd.DataFrame
        Gene panel data
    """
    logger.info("Loading gene panel data...")
    
    gene_panel_df = pd.read_csv(gene_panel_file)
    logger.info(f"Loaded {len(gene_panel_df)} genes")
    
    return gene_panel_df


def apply_mechanistic_logic(gene_panel_df):
    """
    Apply mechanistic reasoning to predict direction of change in astrocyte PKC KO
    
    This function encodes biological knowledge about how PKC regulates
    calcium signaling components
    
    Parameters:
    -----------
    gene_panel_df : pd.DataFrame
        Gene panel with DEG results and expression data
        
    Returns:
    --------
    predictions_df : pd.DataFrame
        Predictions with mechanistic rationale
    """
    logger.info("Applying mechanistic logic...")
    
    predictions = []
    
    for idx, row in gene_panel_df.iterrows():
        gene = row['gene_symbol']
        pathway = row['pathway']
        deg_regulation = row['regulation'] if 'regulation' in row else 'not_found'
        log2fc = row['log2_fold_change'] if 'log2_fold_change' in row else np.nan
        expressed = row.get('expressed_in_astrocytes', True)
        
        # Initialize prediction
        predicted_change = 'unknown'
        confidence = 'low'
        rationale = ''
        
        # PKC isoforms themselves
        if gene in ['PRKCA', 'PRKCB']:
            predicted_change = 'down'
            confidence = 'high'
            rationale = 'Direct KO target - conventional PKC isoforms removed'
        
        # Phospholipase C - upstream of PKC, may be downregulated via feedback
        elif pathway == 'Phospholipase C':
            if deg_regulation == 'downregulated':
                predicted_change = 'down'
                confidence = 'medium'
                rationale = 'Proxy shows downregulation; reduced PKC may decrease PLC via positive feedback'
            elif deg_regulation == 'upregulated':
                predicted_change = 'up'
                confidence = 'medium'
                rationale = 'Proxy shows upregulation; may indicate compensatory response'
            else:
                predicted_change = 'no_change'
                confidence = 'low'
                rationale = 'No significant change in proxy; effect uncertain in astrocytes'
        
        # IP3 receptors - PKC can phosphorylate and modulate IP3R
        elif pathway == 'IP3 Receptor':
            if deg_regulation == 'upregulated':
                predicted_change = 'up'
                confidence = 'medium'
                rationale = 'Proxy shows upregulation; loss of PKC-mediated inhibition may increase IP3R'
            elif deg_regulation == 'downregulated':
                predicted_change = 'down'
                confidence = 'medium'
                rationale = 'Proxy shows downregulation; may indicate complex regulation'
            else:
                predicted_change = 'no_change'
                confidence = 'low'
                rationale = 'No change in proxy; PKC effect on IP3R may be context-dependent'
        
        # SERCA/PMCA pumps - can be regulated by PKC
        elif pathway in ['SERCA Pump', 'PMCA Pump']:
            if deg_regulation == 'upregulated':
                predicted_change = 'up'
                confidence = 'medium'
                rationale = 'Proxy shows upregulation; loss of PKC may reduce inhibitory phosphorylation'
            elif deg_regulation == 'downregulated':
                predicted_change = 'down'
                confidence = 'medium'
                rationale = 'Proxy shows downregulation; may indicate transcriptional regulation'
            else:
                predicted_change = 'no_change'
                confidence = 'low'
                rationale = 'No change in proxy; pump regulation may be post-translational'
        
        # Store-operated calcium entry (SOCE)
        elif pathway == 'SOCE':
            if deg_regulation == 'upregulated':
                predicted_change = 'up'
                confidence = 'medium'
                rationale = 'Proxy shows upregulation; may be compensatory for altered Ca2+ handling'
            elif deg_regulation == 'downregulated':
                predicted_change = 'down'
                confidence = 'medium'
                rationale = 'Proxy shows downregulation; reduced Ca2+ signaling demands'
            else:
                predicted_change = 'no_change'
                confidence = 'low'
                rationale = 'No change in proxy; SOCE may not be directly PKC-regulated'
        
        # Other pathways
        else:
            if deg_regulation == 'upregulated':
                predicted_change = 'up'
                confidence = 'low'
                rationale = f'Proxy shows upregulation; direct PKC mechanism unclear for {pathway}'
            elif deg_regulation == 'downregulated':
                predicted_change = 'down'
                confidence = 'low'
                rationale = f'Proxy shows downregulation; direct PKC mechanism unclear for {pathway}'
            else:
                predicted_change = 'no_change'
                confidence = 'low'
                rationale = 'No change in proxy and no clear mechanistic link to PKC'
        
        # Adjust confidence based on astrocyte expression
        if not expressed:
            confidence = 'very_low'
            rationale += '; Gene not highly expressed in astrocytes'
        
        # Store prediction
        predictions.append({
            'gene': gene,
            'pathway': pathway,
            'proxy_regulation': deg_regulation,
            'proxy_log2fc': log2fc,
            'expressed_in_astrocytes': expressed,
            'predicted_change': predicted_change,
            'confidence': confidence,
            'rationale': rationale
        })
    
    predictions_df = pd.DataFrame(predictions)
    
    # Summary statistics
    logger.info("\nPrediction Summary:")
    logger.info(f"  Predicted upregulated: {(predictions_df['predicted_change'] == 'up').sum()}")
    logger.info(f"  Predicted downregulated: {(predictions_df['predicted_change'] == 'down').sum()}")
    logger.info(f"  Predicted no change: {(predictions_df['predicted_change'] == 'no_change').sum()}")
    logger.info(f"  Unknown: {(predictions_df['predicted_change'] == 'unknown').sum()}")
    
    logger.info("\nConfidence Distribution:")
    for conf in ['high', 'medium', 'low', 'very_low']:
        count = (predictions_df['confidence'] == conf).sum()
        logger.info(f"  {conf}: {count}")
    
    return predictions_df


def save_predictions(predictions_df, output_csv, output_excel):
    """
    Save final predictions to CSV and Excel
    
    Parameters:
    -----------
    predictions_df : pd.DataFrame
        Final predictions
    output_csv : str
        Output CSV file path
    output_excel : str
        Output Excel file path
    """
    logger.info(f"Saving predictions...")
    
    # Sort by pathway and confidence
    confidence_order = {'high': 4, 'medium': 3, 'low': 2, 'very_low': 1}
    predictions_df['confidence_rank'] = predictions_df['confidence'].map(confidence_order)
    predictions_df = predictions_df.sort_values(['pathway', 'confidence_rank'], ascending=[True, False])
    predictions_df = predictions_df.drop('confidence_rank', axis=1)
    
    # Save to CSV
    predictions_df.to_csv(output_csv, index=False)
    logger.info(f"Saved predictions to CSV: {output_csv}")
    
    # Save to Excel with formatting
    try:
        with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
            predictions_df.to_excel(writer, sheet_name='Predictions', index=False)
            
            # Get worksheet
            worksheet = writer.sheets['Predictions']
            
            # Adjust column widths
            for column in worksheet.columns:
                max_length = 0
                column_letter = column[0].column_letter
                for cell in column:
                    try:
                        if len(str(cell.value)) > max_length:
                            max_length = len(str(cell.value))
                    except:
                        pass
                adjusted_width = min(max_length + 2, 50)
                worksheet.column_dimensions[column_letter].width = adjusted_width
        
        logger.info(f"Saved predictions to Excel: {output_excel}")
    
    except Exception as e:
        logger.warning(f"Could not save Excel file: {str(e)}")
        logger.info("CSV file saved successfully")


def generate_summary_stats(predictions_df):
    """
    Generate summary statistics for the predictions
    
    Parameters:
    -----------
    predictions_df : pd.DataFrame
        Final predictions
        
    Returns:
    --------
    summary_df : pd.DataFrame
        Summary statistics by pathway
    """
    logger.info("Generating summary statistics...")
    
    summary_stats = []
    
    for pathway in predictions_df['pathway'].unique():
        pathway_data = predictions_df[predictions_df['pathway'] == pathway]
        
        summary_stats.append({
            'pathway': pathway,
            'total_genes': len(pathway_data),
            'predicted_up': (pathway_data['predicted_change'] == 'up').sum(),
            'predicted_down': (pathway_data['predicted_change'] == 'down').sum(),
            'no_change': (pathway_data['predicted_change'] == 'no_change').sum(),
            'high_confidence': (pathway_data['confidence'] == 'high').sum(),
            'medium_confidence': (pathway_data['confidence'] == 'medium').sum(),
            'low_confidence': (pathway_data['confidence'] == 'low').sum(),
        })
    
    summary_df = pd.DataFrame(summary_stats)
    
    # Save summary
    summary_file = os.path.join(config.TABLES_DIR, 'prediction_summary_by_pathway.csv')
    summary_df.to_csv(summary_file, index=False)
    logger.info(f"Saved summary statistics to: {summary_file}")
    
    return summary_df


def main():
    """
    Main function to generate final predictions
    """
    logger.info("Starting Step 4: Generate Final Predictions")
    
    try:
        # Load gene panel data
        gene_panel_df = load_gene_panel_data(config.OUTPUT_FILES['gene_panel_expression'])
        
        # Apply mechanistic logic
        predictions_df = apply_mechanistic_logic(gene_panel_df)
        
        # Save predictions
        save_predictions(
            predictions_df,
            config.OUTPUT_FILES['final_predictions'],
            config.OUTPUT_FILES['final_predictions_excel']
        )
        
        # Generate summary statistics
        summary_df = generate_summary_stats(predictions_df)
        
        logger.info("\n" + "="*60)
        logger.info("Step 4 completed successfully!")
        logger.info(f"Final predictions saved to:")
        logger.info(f"  CSV: {config.OUTPUT_FILES['final_predictions']}")
        logger.info(f"  Excel: {config.OUTPUT_FILES['final_predictions_excel']}")
        logger.info("="*60)
        
    except Exception as e:
        logger.error(f"Error in Step 4: {str(e)}")
        raise


if __name__ == "__main__":
    main()
