"""
Step 4: Generate Final Predictions
Predict protein-level signaling changes in PKC KO astrocytes
Based on PKC's phosphorylation targets and signaling mechanisms
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
    Predict protein-level signaling activity changes based on PKC's mechanistic roles
    
    Key principle: PKC KO affects PROTEIN ACTIVITY (phosphorylation) not mRNA levels
    - Transcriptional changes are secondary/compensatory (long-term)
    - Immediate effects are on calcium signaling dynamics (protein function)
    
    Parameters:
    -----------
    gene_panel_df : pd.DataFrame
        Gene panel with DEG results and expression data
        
    Returns:
    --------
    predictions_df : pd.DataFrame
        Predictions of protein activity/signaling changes
    """
    logger.info("Applying mechanistic logic for protein-level signaling predictions...")
    
    predictions = []
    
    for idx, row in gene_panel_df.iterrows():
        gene = row['gene_symbol']
        pathway = row['pathway']
        deg_regulation = row['regulation'] if 'regulation' in row else 'not_found'
        log2fc = row['log2_fold_change'] if 'log2_fold_change' in row else np.nan
        expressed = row.get('expressed_in_astrocytes', True)
        
        # Initialize prediction for SIGNALING ACTIVITY (not mRNA)
        predicted_signaling_change = 'unknown'
        predicted_mrna_change = 'minimal'  # Default: no transcriptional change expected
        confidence = 'low'
        rationale = ''
        
        # ===== PKC isoforms - Direct KO targets =====
        if gene in ['PRKCA', 'PRKCB', 'PRKCG']:
            predicted_signaling_change = 'loss_of_function'
            predicted_mrna_change = 'down'  # Gene knocked out
            confidence = 'high'
            rationale = ('Direct KO target: PKC protein absent. '
                        'SIGNALING EFFECT: Complete loss of PKC kinase activity. '
                        'mRNA may show compensatory upregulation in proxy data.')
        
        elif gene in ['PRKCD', 'PRKCE', 'PRKCH', 'PRKCI', 'PRKCQ', 'PRKCZ']:
            predicted_signaling_change = 'no_change'
            predicted_mrna_change = 'minimal'
            confidence = 'medium'
            rationale = ('Non-targeted PKC isoform: Protein remains functional. '
                        'May show compensatory activity changes but not knocked out.')
        
        # ===== IP3 Receptors - PKC phosphorylation targets =====
        elif pathway == 'IP3 Receptor':
            predicted_signaling_change = 'increased_activity'
            predicted_mrna_change = 'minimal'
            confidence = 'high'
            rationale = ('PKC phosphorylates IP3R → reduces IP3R sensitivity/open probability. '
                        'PKC KO → Loss of inhibitory phosphorylation → Enhanced IP3R activity. '
                        'SIGNALING EFFECT: Increased Ca²⁺ release from ER stores. '
                        'mRNA levels likely unchanged (post-translational regulation).')
        
        # ===== Phospholipase C - Upstream of DAG-PKC pathway =====
        elif pathway == 'Phospholipase C':
            predicted_signaling_change = 'potentially_increased'
            predicted_mrna_change = 'minimal'
            confidence = 'medium'
            rationale = ('PLC generates DAG → activates PKC (positive feedback). '
                        'PKC KO → Loss of positive feedback → May reduce PLC activity indirectly. '
                        'SIGNALING EFFECT: Potentially reduced PLC-mediated IP3 production. '
                        'However, compensation via other pathways possible.')
        
        # ===== SERCA Pumps - PKC modulates pump activity =====
        elif pathway == 'SERCA Pump':
            predicted_signaling_change = 'increased_activity'
            predicted_mrna_change = 'minimal'
            confidence = 'medium'
            rationale = ('PKC can phosphorylate SERCA or associated proteins (PLB) → reduces pump efficiency. '
                        'PKC KO → Loss of inhibitory phosphorylation → Enhanced SERCA activity. '
                        'SIGNALING EFFECT: Faster ER Ca²⁺ uptake, reduced cytosolic Ca²⁺ duration. '
                        'mRNA levels likely unchanged.')
        
        # ===== PMCA Pumps - PKC phosphorylation targets =====
        elif pathway == 'PMCA Pump':
            predicted_signaling_change = 'altered_activity'
            predicted_mrna_change = 'minimal'
            confidence = 'medium'
            rationale = ('PKC phosphorylates PMCA → can enhance or inhibit pump activity (isoform-dependent). '
                        'PKC KO → Loss of PKC-mediated modulation → Altered PMCA regulation. '
                        'SIGNALING EFFECT: Changed plasma membrane Ca²⁺ extrusion kinetics. '
                        'Direction depends on specific PMCA isoform and PKC interaction.')
        
        # ===== SOCE (STIM1/ORAI) - PKC regulates SOCE activation =====
        elif pathway == 'SOCE':
            predicted_signaling_change = 'potentially_increased'
            predicted_mrna_change = 'minimal'
            confidence = 'medium'
            rationale = ('PKC can phosphorylate STIM1/ORAI → inhibits store-operated Ca²⁺ entry. '
                        'PKC KO → Loss of inhibitory phosphorylation → Enhanced SOCE. '
                        'SIGNALING EFFECT: Increased Ca²⁺ influx following ER store depletion. '
                        'May compensate for altered internal Ca²⁺ handling.')
        
        # ===== Calcium Channels - PKC modulates channel activity =====
        elif pathway == 'Calcium Channel':
            predicted_signaling_change = 'altered_gating'
            predicted_mrna_change = 'minimal'
            confidence = 'low'
            rationale = ('PKC phosphorylates various Ca²⁺ channels → modulates gating/inactivation. '
                        'PKC KO → Loss of PKC-mediated modulation → Altered channel kinetics. '
                        'SIGNALING EFFECT: Changed Ca²⁺ influx dynamics. '
                        'Direction depends on specific channel type (L-type, T-type, etc.).')
        
        # ===== Calcium Buffers - Post-translational regulation possible =====
        elif pathway == 'Calcium Buffer':
            predicted_signaling_change = 'no_direct_change'
            predicted_mrna_change = 'minimal'
            confidence = 'low'
            rationale = ('Ca²⁺ buffer proteins (calmodulin, calbindin) not direct PKC targets. '
                        'SIGNALING EFFECT: Buffer capacity unchanged at protein level. '
                        'May show secondary changes due to altered Ca²⁺ dynamics.')
        
        # ===== G-proteins - PKC can phosphorylate GPCRs =====
        elif pathway == 'G-protein':
            predicted_signaling_change = 'potentially_altered'
            predicted_mrna_change = 'minimal'
            confidence = 'low'
            rationale = ('PKC phosphorylates some GPCRs → desensitization/internalization. '
                        'PKC KO → Loss of GPCR desensitization → Potentially enhanced G-protein signaling. '
                        'SIGNALING EFFECT: Altered receptor sensitivity to agonists.')
        
        # Default case
        else:
            predicted_signaling_change = 'unknown'
            predicted_mrna_change = 'minimal'
            confidence = 'very_low'
            rationale = 'Unclear mechanistic link between PKC and this pathway component.'
        
        # Adjust confidence based on astrocyte expression
        if not expressed:
            confidence = 'very_low'
            rationale += ' NOTE: Gene not highly expressed in astrocytes - effect may be minimal.'
        
        # Store prediction
        predictions.append({
            'gene': gene,
            'pathway': pathway,
            'predicted_signaling_change': predicted_signaling_change,
            'predicted_mrna_change': predicted_mrna_change,
            'proxy_mrna_regulation': deg_regulation,
            'proxy_log2fc': log2fc,
            'expressed_in_astrocytes': expressed,
            'confidence': confidence,
            'mechanistic_rationale': rationale
        })
    
    predictions_df = pd.DataFrame(predictions)
    
    # Summary statistics
    logger.info("\n=== PROTEIN-LEVEL SIGNALING PREDICTIONS ===")
    logger.info("\nSignaling Activity Changes:")
    for change_type in predictions_df['predicted_signaling_change'].unique():
        count = (predictions_df['predicted_signaling_change'] == change_type).sum()
        logger.info(f"  {change_type}: {count}")
    
    logger.info("\nmRNA Changes (transcriptional):")
    for change_type in predictions_df['predicted_mrna_change'].unique():
        count = (predictions_df['predicted_mrna_change'] == change_type).sum()
        logger.info(f"  {change_type}: {count}")
    
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
    Generate summary statistics for the signaling predictions
    
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
        
        # Count signaling activity changes
        increased_activity = pathway_data['predicted_signaling_change'].str.contains(
            'increased|enhanced', case=False, na=False).sum()
        decreased_activity = pathway_data['predicted_signaling_change'].str.contains(
            'decreased|reduced|loss', case=False, na=False).sum()
        altered_activity = pathway_data['predicted_signaling_change'].str.contains(
            'altered|changed', case=False, na=False).sum()
        no_change = (pathway_data['predicted_signaling_change'] == 'no_change').sum()
        
        # Count mRNA changes
        mrna_down = (pathway_data['predicted_mrna_change'] == 'down').sum()
        mrna_minimal = (pathway_data['predicted_mrna_change'] == 'minimal').sum()
        
        summary_stats.append({
            'pathway': pathway,
            'total_genes': len(pathway_data),
            'signaling_increased': increased_activity,
            'signaling_decreased': decreased_activity,
            'signaling_altered': altered_activity,
            'signaling_no_change': no_change,
            'mrna_down': mrna_down,
            'mrna_minimal': mrna_minimal,
            'high_confidence': (pathway_data['confidence'] == 'high').sum(),
            'medium_confidence': (pathway_data['confidence'] == 'medium').sum(),
            'low_confidence': (pathway_data['confidence'] == 'low').sum(),
        })
    
    summary_df = pd.DataFrame(summary_stats)
    
    # Save summary
    summary_file = os.path.join(config.TABLES_DIR, 'signaling_prediction_summary_by_pathway.csv')
    summary_df.to_csv(summary_file, index=False)
    logger.info(f"Saved summary statistics to: {summary_file}")
    
    return summary_df


def main():
    """
    Main function to generate signaling-level predictions
    """
    logger.info("="*60)
    logger.info("Step 4: Protein-Level Signaling Predictions")
    logger.info("="*60)
    logger.info("NOTE: PKC KO primarily affects PROTEIN ACTIVITY (phosphorylation),")
    logger.info("      not mRNA levels. Transcriptional changes are secondary/compensatory.")
    logger.info("="*60)
    
    try:
        # Load gene panel data
        gene_panel_df = load_gene_panel_data(config.OUTPUT_FILES['gene_panel_expression'])
        
        # Apply mechanistic logic for signaling predictions
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
        logger.info("="*60)
        logger.info("INTERPRETATION GUIDE:")
        logger.info("  - 'predicted_signaling_change': Protein activity/function changes")
        logger.info("  - 'predicted_mrna_change': Expected transcriptional changes")
        logger.info("  - Proxy mRNA data shows compensatory responses, not direct PKC effects")
        logger.info("="*60)
        logger.info(f"Signaling predictions saved to:")
        logger.info(f"  CSV: {config.OUTPUT_FILES['final_predictions']}")
        logger.info(f"  Excel: {config.OUTPUT_FILES['final_predictions_excel']}")
        logger.info("="*60)
        
    except Exception as e:
        logger.error(f"Error in Step 4: {str(e)}")
        raise


if __name__ == "__main__":
    main()
