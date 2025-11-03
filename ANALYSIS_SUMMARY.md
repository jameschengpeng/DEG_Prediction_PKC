# PKC Knockout Analysis: Protein-Level Signaling Predictions

## Overview

This analysis predicts the effects of PKC (Protein Kinase C) knockout in mouse astrocytes on calcium signaling pathways. **Critically, PKC primarily functions at the post-translational level (phosphorylation), not as a transcription factor.**

## Key Biological Insight

### What PKC Does
- **Primary function**: Phosphorylates target proteins (IP3 receptors, calcium channels, pumps, etc.)
- **Timescale**: Seconds to minutes (immediate protein activity changes)
- **Mechanism**: Post-translational modification

### What PKC Does NOT Do (Primarily)
- **Not a transcription factor**: Doesn't directly regulate gene expression
- **Minimal mRNA changes expected**: Transcriptional effects are secondary/compensatory
- **Long-term effects**: Days to weeks for transcriptional compensation

## Analysis Approach

### 1. What We're Predicting
**Protein-level signaling activity changes** (NOT mRNA changes)

### 2. Proxy Data Used
- **Dataset**: GSE43217 (Human HT1080 cells with PKC inhibitor Go6976)
- **Purpose**: Shows compensatory transcriptional responses to PKC loss
- **Limitation**: Small sample size (2 control vs 2 treated), no significant DEGs after FDR correction
- **Interpretation**: mRNA data reflects secondary compensatory changes, not direct PKC effects

### 3. Prediction Categories

#### High Confidence Predictions (N=6)

**PKC Isoforms (PRKCA, PRKCB, PRKCG)**
- **Signaling**: Loss of function (protein absent)
- **mRNA**: Down (knockout targets)
- **Mechanism**: Direct knockout

**IP3 Receptors (ITPR1, ITPR2, ITPR3)**
- **Signaling**: Increased activity
- **mRNA**: Minimal change
- **Mechanism**: PKC phosphorylates IP3R → reduces sensitivity/open probability
  - PKC KO → Loss of inhibitory phosphorylation → Enhanced Ca²⁺ release from ER

#### Medium Confidence Predictions (N=37)

**SERCA Pumps (ATP2A1/2/3)**
- **Signaling**: Increased activity
- **mRNA**: Minimal change
- **Mechanism**: PKC phosphorylates SERCA/PLB → reduces pump efficiency
  - PKC KO → Enhanced SERCA activity → Faster ER Ca²⁺ uptake

**PMCA Pumps (ATP2B1/2/3/4)**
- **Signaling**: Altered activity (isoform-dependent)
- **mRNA**: Minimal change
- **Mechanism**: PKC modulates PMCA activity → Changed plasma membrane Ca²⁺ extrusion

**SOCE Components (STIM1/2, ORAI-like)**
- **Signaling**: Potentially increased
- **mRNA**: Minimal change
- **Mechanism**: PKC inhibits STIM1/ORAI → KO enhances store-operated Ca²⁺ entry

**Phospholipase C (PLCB/D/G isoforms)**
- **Signaling**: Potentially increased
- **mRNA**: Minimal change
- **Mechanism**: Loss of positive feedback from PKC → May alter PLC activity

#### Low Confidence Predictions (N=19)

**Calcium Channels, Buffers, G-proteins**
- **Signaling**: Variable (altered gating, no direct change, potentially altered)
- **mRNA**: Minimal change
- **Mechanism**: Context-dependent or indirect PKC effects

## Expected Outcomes in Real PKC KO

### Immediate Effects (Minutes-Hours)
1. **Increased IP3R activity** → Enhanced ER Ca²⁺ release
2. **Enhanced SERCA/PMCA function** → Faster Ca²⁺ clearance
3. **Increased SOCE** → Enhanced Ca²⁺ influx after store depletion
4. **Altered Ca²⁺ dynamics**: Potentially faster, larger Ca²⁺ transients with quicker recovery

### Long-term Effects (Days-Weeks)
1. **Compensatory transcriptional changes**: Cells may upregulate alternative pathways
2. **Stress response activation**: Possible induction of other kinases
3. **Homeostatic adjustments**: Changes in Ca²⁺ handling protein expression

## Output Files

### Predictions
- `predictions/pkc_ko_predictions.csv` - Full predictions with mechanistic rationale
- `predictions/pkc_ko_predictions.xlsx` - Excel version with formatting

### Summary Tables
- `results/tables/signaling_prediction_summary_by_pathway.csv` - Summary by pathway

### Figures
- `results/figures/prediction_summary_by_pathway.png` - Bar chart by pathway
- `results/figures/predictions_heatmap.png` - Heatmap of signaling changes
- `results/figures/confidence_distribution.png` - Confidence levels
- `results/figures/volcano_plot.png` - Proxy mRNA vs predicted signaling

## Interpretation Guide

### Column Descriptions

1. **gene**: Gene symbol
2. **pathway**: Functional category (IP3 Receptor, SERCA Pump, PKC Isoform, etc.)
3. **predicted_signaling_change**: Predicted change in PROTEIN ACTIVITY
   - `loss_of_function`: Protein absent/inactive
   - `increased_activity`: Enhanced protein function
   - `decreased_activity`: Reduced protein function
   - `altered_activity`: Changed regulation (direction unclear)
   - `no_change`: No direct PKC effect

4. **predicted_mrna_change**: Expected transcriptional change
   - `down`: Gene expression reduced (KO targets only)
   - `minimal`: No significant mRNA change expected

5. **proxy_mrna_regulation**: What proxy dataset shows (compensatory, not direct)
6. **confidence**: Prediction confidence (high/medium/low)
7. **mechanistic_rationale**: Biological reasoning for prediction

### Key Points
- Focus on **predicted_signaling_change** for understanding functional effects
- **proxy_mrna_regulation** shows secondary compensatory responses, not direct PKC effects
- Most genes show **minimal mRNA change** despite significant signaling changes
- This is expected: PKC acts at protein level, not transcriptional level

## Validation Approaches

To validate these predictions experimentally:

1. **Calcium imaging**: Measure Ca²⁺ transients in PKC KO astrocytes
   - Expected: Enhanced IP3-evoked Ca²⁺ release
   - Expected: Faster Ca²⁺ clearance (enhanced pumps)

2. **Phosphorylation assays**: Western blot for phospho-IP3R, phospho-PLB
   - Expected: Reduced phosphorylation of PKC targets

3. **Functional assays**: Patch-clamp, FURA-2 Ca²⁺ measurements
   - Expected: Altered Ca²⁺ channel kinetics
   - Expected: Changed SOCE dynamics

4. **Protein levels**: NOT expected to change for most targets
   - Exception: PKC itself (knocked out)

## Limitations

1. **Small sample size**: Proxy dataset (2 vs 2 samples) has limited statistical power
2. **Species difference**: Human proxy data, mouse astrocyte context
3. **Cell type**: HT1080 cells ≠ astrocytes (but core PKC mechanisms conserved)
4. **Acute vs chronic**: Proxy uses inhibitor (acute), KO is chronic (compensation possible)
5. **Isoform specificity**: Different PKC isoforms have distinct targets

## Conclusion

PKC knockout primarily affects **calcium signaling dynamics at the protein level**, not mRNA expression. The predictions focus on post-translational changes in protein activity, which are the direct and immediate consequences of PKC loss. Transcriptional changes observed in proxy data represent secondary compensatory responses, not the primary biological effect.

This analysis framework correctly separates:
1. **Direct signaling effects**: Protein phosphorylation/activity (main predictions)
2. **Compensatory transcriptional responses**: mRNA changes (secondary, long-term)

---

**Analysis completed**: November 3, 2025  
**Pipeline version**: 1.0 (Protein-level signaling focus)
