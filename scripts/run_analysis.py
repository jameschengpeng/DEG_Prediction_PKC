"""
Main runner script for DEG Prediction PKC project
Executes all analysis steps in sequence
"""

import os
import sys
import logging
from datetime import datetime

# Add scripts directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Import step modules
import step1_download_preprocess
import step2_deg_analysis
import step3_map_gene_panel
import step4_generate_predictions
import step5_visualize

# Set up logging
log_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'results')
os.makedirs(log_dir, exist_ok=True)

log_file = os.path.join(log_dir, f'analysis_log_{datetime.now().strftime("%Y%m%d_%H%M%S")}.txt')

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler(sys.stdout)
    ]
)

logger = logging.getLogger(__name__)


def run_all_steps():
    """
    Run all analysis steps in sequence
    """
    logger.info("="*80)
    logger.info("Starting DEG Prediction PKC Analysis Pipeline")
    logger.info("="*80)
    
    start_time = datetime.now()
    
    steps = [
        ("Step 1: Download and Preprocess", step1_download_preprocess.main),
        ("Step 2: DEG Analysis", step2_deg_analysis.main),
        ("Step 3: Map Gene Panel", step3_map_gene_panel.main),
        ("Step 4: Generate Predictions", step4_generate_predictions.main),
        ("Step 5: Visualize Results", step5_visualize.main),
    ]
    
    completed_steps = []
    failed_step = None
    
    for step_name, step_func in steps:
        logger.info("\n" + "="*80)
        logger.info(f"Running {step_name}")
        logger.info("="*80)
        
        try:
            step_func()
            completed_steps.append(step_name)
            logger.info(f"✓ {step_name} completed successfully")
            
        except Exception as e:
            failed_step = step_name
            logger.error(f"✗ {step_name} failed with error: {str(e)}")
            logger.error("Pipeline stopped due to error")
            break
    
    # Print summary
    end_time = datetime.now()
    duration = end_time - start_time
    
    logger.info("\n" + "="*80)
    logger.info("PIPELINE SUMMARY")
    logger.info("="*80)
    logger.info(f"Start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"End time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Duration: {duration}")
    logger.info(f"\nCompleted steps ({len(completed_steps)}/{len(steps)}):")
    
    for step in completed_steps:
        logger.info(f"  ✓ {step}")
    
    if failed_step:
        logger.error(f"\n✗ Failed at: {failed_step}")
        logger.error("Please check the error messages above and fix the issue.")
        return False
    else:
        logger.info("\n✓ All steps completed successfully!")
        logger.info(f"\nLog saved to: {log_file}")
        return True


def run_single_step(step_number):
    """
    Run a single step of the analysis
    
    Parameters:
    -----------
    step_number : int
        Step number to run (1-5)
    """
    steps = {
        1: ("Step 1: Download and Preprocess", step1_download_preprocess.main),
        2: ("Step 2: DEG Analysis", step2_deg_analysis.main),
        3: ("Step 3: Map Gene Panel", step3_map_gene_panel.main),
        4: ("Step 4: Generate Predictions", step4_generate_predictions.main),
        5: ("Step 5: Visualize Results", step5_visualize.main),
    }
    
    if step_number not in steps:
        logger.error(f"Invalid step number: {step_number}")
        logger.error("Please provide a step number between 1 and 5")
        return False
    
    step_name, step_func = steps[step_number]
    
    logger.info("="*80)
    logger.info(f"Running {step_name}")
    logger.info("="*80)
    
    try:
        step_func()
        logger.info(f"\n✓ {step_name} completed successfully")
        logger.info(f"Log saved to: {log_file}")
        return True
        
    except Exception as e:
        logger.error(f"\n✗ {step_name} failed with error: {str(e)}")
        logger.error("Please check the error messages above")
        return False


def print_usage():
    """
    Print usage instructions
    """
    print("\n" + "="*80)
    print("DEG Prediction PKC - Main Runner Script")
    print("="*80)
    print("\nUsage:")
    print("  python run_analysis.py [step_number]")
    print("\nOptions:")
    print("  No arguments      : Run all steps in sequence")
    print("  step_number (1-5) : Run only the specified step")
    print("\nSteps:")
    print("  1. Download and preprocess GSE43217 dataset")
    print("  2. Perform DEG analysis on proxy data")
    print("  3. Map gene panel and integrate astrocyte expression")
    print("  4. Generate final predictions with mechanistic reasoning")
    print("  5. Create visualizations")
    print("\nExamples:")
    print("  python run_analysis.py       # Run all steps")
    print("  python run_analysis.py 1     # Run only step 1")
    print("  python run_analysis.py 5     # Run only step 5 (visualization)")
    print("="*80 + "\n")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        # No arguments - run all steps
        success = run_all_steps()
        sys.exit(0 if success else 1)
        
    elif len(sys.argv) == 2:
        if sys.argv[1] in ['-h', '--help', 'help']:
            print_usage()
            sys.exit(0)
        
        try:
            step_num = int(sys.argv[1])
            success = run_single_step(step_num)
            sys.exit(0 if success else 1)
            
        except ValueError:
            print(f"Error: Invalid step number '{sys.argv[1]}'")
            print_usage()
            sys.exit(1)
    
    else:
        print("Error: Too many arguments")
        print_usage()
        sys.exit(1)
