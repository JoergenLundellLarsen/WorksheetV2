import os
import sys
import pytest

# Add aparent directory to the system path so moduless can be imported
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.ljts.orchestrator import Orchestrator

def test_orchestrator_reads_json():
    """
    Test that the Orchestrator correctly reads values from a json configuration file.
    """
    orchestrator = Orchestrator("config/test_config.json")
    
    # Sanity checks on values to ensure they were properly read
    assert orchestrator.temperature > 0
    assert orchestrator.Lx > 0
    assert orchestrator.max_displacement > 0


def test_output_files_created():
    """
    Run the simulation and test that key output files are created.
    """
    # Run main.py with a test config
    exit_code = os.system("python main.py config/test_config.json")
    assert exit_code == 0, "main.py failed to run"

    # Check if important output files were generated
    assert os.path.exists("data/results_log.txt"), "results_log.txt not found"
    assert os.path.exists("data/config_init.xyz"), "config_init.xyz not found"

