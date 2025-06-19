import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from src.ljts.orchestrator import Orchestrator

def test_orchestrator_reads_json():
    orchestrator = Orchestrator("config/test_config.json")
    assert orchestrator.temperature > 0
    assert orchestrator.Lx > 0
    assert orchestrator.max_displacement > 0


def test_output_files_created():
    os.system("python main.py config/test_config.json")
    assert os.path.exists("data/results_log.txt")
    assert os.path.exists("data/config_init.xyz")
