# orchestrator.py
import json

class Orchestrator:
    def __init__(self, json_path):
        with open(json_path, 'r') as f:
            config = json.load(f)

        self.Lx = config["Lx"]
        self.Ly = config["Ly"]
        self.Lz = config["Lz"]
        self.density_liquid = config["density_liquid"]
        self.density_vapor = config["density_vapor"]
        self.temperature = config["temperature"]
        self.total_sweeps = config["total_sweeps"]
        self.equilibration_sweeps = config["equilibration_sweeps"]
        self.log_interval = config["log_interval"]
        self.trajectory_interval = config["trajectory_interval"]
        self.max_displacement = config["max_displacement"]
        self.zeta = config["zeta"]
