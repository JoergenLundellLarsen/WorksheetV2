import json

class Orchestrator:
    """
    Class responsible for loading and storing simulation parameters from a JSON configuration file.
    
    """

    def __init__(self, json_path):
        """
        Initializes the Orchestrator by reading simulation parameters from a JSON file.

        Parameters:
            json_path (str): Path to the JSON configuration file containing simulation parameters.

        """
        with open(json_path, 'r') as f:
            config = json.load(f)

        # Required simulation parameters from JSON config
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
