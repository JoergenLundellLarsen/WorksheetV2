import json

class Orchestrator:
    def __init__(self, json_path):
        with open(json_path, 'r') as f:
            self.config = json.load(f)

        # Extract setup
        setup = self.config["setup"]
        self.Lx = setup["Lx"]
        self.Ly = setup["Ly"]
        self.Lz = setup["Lz"]
        self.compartments = setup["compartments"]

        # Extract temperature
        self.temperature = self.config["temperature"]

        # Extract steps
        steps = self.config["steps"]
        self.total_steps = steps["total"]
        self.reset_sampling_at = steps.get("reset_sampling_at", [])

        # Output frequencies
        self.console_output_frequency = self.config["console_output"]["frequency"]
        self.trajectory_output_frequency = self.config["trajectory_output"]["frequency"]

        # Output files
        self.initial_config_file = self.config["configuration_output"]["initial"]
        self.final_config_file = self.config["configuration_output"]["final"]
        self.trajectory_file = self.config["trajectory_output"]["file"]

        # Control parameters
        control = self.config["control_parameters"]
        self.max_displacement = control["maximum_displacement"]

        # Convert the sx, sy, sz distortions into zeta values
        distortions = control["test_area_distortion"]
        self.zeta_1 = distortions[0]["sy"]
        self.zeta_2 = distortions[1]["sy"]

    def summary(self):
        print("Simulation Box:", self.Lx, self.Ly, self.Lz)
        print("Compartments:", self.compartments)
        print("Temperature:", self.temperature)
        print("Total Steps:", self.total_steps)
        print("Reset Sampling At:", self.reset_sampling_at)
        print("Console Output Every:", self.console_output_frequency)
        print("Trajectory Output Every:", self.trajectory_output_frequency)
        print("Initial Config File:", self.initial_config_file)
        print("Final Config File:", self.final_config_file)
        print("Trajectory File:", self.trajectory_file)
        print("Max Displacement:", self.max_displacement)
        print("Zeta Values:", self.zeta_1, self.zeta_2)
