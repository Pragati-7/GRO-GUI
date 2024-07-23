import tkinter as tk
from tkinter import filedialog
import subprocess
import os

class Equilibration(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Equilibration")
        self.create_widgets()
    def create_widgets(self):
        # Input file selection
        self.label_input_gro = tk.Label(self, text="Input EM Gro File:")
        self.label_input_gro.grid(row=0, column=0)
        self.entry_input_gro = tk.Entry(self, width=50)
        self.entry_input_gro.grid(row=0, column=1)
        self.button_browse_gro = tk.Button(self, text="Browse", command=self.select_input_gro)
        self.button_browse_gro.grid(row=0, column=2)

        # Execute equilibration
        self.button_equilibrate = tk.Button(self, text="Perform Equilibration", command=self.perform_equilibration)
        self.button_equilibrate.grid(row=1, column=1)
    def select_input_gro(self):
        file_path = filedialog.askopenfilename(filetypes=[("EM Gro Files", "*.gro")])
        self.entry_input_gro.delete(0, tk.END)
        self.entry_input_gro.insert(0, file_path)
    def perform_equilibration(self):
        input_gro = self.entry_input_gro.get()

        if input_gro:
            self.run_nvt_equilibration(input_gro)
            self.run_npt_equilibration()
        else:
            print("Please select the input file.")
    def run_nvt_equilibration(self, input_gro):
        try:
            output_dir = os.path.dirname(input_gro)
            mdp_file = os.path.join(output_dir, "nvt.mdp")
            topol_file = os.path.join(output_dir, "topol.top")
            nvt_tpr = os.path.join(output_dir, "nvt.tpr")

            # Generate nvt.tpr with grompp
            grompp_command = f"gmx grompp -f {mdp_file} -c {input_gro} -r {input_gro} -p {topol_file} -o {nvt_tpr}"
            subprocess.run(grompp_command, shell=True, check=True)
            
            # Run mdrun for NVT equilibration
            os.chdir(output_dir)
            mdrun_command = "gmx mdrun -deffnm nvt"
            subprocess.run(mdrun_command, shell=True, check=True)

            # Analyze temperature
            echo_process = subprocess.Popen(['echo', '16 0'], stdout=subprocess.PIPE)
            energy_command = f"gmx energy -f nvt.edr -o temperature.xvg"
            subprocess.run(energy_command, stdin=echo_process.stdout, shell=True, check=True)
            print("NVT equilibration completed successfully.")

            # Plot temperature
            plot_temperature_command = f"xmgrace {os.path.join(output_dir, 'temperature.xvg')}"
            subprocess.Popen(plot_temperature_command, shell=True)

        except Exception as e:
            print("Error during NVT equilibration:", e)
    def run_npt_equilibration(self):
        try:
            output_dir = os.path.dirname(self.entry_input_gro.get())
            npt_mdp_file = os.path.join(output_dir, "npt.mdp")
            nvt_gro_file = os.path.join(output_dir, "nvt.gro")
            nvt_cpt_file = os.path.join(output_dir, "nvt.cpt")
            topol_file = os.path.join(output_dir, "topol.top")
            npt_tpr = os.path.join(output_dir, "npt.tpr")
            # Check if required files exist
            if not all(map(os.path.exists, [npt_mdp_file, nvt_gro_file, nvt_cpt_file, topol_file])):
                raise FileNotFoundError("One or more required files not found.")
            # Generate npt.tpr with grompp
            grompp_command = f"gmx grompp -f {npt_mdp_file} -c {nvt_gro_file} -r {nvt_gro_file} -t {nvt_cpt_file} -p {topol_file} -o {npt_tpr}"
            subprocess.run(grompp_command, shell=True, check=True)
            # Run mdrun for NPT equilibration
            os.chdir(output_dir)
            mdrun_command = "gmx mdrun -deffnm npt"
            subprocess.run(mdrun_command, shell=True, check=True)
            # Analyze pressure
            echo_process = subprocess.Popen(['echo', '18 0'], stdout=subprocess.PIPE)
            pressure_command = f"gmx energy -f npt.edr -o pressure.xvg"
            subprocess.run(pressure_command, stdin=echo_process.stdout, shell=True, check=True)
            plot_pressure_command = f"xmgrace {os.path.join(output_dir, 'pressure.xvg')}"
            subprocess.Popen(plot_pressure_command, shell=True)
            # Analyze density
            echo_process = subprocess.Popen(['echo', '24 0'], stdout=subprocess.PIPE)
            density_command = f"gmx energy -f npt.edr -o density.xvg"
            subprocess.run(density_command, stdin=echo_process.stdout, shell=True, check=True)
            plot_density_command = f"xmgrace {os.path.join(output_dir, 'density.xvg')} -title 'title'"
            subprocess.Popen(plot_density_command, shell=True)
            print("NPT equilibration completed successfully.")
        except subprocess.CalledProcessError as e:
            print("Error during NPT equilibration. Command returned non-zero exit status:", e)
        except Exception as e:
            print("An unexpected error occurred during NPT equilibration:", e)
if __name__ == "__main__":
    app = Equilibration()
    app.mainloop()
