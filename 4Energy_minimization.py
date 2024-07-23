import tkinter as tk
from tkinter import filedialog
import subprocess
import os
class EnergyMinimization(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Energy Minimization")
        self.create_widgets()
    def create_widgets(self):
        
        # Input file selection
        self.label_input_gro = tk.Label(self, text="Input Solvated and Ions Gro File:")
        self.label_input_gro.grid(row=0, column=0)
        self.entry_input_gro = tk.Entry(self, width=50)
        self.entry_input_gro.grid(row=0, column=1)
        self.button_browse_gro = tk.Button(self, text="Browse", command=self.select_input_gro)
        self.button_browse_gro.grid(row=0, column=2)

        # Execute grompp and mdrun
        self.button_execute = tk.Button(self, text="Execute Grompp and Mdrun", command=self.execute_grompp_and_mdrun)
        self.button_execute.grid(row=1, column=1)
    def select_input_gro(self):
        file_path = filedialog.askopenfilename(filetypes=[("Solvated and Ions Gro Files", "*_solv_ions.gro")])
        self.entry_input_gro.delete(0, tk.END)
        self.entry_input_gro.insert(0, file_path)
    def execute_grompp_and_mdrun(self):
        input_gro = self.entry_input_gro.get()

        if input_gro:
            self.perform_grompp_and_mdrun(input_gro)
        else:
            print("Please select the input file.")

    def perform_grompp_and_mdrun(self, input_gro):
        try:
            output_dir = os.path.dirname(input_gro)
            output_prefix = os.path.splitext(os.path.basename(input_gro))[0]
            # Locate minim.mdp and topol.top in the same directory as input_gro
            mdp_file = os.path.join(output_dir, "minim.mdp")
            topol_file = os.path.join(output_dir, "topol.top")
            # Generate em.tpr with grompp
            em_tpr = os.path.join(output_dir, "em.tpr")
            grompp_command = f"gmx grompp -f {mdp_file} -c {input_gro} -p {topol_file} -o {em_tpr}"
            subprocess.run(grompp_command, shell=True, check=True)
            # Run mdrun
            os.chdir(output_dir)
            mdrun_command = f"gmx mdrun -v -deffnm em"
            subprocess.run(mdrun_command, shell=True, check=True)
            # Move output files to the same directory as em.tpr
            output_files = ["em.log", "em.edr", "em.trr", "em.gro"]
            for file in output_files:
                source = os.path.join(output_dir, file)
                destination = os.path.join(os.path.dirname(em_tpr), file)
                os.replace(source, destination)
            # Analyze em.edr file with gmx energy
            echo_process= subprocess.Popen(['echo','10 0'], stdout=subprocess.PIPE)
            energy_command = f"gmx energy -f em.edr -o {os.path.join(output_dir, 'potential.xvg')}"
            subprocess.run(energy_command, stdin=echo_process.stdout,shell=True, check=True)
            print("Grompp and Mdrun executed successfully.")
            plot_command = f"xmgrace {os.path.join(os.path.dirname(em_tpr), 'potential.xvg')}"
            subprocess.Popen(plot_command, shell=True)
        except Exception as e:
            print("Error:", e)
if __name__ == "__main__":
    app = EnergyMinimization()
    app.mainloop()
