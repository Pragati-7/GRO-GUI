import tkinter as tk
from tkinter import filedialog
import subprocess
import os

class ProductionMD(tk.Tk):

    def __init__(self):
        super().__init__()
        self.title("Production Molecular Dynamics")
        self.create_widgets()

    def create_widgets(self):
        # Input file selection
        self.label_input_gro = tk.Label(self, text="Input NPT Gro File:")
        self.label_input_gro.grid(row=0, column=0)
        self.entry_input_gro = tk.Entry(self, width=50)
        self.entry_input_gro.grid(row=0, column=1)
        self.button_browse_gro = tk.Button(self, text="Browse", command=self.select_input_gro)
        self.button_browse_gro.grid(row=0, column=2)

        # Execute production MD
        self.button_produce_md = tk.Button(self, text="Run Production MD", command=self.run_production_md)
        self.button_produce_md.grid(row=1, column=1)

    def select_input_gro(self):
        file_path = filedialog.askopenfilename(filetypes=[("NPT Gro Files", "*.gro")])
        self.entry_input_gro.delete(0, tk.END)
        self.entry_input_gro.insert(0, file_path)

    def run_production_md(self):
        input_gro = self.entry_input_gro.get()

        if input_gro:
            self.perform_production_md(input_gro)
        else:
            print("Please select the input file.")

    def perform_production_md(self, input_gro):
        try:
            output_dir = os.path.dirname(input_gro)
            md_tpr = os.path.join(output_dir, "md_0_1.tpr")
            npt_cpt = os.path.join(output_dir, "npt.cpt")
            topol_top = os.path.join(output_dir, "topol.top")

            # Generate md_0_1.tpr with grompp
            os.chdir(output_dir)
            grompp_command = ["gmx", "grompp", "-f", "md.mdp", "-c", input_gro, "-t", npt_cpt, "-p", topol_top, "-o", md_tpr]
            subprocess.run(grompp_command, check=True, stderr=subprocess.PIPE)
            
            # Run mdrun for production MD
            os.chdir(output_dir)
            mdrun_command = ["gmx", "mdrun", "-deffnm", "md_0_1"]
            subprocess.run(mdrun_command, check=True, stderr=subprocess.PIPE)

            print("Production MD completed successfully.")

        except subprocess.CalledProcessError as e:
            print(f"Error during production MD: {e}")
            print(f"Command returned non-zero exit status {e.returncode}.")
            print("Error output:\n", e.stderr.decode())
        except Exception as e:
            print("An unexpected error occurred during production MD:", e)
if __name__ == "__main__":
    app = ProductionMD()
    app.mainloop()

