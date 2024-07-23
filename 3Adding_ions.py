import tkinter as tk
from tkinter import filedialog
import subprocess
import os

class AddingIons(tk.Tk):

    def __init__(self):
        super().__init__()
        self.title("Adding Ions")
        self.create_widgets()

    def create_widgets(self):
        # Input file selection
        self.label_input_gro = tk.Label(self, text="Input Solvated Gro File:")
        self.label_input_gro.grid(row=0, column=0)
        self.entry_input_gro = tk.Entry(self, width=50)
        self.entry_input_gro.grid(row=0, column=1)
        self.button_browse_gro = tk.Button(self, text="Browse", command=self.select_input_gro)
        self.button_browse_gro.grid(row=0, column=2)

        # Execute grompp and genion
        self.button_execute = tk.Button(self, text="Execute Grompp and Genion", command=self.execute_grompp_and_genion)
        self.button_execute.grid(row=1, column=1)

    def select_input_gro(self):
        file_path = filedialog.askopenfilename(filetypes=[("Solvated Gro Files", "*_solv.gro")])
        self.entry_input_gro.delete(0, tk.END)
        self.entry_input_gro.insert(0, file_path)

    def execute_grompp_and_genion(self):
        input_gro = self.entry_input_gro.get()

        if input_gro:
            self.perform_grompp_and_genion(input_gro)
        else:
            print("Please select the input file.")

    def perform_grompp_and_genion(self, input_gro):
        try:
            output_dir = os.path.dirname(input_gro)
            output_prefix = os.path.splitext(os.path.basename(input_gro))[0]
            # Locate topol.top file
            input_top = os.path.join(output_dir, "topol.top")
            # Full path to ions.mdp
            ions_mdp = os.path.join(output_dir, "ions.mdp")
            # Generate ions.tpr with grompp
            ions_tpr = os.path.join(output_dir, output_prefix + "_ions.tpr")
            grompp_command = f"gmx grompp -f {ions_mdp} -c {input_gro} -p {input_top} -o {ions_tpr}"
            subprocess.run(grompp_command, shell=True, check=True)
            # Run genion
            output_solv_ions = os.path.join(output_dir, output_prefix + "_solv_ions.gro")
            echo_process=subprocess.Popen(['echo', '13'], stdout=subprocess.PIPE)
            genion_command = f"gmx genion -s {ions_tpr} -o {output_solv_ions} -p {input_top} -pname NA -nname CL -neutral"
            subprocess.run(genion_command, stdin=echo_process.stdout, shell=True, check=True)
            print("Grompp and Genion executed successfully.")
            # Open output GRO file in VMD
            self.open_in_vmd(output_solv_ions)

        except Exception as e:
            print("Error:", e)
    def open_in_vmd(self, file_path):
        if file_path:
            vmd_command= ["vmd", file_path]
            subprocess.run(vmd_command)
        else:
            print("File path is empty.")
if __name__ == "__main__":
    app = AddingIons()
    app.mainloop()
