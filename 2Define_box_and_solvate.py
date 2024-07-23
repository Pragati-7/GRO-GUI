import tkinter as tk
from tkinter import filedialog
import subprocess
import os

class DefineBoxAndSolvate(tk.Tk):

    def __init__(self):
        super().__init__()
        self.title("Define Box and Solvate")
        self.create_widgets()

    def create_widgets(self):
        # Input file selection
        self.label_input = tk.Label(self, text="Input Processed Gro File:")
        self.label_input.grid(row=0, column=0)
        self.entry_input = tk.Entry(self, width=50)
        self.entry_input.grid(row=0, column=1)
        self.button_browse_input = tk.Button(self, text="Browse", command=self.select_input_file)
        self.button_browse_input.grid(row=0, column=2)

        # Execute box definition and solvation
        self.button_execute = tk.Button(self, text="Execute Box Definition and Solvation", command=self.execute_box_definition_and_solvation)
        self.button_execute.grid(row=1, column=1)

    def select_input_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("Processed Gro Files", "*_processed.gro")])
        self.entry_input.delete(0, tk.END)
        self.entry_input.insert(0, file_path)

    def execute_box_definition_and_solvation(self):
        input_file = self.entry_input.get()
        if input_file:
            self.define_box_and_solvate(input_file)
        else:
            print("Please select an input file.")

    def define_box_and_solvate(self, input_file):
        try:
            output_dir = os.path.dirname(input_file)
            output_prefix = os.path.splitext(input_file)[0]
            output_box = os.path.join(output_dir, output_prefix + "_newbox.gro")
            output_solv = os.path.join(output_dir, output_prefix + "_solv.gro")
            top_file = os.path.join(output_dir, "topol.top")
            # Define box
            box_command = f"gmx editconf -f {input_file} -o {output_box} -c -d 1.0 -bt cubic"
            subprocess.run(box_command, shell=True, check=True)
            # Visualize box
            vmd_command = ["vmd", output_box]
            subprocess.run(vmd_command)
            # Solvate system
            solvate_command = f"gmx solvate -cp {output_box} -cs spc216.gro -o {output_solv} -p {top_file}"
            subprocess.run(solvate_command, shell=True, check=True)
            # Visualize solvated system
            vmd_command = ["vmd", output_solv]
            subprocess.run(vmd_command)
            print("Box definition and solvation completed successfully.")
        except Exception as e:
            print("Error:", e)

if __name__ == "__main__":
    app = DefineBoxAndSolvate()
    app.mainloop()
