import tkinter as tk
from tkinter import filedialog
import subprocess
import os

class TopologyGenerator(tk.Tk):

    def __init__(self):
        super().__init__()
        self.title("Generate Topology")
        self.create_widgets()

    def create_widgets(self):

        # Input file selection
        self.label_input = tk.Label(self, text="Input PDB File:")
        self.label_input.grid(row=0, column=0)
        self.entry_input = tk.Entry(self, width=50)
        self.entry_input.grid(row=0, column=1)
        self.button_browse_input = tk.Button(self, text="Browse", command=self.select_input_file)
        self.button_browse_input.grid(row=0, column=2)

        # Execute filtering, visualization, and topology generation
        self.button_execute = tk.Button(self, text="Generate Topology", command=self.execute_filter)
        self.button_execute.grid(row=1, column=1)

    def filter_pdb(self, input_file):
        try:
            output_file = os.path.splitext(input_file)[0] + "_clean.pdb"
            with open(output_file, 'w') as f_out:
                grep_command = ["grep", "-v", "HOH", input_file]
                subprocess.run(grep_command, stdout=f_out, check=True)
            print("Filtering completed successfully. Output file:", output_file)
            return output_file
        except Exception as e:
            print("Error:", e)
            return None

    def select_input_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("PDB Files", "*.pdb")])
        self.entry_input.delete(0, tk.END)
        self.entry_input.insert(0, file_path)

    def execute_filter(self):
        input_file = self.entry_input.get()
        if input_file:
            filtered_file = self.filter_pdb(input_file)
            if filtered_file:
                self.visualize_and_generate_topology(filtered_file)
        else:
            print("Please select an input file.")

    def visualize_and_generate_topology(self, file_path):
        try:
            vmd_command = ["vmd", file_path]
            subprocess.run(vmd_command)
            # After VMD closes, proceed to generate topology
            self.generate_topology(file_path)
        except Exception as e:
            print("Error launching VMD:", e)

    def generate_topology(self, input_file):
        
        try:
            output_dir= os.path.dirname(input_file)
            output_prefix = os.path.splitext(input_file)[0]
            output_gro = os.path.join(output_dir,output_prefix + "_processed.gro")
            output_top = os.path.join(output_dir,"topol.top")
            output_posre = os.path.join(output_dir,"posre.itp")
            echo_process= subprocess.Popen(['echo', '15'], stdout=subprocess.PIPE)
            pdb2gmx_command = f"gmx pdb2gmx -f {input_file} -o {output_gro} -p {output_top} -i {output_posre} -water spce"
            subprocess.run(pdb2gmx_command, stdin=echo_process.stdout, shell=True)
            print("Topology generation completed successfully.")
            vmd_command= ["vmd", output_gro]
            subprocess.run(vmd_command)

        except Exception as e:
            print("Error generating topology:", e)

if __name__ == "__main__":
    app = TopologyGenerator()
    app.mainloop()