import tkinter as tk
from tkinter import filedialog
import subprocess
import os

class Analysis(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("System Analysis")
        self.create_widgets()
    def create_widgets(self):
        # Input XTC file selection
        self.label_input_xtc = tk.Label(self, text="Input XTC File:")
        self.label_input_xtc.grid(row=0, column=0)
        self.entry_input_xtc = tk.Entry(self, width=50)
        self.entry_input_xtc.grid(row=0, column=1)
        self.button_browse_xtc = tk.Button(self, text="Browse", command=self.select_input_xtc)
        self.button_browse_xtc.grid(row=0, column=2)

        # Execute trjconv, rms, and gyrate
        self.button_analyze = tk.Button(self, text="Analyze System", command=self.analyze_system)
        self.button_analyze.grid(row=1, column=1)
    def select_input_xtc(self):
        file_path = filedialog.askopenfilename(filetypes=[("XTC Files", "*.xtc")])
        self.entry_input_xtc.delete(0, tk.END)
        self.entry_input_xtc.insert(0, file_path)
    def analyze_system(self):
        input_xtc = self.entry_input_xtc.get()

        if input_xtc:
            output_dir = os.path.dirname(input_xtc)
            tpr_file = os.path.join(output_dir, "md_0_1.tpr")
            output_xtc = os.path.join(output_dir, "md_0_1_noPBC.xtc")
            output_rmsd = os.path.join(output_dir, "rmsd.xvg")
            output_rmsd_xtal = os.path.join(output_dir, "rmsd_xtal.xvg")
            output_gyrate = os.path.join(output_dir, "gyrate.xvg")

            self.execute_trjconv(tpr_file, input_xtc, output_xtc)
            self.execute_rms(tpr_file, output_xtc, output_rmsd)
            self.execute_rms_xtal(os.path.join(output_dir, "em.tpr"), output_xtc, output_rmsd_xtal)

            # Visualize RMSD files
            self.visualize_xvg(output_rmsd, output_rmsd_xtal)

            # Execute gyrate
            self.execute_gyrate(tpr_file, output_xtc, output_gyrate)

            # Visualize gyrate.xvg
            self.visualize_xvg(output_gyrate)
        else:
            print("Please select the input XTC file.")
    def execute_trjconv(self, tpr_file, input_xtc, output_xtc):
        try:
            trjconv_command = f"gmx trjconv -s {tpr_file} -f {input_xtc} -o {output_xtc} -pbc mol -center <<< $'1\n0'"
            subprocess.run(trjconv_command, shell=True, check=True)
            print("trjconv executed successfully.")
        except Exception as e:
            print("Error during trjconv execution:", e)
    def execute_rms(self, tpr_file, input_xtc, output_rmsd):
        try:
            echo_process = subprocess.Popen(['echo', '4', '4'], stdout=subprocess.PIPE)
            rms_command = f"gmx rms -s {tpr_file} -f {input_xtc} -o {output_rmsd} -tu ns"
            subprocess.run(rms_command, stdin=echo_process.stdout, shell=True, check=True)
            print("RMS analysis executed successfully.")
        except Exception as e:
            print("Error during RMS analysis:", e)
    def execute_rms_xtal(self, em_tpr_file, input_xtc, output_rmsd_xtal):
        try:
            echo_process = subprocess.Popen(['echo', '4', '4'], stdout=subprocess.PIPE)
            rms_command = f"gmx rms -s {em_tpr_file} -f {input_xtc} -o {output_rmsd_xtal} -tu ns"
            subprocess.run(rms_command, stdin=echo_process.stdout, shell=True, check=True)
            print("RMS analysis relative to crystal structure executed successfully.")
        except Exception as e:
            print("Error during RMS analysis relative to crystal structure:", e)
    def execute_gyrate(self, tpr_file, input_xtc, output_gyrate):
        try:
            echo_process = subprocess.Popen(['echo', '1'], stdout=subprocess.PIPE)
            gyrate_command = f"gmx gyrate -s {tpr_file} -f {input_xtc} -o {output_gyrate}"
            subprocess.run(gyrate_command,stdin=echo_process.stdout, shell=True, check=True)
            print("gyrate executed successfully.")
        except Exception as e:
            print("Error during gyrate execution:", e)
    def visualize_xvg(self, *xvg_files):
        try:
            subprocess.Popen(["xmgrace"] + list(xvg_files))
            print("XVG files visualized successfully.")
        except Exception as e:
            print("Error during XVG files visualization:", e)
if __name__ == "__main__":
    app = Analysis()
    app.mainloop()
