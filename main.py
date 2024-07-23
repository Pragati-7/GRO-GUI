import tkinter as tk
from tkinter import messagebox
import subprocess
import logging
from PIL import Image, ImageTk

class StepRunner(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("GUINEVERE")
        
        self.create_widgets()

    def create_widgets(self):
        steps = [
            ("Generate Topology", "1Generate_topology.py", "1topology.png"),
            ("Define Box and Solvate", "2Define_box_and_solvate.py", "2defineboxandsolvent.png"),
            ("Add Ions", "3Adding_ions.py", "3addions.png"),
            ("Energy Minimization", "4Energy_minimization.py", "4energy.png"),
            ("Equilibration", "5Equilibration.py", "5equilibration.png"),
            ("Production MD", "6Production_MD.py", "6prodmd.png"),
            ("Analysis", "7Analysis.py", "7analysis.png")
        ]
        #icon size
        icon_size = (50, 50)
        max_text_length = max(len(step_name) for step_name, _, _ in steps)
        for step_name, script_name, icon_file in steps:
            icon_img = Image.open(f"./images/{icon_file}")
            icon_img = icon_img.resize(icon_size)  
            icon_img = ImageTk.PhotoImage(icon_img)
            button = tk.Button(
                self, text=step_name,
                command=lambda s=script_name: self.run_step(s),
                compound=tk.LEFT,  
                image=icon_img,  
                activebackground="#8c52ff",
                activeforeground="black",
                anchor="center",
                bd=3,
                bg="lightgray",
                cursor="hand2",
                disabledforeground="gray",
                fg="black",
                font=("Dante", 12),
                height=104,  #ADJUST AS REQUIRED 
                highlightbackground="black",
                highlightcolor="green",
                justify="center",
                width=max_text_length * 10  # Adjust width based on text length
            )
            button.image = icon_img  
            button.pack(anchor="nw") 
        logo_img = ImageTk.PhotoImage(file="./images/logo.png")
        logo_widget = tk.Label(image=logo_img, borderwidth=0, highlightthickness=0)
        logo_widget.image = logo_img
        logo_widget.place(relx=0.6, rely=0.25, anchor="center")
        description_text = (
            "GUINEVERE is a fully integrated and efficient\n"
            "Graphical User Interface (GUI) to the updated\n"
            "molecular dynamic GROMACS version 4. \nGUINEVERE"
            "is a cross platform tkinter/python based GUI\n" 
            "application designed to break the command\n"
            "line barrier and introduces a new user-freindly\n"
            "environment to run molecular dynamics\n" 
            "simulation through GROMACS.\n"
                          )
        description_label = tk.Label(
            self,
            text=description_text,
            font=("Dante", 17),
            fg="black",
            bg="#8c52ff",
            justify="center",
        )
        description_label.place(relx=0.6, rely=0.65, anchor="center")
        footer_text=(
        "KEY WORDS: Molecular dyanamics, Simulation, GROMACS, Cross- platform, GUI, GUINEVERE"
             )
        footer_label= tk.Label(
            self, 
            text=footer_text,
            font=("Dante",15),
            fg="black",
            bg="#8c52ff",
            justify="left"
        )
        footer_label.place(relx=0.24, rely=0.9, anchor="w")
    def run_step(self, script_name):
        try:
            subprocess.Popen(["python", script_name])
        except Exception as e:
            logging.error(f"Error running {script_name}: {e}")
            messagebox.showerror("Error", f"Failed to run {script_name}")
if __name__ == "__main__":
    app = StepRunner()
    app.state('zoomed')  # Maximizing the window
    # app.attributes('-zoomed', True)
    app.iconbitmap("./images/favicon.ico")
    app.configure(bg="#8c52ff")
    app.mainloop()
