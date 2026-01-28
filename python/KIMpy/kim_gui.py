import subprocess
import threading
import tkinter as tk
from queue import Queue
from tkinter import filedialog, scrolledtext, ttk

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk


class KimGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Kim's GUI")

        # Set the initial size of the window
        self.root.geometry("800x600")  # Width x Height

        # Set the theme based on system appearance (light or dark mode)
        self.set_theme()

        # Create a Notebook (tabbed layout)
        self.notebook = ttk.Notebook(self.root)

        # Create tabs
        self.configuration_tab = ttk.Frame(self.notebook, style="My.TFrame")
        self.run_tab = ttk.Frame(self.notebook, style="My.TFrame")
        self.output_tab = ttk.Frame(self.notebook, style="My.TFrame")

        # Add tabs to the Notebook
        self.notebook.add(self.configuration_tab, text="Configuration")
        self.notebook.add(self.run_tab, text="Run")
        self.notebook.add(self.output_tab, text="Output")

        # Populate tabs with content
        self.create_configuration_tab()
        self.create_run_tab()
        self.create_output_tab()

        # Pack the Notebook into the main window
        self.notebook.pack(expand=True, fill="both")

        # Variables to store the process and queue for output
        self.process = None
        self.queue = Queue()

        self.update_output_flag = True  # Flag to control the update_output function

    def set_theme(self):
        # The theme setting remains the same as in the previous examples
        pass

    def create_configuration_tab(self):
        # Code for Configuration tab
        frame = ttk.Frame(self.configuration_tab, style="My.TFrame")
        frame.pack(padx=10, pady=10, fill="both", expand=True)

        # label = tk.Label(frame, text="Configuration settings", foreground='#FFFFFF', background='#303030', font=('Helvetica', 14))
        # label.grid(row=0, column=0, columnspan=3, pady=(0, 10))

        default_button = tk.Button(
            frame, text="Default", command=self.load_default_config, bg="#4CAF50", fg="#000000"
        )
        default_button.grid(row=1, column=0, pady=5)

        load_config_button = tk.Button(
            frame, text="Load Config", command=self.load_config, bg="#4CAF50", fg="#000000"
        )
        load_config_button.grid(row=1, column=1, pady=5)

        save_config_button = tk.Button(
            frame, text="Save Config", command=self.save_config, bg="#4CAF50", fg="#000000"
        )
        save_config_button.grid(row=1, column=2, pady=5)

        # Textboxes for displaying and manipulating the loaded name list
        self.name_list_text = scrolledtext.ScrolledText(
            frame, wrap=tk.WORD, width=70, height=30, bg="#000000", fg="#FFFFFF"
        )
        self.name_list_text.grid(row=2, column=0, columnspan=3, padx=10, pady=(0, 10))

        # Textbox for specifying the save location
        save_location_label = tk.Label(
            frame, text="Save Location:", foreground="#FFFFFF", background="#303030"
        )
        save_location_label.grid(row=3, column=0, pady=(10, 0), sticky="e")

        self.save_location_text = tk.Entry(frame, bg="#303030", fg="#000000")
        self.save_location_text.grid(row=3, column=1, pady=(10, 0), padx=5, sticky="w")

        save_location_button = tk.Button(
            frame, text="Browse", command=self.browse_save_location, bg="#4CAF50", fg="#000000"
        )
        save_location_button.grid(row=3, column=2, pady=(10, 0), padx=5, sticky="w")

    def browse_save_location(self):
        # Open a file dialog to select a save location
        save_location = filedialog.asksaveasfilename(
            defaultextension=".txt", filetypes=[("Text files", "*.nml"), ("All files", "*.*")]
        )

        if save_location:
            # Set the selected save location in the entry widget
            self.save_location_text.delete(0, tk.END)
            self.save_location_text.insert(tk.END, save_location)

    def load_default_config(self):
        # Load a default name list
        self.name_list = ["Name1", "Name2", "Name3"]
        self.display_name_list()

    def load_config(self):
        # Open a file dialog to select a name list file
        file_path = filedialog.askopenfilename(
            filetypes=[("Text files", "*.nml"), ("All files", "*.*")]
        )

        if file_path:
            # Read the name list from the file
            with open(file_path, "r") as file:
                self.name_list = [line.strip() for line in file]

            self.display_name_list()

    def save_config(self):
        # Open a file dialog to select a save location
        save_path = filedialog.asksaveasfilename(
            defaultextension=".txt", filetypes=[("Text files", "*.nml"), ("All files", "*.*")]
        )

        if save_path:
            # Write the name list to the specified location
            with open(save_path, "w") as file:
                for name in self.name_list:
                    file.write(name + "\n")

    def display_name_list(self):
        # Display the name list in the textbox
        self.name_list_text.delete(1.0, tk.END)  # Clear previous content
        for name in self.name_list:
            self.name_list_text.insert(tk.END, name + "\n")

    def create_run_tab(self):
        # Code for Run tab
        # label = tk.Label(self.run_tab, text="Run settings", foreground='#FFFFFF', background='#303030')  # Text and background colors
        # label.pack(padx=10, pady=10)

        # Add a button to run an executable
        run_button = tk.Button(
            self.run_tab,
            text="Run Executable",
            command=self.start_executable,
            bg="#4CAF50",
            fg="#FFFFFF",
        )  # Button colors
        run_button.pack(pady=10)

        # Text box to display command line output
        self.output_text = scrolledtext.ScrolledText(
            self.run_tab, wrap=tk.WORD, width=60, height=10, bg="#303030", fg="#FFFFFF"
        )  # Text box colors
        self.output_text.pack(padx=10, pady=10)

    def create_output_tab(self):
        # Code for Output tab
        # label = tk.Label(self.output_tab, text="Output display", foreground='#FFFFFF', background='#303030')  # Text and background colors
        # label.pack(padx=10, pady=10)
        frame = ttk.Frame(self.output_tab, style="My.TFrame")
        frame.pack(padx=10, pady=10, fill="both", expand=True)

        # label = tk.Label(frame, text="Output display", foreground='#FFFFFF', background='#303030', font=('Helvetica', 14))
        # label.grid(row=0, column=0, pady=(0, 10))

        # Button to plot data using Matplotlib
        plot_button = tk.Button(
            frame, text="Plot Data", command=self.plot_data, bg="#4CAF50", fg="#000000"
        )
        plot_button.grid(row=1, column=0, pady=10)

        # Text box to display command line output
        self.output_text = scrolledtext.ScrolledText(
            frame, wrap=tk.WORD, width=60, height=10, bg="#303030", fg="#000000"
        )
        self.output_text.grid(row=2, column=0, padx=10, pady=10)

    def plot_data(self):
        # Example data, replace with your own data
        x_data = [1, 2, 3, 4, 5]
        y_data = [10, 12, 8, 15, 7]

        # Create a figure and axis
        fig, ax = plt.subplots()

        # Plot the data
        ax.plot(x_data, y_data, label="Data Plot")

        # Add labels and title
        ax.set_xlabel("X-axis Label")
        ax.set_ylabel("Y-axis Label")
        ax.set_title("Matplotlib Plot")

        # Add a legend
        ax.legend()

        # Create a Tkinter window to embed the Matplotlib plot
        plot_window = tk.Toplevel(self.root)
        plot_window.title("Matplotlib Plot")

        # Embed the Matplotlib plot in the Tkinter window
        canvas = FigureCanvasTkAgg(fig, master=plot_window)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Navigation toolbar (optional)
        toolbar = NavigationToolbar2Tk(canvas, plot_window)
        toolbar.update()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    def start_executable(self):
        # Start the executable in a separate thread
        threading.Thread(target=self.run_executable).start()

    def run_executable(self):
        # Replace 'your_executable_command' with the actual command to run your executable
        command = "python test.py"

        try:
            # Start the process
            self.process = subprocess.Popen(
                command,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                universal_newlines=True,
            )
            # p = subprocess.Popen(command, stdout=subprocess.PIPE, bufsize=1)
            # for line in iter(p.stdout.readline, b''):
            #    print(line)
            # p.stdout.close()
            # p.wait()

            # Continuously read the output and put it in the queue
            # for line in iter(self.process.stdout.readline, ''):
            while True:
                line = self.process.stdout.readline()
                if not line:
                    break
                self.queue.put(line)
                self.output_text.insert(tk.END, line)
                self.output_text.yview(tk.END)  # Auto-scroll to the bottom
            self.process.stdout.close()
            # Wait for the process to complete
            self.process.wait()

        except Exception as e:
            # Display an error message if the execution fails
            self.output_text.insert(tk.END, f"Error: {str(e)}")

    def update_output(self):
        # Check if the Tkinter window is still alive before updating the output
        if self.root.winfo_exists() and self.update_output_flag:
            # Check if there is new output in the queue
            while not self.queue.empty():
                line = self.queue.get()
                self.output_text.insert(tk.END, line)
                self.output_text.yview(tk.END)  # Auto-scroll to the bottom

            # Schedule the update after 100ms (adjust as needed)
            self.root.after(100, self.update_output)

    def run(self):
        # Start the main loop
        self.root.after(100, self.update_output)
        self.root.mainloop()


# Create an instance of the KimGUI class and run the GUI
root = tk.Tk()
kim_gui = KimGUI(root)
kim_gui.run()
