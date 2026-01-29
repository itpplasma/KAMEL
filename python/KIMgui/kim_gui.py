#!/usr/bin/env python3
"""
KIM GUI - A graphical interface for KIM (KiLCA Integral Model)
Allows configuration of namelist variables, running KIM code, and plotting results.
"""

import os
import subprocess
import sys
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import f90nml
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

# Add the KIMpy directory to the path
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "KIMpy"))
from kimpy import KIMpy


class KIMConfigGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("KIM Configuration GUI")

        # Initialize variables
        self.runpath = tk.StringVar(value="./kim_runs/")
        self.config_file = None
        self.kim_runner = None
        self.running = False
        self.use_gui_config = tk.BooleanVar(value=True)

        # Create main interface
        self.create_widgets()

        # Load default configuration
        self.load_default_config()

        # Bring window to front and set size to fit content
        self.setup_window()

    def create_widgets(self):
        # Create notebook for tabs
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill="both", expand=True, padx=10, pady=10)

        # Run tab (first tab)
        self.run_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.run_frame, text="Run & Results")

        # Configuration tab (second tab)
        self.config_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.config_frame, text="Configuration")

        # Create configuration widgets
        self.create_config_widgets()

        # Create run widgets
        self.create_run_widgets()

    def create_config_widgets(self):
        # Main frame with scrollbar
        canvas = tk.Canvas(self.config_frame)
        scrollbar = ttk.Scrollbar(self.config_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)

        scrollable_frame.bind(
            "<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )

        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        # Pack scrollbar and canvas
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        # Path selection
        path_frame = ttk.LabelFrame(scrollable_frame, text="Run Path", padding="10")
        path_frame.pack(fill="x", padx=5, pady=5)

        ttk.Label(path_frame, text="Run Directory:").pack(anchor="w")
        path_entry_frame = ttk.Frame(path_frame)
        path_entry_frame.pack(fill="x", pady=5)

        ttk.Entry(path_entry_frame, textvariable=self.runpath, width=50).pack(
            side="left", fill="x", expand=True
        )
        ttk.Button(path_entry_frame, text="Browse", command=self.browse_runpath).pack(
            side="right", padx=(5, 0)
        )

        # Create two-column layout for configuration sections
        columns_frame = ttk.Frame(scrollable_frame)
        columns_frame.pack(fill="x", padx=5, pady=5)

        # Left column
        left_column = ttk.Frame(columns_frame)
        left_column.pack(side="left", fill="both", expand=True, padx=(0, 5))

        # Right column
        right_column = ttk.Frame(columns_frame)
        right_column.pack(side="right", fill="both", expand=True, padx=(5, 0))

        # Configuration sections
        self.config_vars = {}

        # KIM_CONFIG section (left column)
        config_frame = ttk.LabelFrame(left_column, text="KIM Configuration", padding="10")
        config_frame.pack(fill="x", pady=5)

        self.config_vars["kim_config"] = {}
        config_fields = [
            ("profile_location", "Profile Location:", str, "./profiles/"),
            ("output_path", "Output Path:", str, "./out/"),
            ("hdf5_input", "HDF5 Input:", bool, False),
            ("hdf5_output", "HDF5 Output:", bool, False),
            ("fdebug", "Debug Level:", int, 1),
            ("fstatus", "Status Level:", int, 1),
            ("number_of_ion_species", "Number of Ion Species:", int, 1),
            ("read_species_from_namelist", "Read Species from Namelist:", bool, False),
            ("type_of_run", "Type of Run:", str, "electrostatic"),
            ("collision_model", "Collision Model:", str, "Krook"),
            ("artificial_debye_case", "Artificial Debye Case:", bool, False),
        ]

        self.create_config_fields(config_frame, config_fields, "kim_config")

        # KIM_SETUP section (left column)
        setup_frame = ttk.LabelFrame(left_column, text="KIM Setup", padding="10")
        setup_frame.pack(fill="x", pady=5)

        self.config_vars["kim_setup"] = {}
        setup_fields = [
            ("btor", "Toroidal Magnetic Field:", float, -17977.413),
            ("r0", "Major Radius (cm):", float, 165.0),
            ("m_mode", "Poloidal Mode Number:", int, 7),
            ("n_mode", "Toroidal Mode Number:", int, 2),
            ("omega", "Frequency:", float, 0.0),
            ("spline_base", "Spline Base:", int, 1),
            ("type_br_field", "Type of Br Field:", int, 12),
            ("collisions_off", "Collisions Off:", bool, True),
            ("eps_reg", "Regularization Parameter:", float, 0.01),
            ("set_profiles_constant", "Set Profiles Constant:", int, 0),
        ]

        self.create_config_fields(setup_frame, setup_fields, "kim_setup")

        # KIM_GRID section (right column)
        grid_frame = ttk.LabelFrame(right_column, text="KIM Grid", padding="10")
        grid_frame.pack(fill="x", pady=5)

        self.config_vars["kim_grid"] = {}
        grid_fields = [
            ("r_plas", "Plasma Radius (cm):", float, 67.0),
            ("r_min", "Minimum Radius (cm):", float, 3.0),
            ("width_res", "Resonance Width:", float, 0.2),
            ("ampl_res", "Resonance Amplitude:", float, 15.0),
            ("hrmax_scaling", "Hr Max Scaling:", float, 1.0),
            ("k_space_dim", "K Space Dimension:", int, 100),
            ("l_space_dim", "L Space Dimension:", int, 1000),
            ("reduce_r", "Reduce R:", bool, False),
            ("reduced_rg_dim", "Reduced RG Dimension:", int, 1000),
            ("grid_spacing", "Grid Spacing:", int, 3),
            ("num_gengrid_points", "Number of Generated Grid Points:", int, 64),
            ("kr_grid_width_res", "Kr Grid Width Resonance:", float, 5.0),
            ("kr_grid_ampl_res", "Kr Grid Amplitude Resonance:", float, 1.0),
            ("delta_l_max", "Delta L Max:", int, 5),
            ("gauss_int_nodes_Nx", "Gauss Integration Nodes Nx:", int, 7),
            ("gauss_int_nodes_Nxp", "Gauss Integration Nodes Nxp:", int, 7),
            ("gauss_int_nodes_Ntheta", "Gauss Integration Nodes Ntheta:", int, 7),
        ]

        self.create_config_fields(grid_frame, grid_fields, "kim_grid")

        # Buttons
        button_frame = ttk.Frame(scrollable_frame)
        button_frame.pack(fill="x", padx=5, pady=10)

        ttk.Button(button_frame, text="Load Config", command=self.load_config).pack(
            side="left", padx=5
        )
        ttk.Button(button_frame, text="Save Config", command=self.save_config).pack(
            side="left", padx=5
        )
        ttk.Button(button_frame, text="Reset to Default", command=self.load_default_config).pack(
            side="left", padx=5
        )

    def create_config_fields(self, parent, fields, section_name):
        for field_name, label, field_type, default_value in fields:
            frame = ttk.Frame(parent)
            frame.pack(fill="x", pady=2)

            ttk.Label(frame, text=label, width=20).pack(side="left")

            if field_type == bool:
                var = tk.BooleanVar(value=default_value)
                ttk.Checkbutton(frame, variable=var).pack(side="left")
            elif field_type in [int, float]:
                var = tk.StringVar(value=str(default_value))
                ttk.Entry(frame, textvariable=var, width=15).pack(side="left")
            elif field_name == "collision_model":
                var = tk.StringVar(value=default_value)
                combo = ttk.Combobox(
                    frame, textvariable=var, values=["Krook", "FokkerPlanck"], width=12
                )
                combo.pack(side="left")
            elif field_name == "type_of_run":
                var = tk.StringVar(value=default_value)
                combo = ttk.Combobox(
                    frame, textvariable=var, values=["electrostatic", "electromagnetic"], width=12
                )
                combo.pack(side="left")
            else:
                var = tk.StringVar(value=str(default_value))
                ttk.Entry(frame, textvariable=var, width=15).pack(side="left")

            self.config_vars[section_name][field_name] = var

    def create_run_widgets(self):
        # Run controls
        control_frame = ttk.LabelFrame(self.run_frame, text="Run Controls", padding="10")
        control_frame.pack(fill="x", padx=5, pady=5)

        # Profile generation
        profile_frame = ttk.Frame(control_frame)
        profile_frame.pack(fill="x", pady=5)

        ttk.Label(profile_frame, text="Profile Generation:").pack(side="left")
        self.gen_profiles_btn = ttk.Button(
            profile_frame, text="Generate Constant Profiles", command=self.generate_profiles
        )
        self.gen_profiles_btn.pack(side="left", padx=10)

        # Run button
        run_frame = ttk.Frame(control_frame)
        run_frame.pack(fill="x", pady=5)

        self.run_btn = ttk.Button(run_frame, text="Run KIM", command=self.run_kim)
        self.run_btn.pack(side="left")

        self.stop_btn = ttk.Button(run_frame, text="Stop", command=self.stop_kim, state="disabled")
        self.stop_btn.pack(side="left", padx=10)

        # Checkbox for using GUI configuration
        self.use_config_cb = ttk.Checkbutton(
            run_frame, text="Use GUI Configuration", variable=self.use_gui_config
        )
        self.use_config_cb.pack(side="left", padx=10)

        # Progress bar
        self.progress = ttk.Progressbar(run_frame, mode="indeterminate")
        self.progress.pack(side="left", fill="x", expand=True, padx=10)

        # Output text
        output_frame = ttk.LabelFrame(self.run_frame, text="Output", padding="10")
        output_frame.pack(fill="both", expand=True, padx=5, pady=5)

        self.output_text = tk.Text(output_frame, height=10, wrap="word")
        output_scrollbar = ttk.Scrollbar(
            output_frame, orient="vertical", command=self.output_text.yview
        )
        self.output_text.configure(yscrollcommand=output_scrollbar.set)

        self.output_text.pack(side="left", fill="both", expand=True)
        output_scrollbar.pack(side="right", fill="y")

        # Results frame
        results_frame = ttk.LabelFrame(self.run_frame, text="Results", padding="10")
        results_frame.pack(fill="x", padx=5, pady=5)

        ttk.Button(results_frame, text="Plot Results", command=self.plot_results).pack(side="left")
        ttk.Button(results_frame, text="Open Output Directory", command=self.open_output_dir).pack(
            side="left", padx=10
        )

    def browse_runpath(self):
        directory = filedialog.askdirectory()
        if directory:
            self.runpath.set(directory)

    def load_default_config(self):
        # Load default values (already set in create_config_fields)
        pass

    def load_config(self):
        filename = filedialog.askopenfilename(
            title="Select KIM Configuration File",
            filetypes=[("Namelist files", "*.nml"), ("All files", "*.*")],
        )
        if filename:
            try:
                nml = f90nml.read(filename)
                self.load_config_from_nml(nml)
                messagebox.showinfo("Success", "Configuration loaded successfully!")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load configuration: {str(e)}")

    def load_config_from_nml(self, nml):
        for section_name, section_vars in self.config_vars.items():
            if section_name in nml:
                for var_name, var_obj in section_vars.items():
                    if var_name in nml[section_name]:
                        value = nml[section_name][var_name]
                        if isinstance(var_obj, tk.BooleanVar):
                            var_obj.set(bool(value))
                        else:
                            var_obj.set(str(value))

    def save_config(self):
        filename = filedialog.asksaveasfilename(
            title="Save KIM Configuration",
            defaultextension=".nml",
            filetypes=[("Namelist files", "*.nml"), ("All files", "*.*")],
        )
        if filename:
            try:
                nml_dict = {}
                for section_name, section_vars in self.config_vars.items():
                    nml_dict[section_name] = {}
                    for var_name, var_obj in section_vars.items():
                        value = var_obj.get()
                        if isinstance(var_obj, tk.BooleanVar):
                            nml_dict[section_name][var_name] = value
                        else:
                            try:
                                if "." in str(value):
                                    nml_dict[section_name][var_name] = float(value)
                                else:
                                    nml_dict[section_name][var_name] = int(value)
                            except ValueError:
                                nml_dict[section_name][var_name] = str(value)

                # Add kim_species sections based on number_of_ion_species
                num_species = int(self.config_vars["kim_config"]["number_of_ion_species"].get())
                kim_species_sections = []

                # Create empty kim_species sections (first two are always empty)
                kim_species_sections.append({})
                kim_species_sections.append({})

                # Add populated kim_species section for the species (default deuterium)
                if num_species >= 1:
                    kim_species_sections.append({"zi": 1, "ai": 2})

                # Add any additional species sections if needed
                for _ in range(3, num_species + 2):
                    kim_species_sections.append({})

                # Write the namelist using the lower-level approach to handle multiple kim_species sections
                with open(filename, "w") as f:
                    # Write kim_config section
                    f.write("&kim_config\n")
                    for key, value in nml_dict["kim_config"].items():
                        if isinstance(value, bool):
                            f.write(f"    {key} = .{str(value).lower()}.\n")
                        elif isinstance(value, str):
                            f.write(f"    {key} = '{value}'\n")
                        else:
                            f.write(f"    {key} = {value}\n")
                    f.write("/\n\n")

                    # Write kim_setup section
                    f.write("&kim_setup\n")
                    for key, value in nml_dict["kim_setup"].items():
                        if isinstance(value, bool):
                            f.write(f"    {key} = .{str(value).lower()}.\n")
                        elif isinstance(value, str):
                            f.write(f"    {key} = '{value}'\n")
                        else:
                            f.write(f"    {key} = {value}\n")
                    f.write("/\n\n")

                    # Write kim_grid section
                    f.write("&kim_grid\n")
                    for key, value in nml_dict["kim_grid"].items():
                        if isinstance(value, bool):
                            f.write(f"    {key} = .{str(value).lower()}.\n")
                        elif isinstance(value, str):
                            f.write(f"    {key} = '{value}'\n")
                        else:
                            f.write(f"    {key} = {value}\n")
                    f.write("/\n\n")

                    # Write kim_species sections
                    for species_data in kim_species_sections:
                        f.write("&kim_species\n")
                        for key, value in species_data.items():
                            f.write(f"    {key} = {value}\n")
                        f.write("/\n\n")

                messagebox.showinfo("Success", "Configuration saved successfully!")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save configuration: {str(e)}")

    def generate_profiles(self):
        try:
            runpath = self.runpath.get()
            if not os.path.exists(runpath):
                os.makedirs(runpath)
            if not os.path.exists(os.path.join(runpath, "profiles")):
                os.makedirs(os.path.join(runpath, "profiles"))

            # Create a temporary KIMpy instance to generate profiles
            temp_kim = KIMpy(runpath)
            temp_kim.generate_constant_profiles()

            self.output_text.insert(tk.END, "Constant profiles generated successfully!\n")
            self.output_text.see(tk.END)
            messagebox.showinfo("Success", "Constant profiles generated!")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to generate profiles: {str(e)}")

    def run_kim(self):
        if self.running:
            messagebox.showwarning("Warning", "KIM is already running!")
            return

        try:
            runpath = self.runpath.get()
            if not os.path.exists(runpath):
                os.makedirs(runpath)

            # Check if KIM.x exists in runpath, if not create softlink
            kim_exe_path = os.path.join(runpath, "KIM.x")
            if not os.path.exists(kim_exe_path):
                # Get the KIM executable path from environment
                code_path = os.environ.get("CODE")
                if code_path:
                    source_kim_path = os.path.join(code_path, "KAMEL", "KIM", "build", "KIM.x")
                    if os.path.exists(source_kim_path):
                        try:
                            os.symlink(source_kim_path, kim_exe_path)
                            self.output_text.insert(
                                tk.END, f"Created softlink to KIM.x from {source_kim_path}\n"
                            )
                            self.output_text.see(tk.END)
                        except OSError as e:
                            raise Exception(f"Failed to create softlink to KIM.x: {str(e)}")
                    else:
                        raise Exception(
                            f"KIM.x not found at {source_kim_path}. Please build KIM first with 'make KIM'."
                        )
                else:
                    raise Exception(
                        "CODE environment variable not set. Please set CODE to point to your KAMEL installation."
                    )

            # Save current configuration to namelist file (if checkbox is checked)
            if self.use_gui_config.get():
                config_path = os.path.join(runpath, "KIM_config.nml")
                nml_dict = {}
                for section_name, section_vars in self.config_vars.items():
                    nml_dict[section_name] = {}
                    for var_name, var_obj in section_vars.items():
                        value = var_obj.get()
                        if isinstance(var_obj, tk.BooleanVar):
                            nml_dict[section_name][var_name] = value
                        else:
                            try:
                                if "." in str(value):
                                    nml_dict[section_name][var_name] = float(value)
                                else:
                                    nml_dict[section_name][var_name] = int(value)
                            except ValueError:
                                nml_dict[section_name][var_name] = str(value)

                # Add kim_species sections based on number_of_ion_species
                num_species = int(self.config_vars["kim_config"]["number_of_ion_species"].get())
                kim_species_sections = []

                # Create empty kim_species sections (first two are always empty)
                kim_species_sections.append({})
                kim_species_sections.append({})

                # Add populated kim_species section for the species (default deuterium)
                if num_species >= 1:
                    kim_species_sections.append({"zi": 1, "ai": 2})

                # Add any additional species sections if needed
                for _ in range(3, num_species + 2):
                    kim_species_sections.append({})

                # Write the namelist using the lower-level approach to handle multiple kim_species sections
                try:
                    with open(config_path, "w") as f:
                        # Write kim_config section
                        f.write("&kim_config\n")
                        for key, value in nml_dict["kim_config"].items():
                            if isinstance(value, bool):
                                f.write(f"    {key} = .{str(value).lower()}.\n")
                            elif isinstance(value, str):
                                f.write(f"    {key} = '{value}'\n")
                            else:
                                f.write(f"    {key} = {value}\n")
                        f.write("/\n\n")

                        # Write kim_setup section
                        f.write("&kim_setup\n")
                        for key, value in nml_dict["kim_setup"].items():
                            if isinstance(value, bool):
                                f.write(f"    {key} = .{str(value).lower()}.\n")
                            elif isinstance(value, str):
                                f.write(f"    {key} = '{value}'\n")
                            else:
                                f.write(f"    {key} = {value}\n")
                        f.write("/\n\n")

                        # Write kim_grid section
                        f.write("&kim_grid\n")
                        for key, value in nml_dict["kim_grid"].items():
                            if isinstance(value, bool):
                                f.write(f"    {key} = .{str(value).lower()}.\n")
                            elif isinstance(value, str):
                                f.write(f"    {key} = '{value}'\n")
                            else:
                                f.write(f"    {key} = {value}\n")
                        f.write("/\n\n")

                        # Write kim_species sections
                        for species_data in kim_species_sections:
                            f.write("&kim_species\n")
                            for key, value in species_data.items():
                                f.write(f"    {key} = {value}\n")
                            f.write("/\n\n")

                    self.output_text.insert(
                        tk.END, "Updated KIM_config.nml with GUI configuration\n"
                    )
                    self.output_text.see(tk.END)
                except Exception as e:
                    raise Exception(f"Failed to write namelist file: {str(e)}")
            else:
                # Check if existing config file exists
                config_path = os.path.join(runpath, "KIM_config.nml")
                if not os.path.exists(config_path):
                    raise Exception(
                        f"No existing KIM_config.nml found at {config_path}. Please create one first or enable 'Use GUI Configuration'."
                    )
                self.output_text.insert(tk.END, "Using existing KIM_config.nml file\n")
                self.output_text.see(tk.END)

            # Initialize KIMpy
            self.kim_runner = KIMpy(runpath)

            # Start KIM in a separate thread
            self.running = True
            self.run_btn.config(state="disabled")
            self.stop_btn.config(state="normal")
            self.progress.start()

            self.output_text.insert(tk.END, "Starting KIM run...\n")
            self.output_text.see(tk.END)

            self.run_thread = threading.Thread(target=self.run_kim_thread)
            self.run_thread.start()

        except Exception as e:
            messagebox.showerror("Error", f"Failed to start KIM: {str(e)}")
            self.running = False
            self.run_btn.config(state="normal")
            self.stop_btn.config(state="disabled")
            self.progress.stop()

    def run_kim_thread(self):
        try:
            self.kim_runner.run()
            self.root.after(0, self.run_completed)
        except Exception as e:
            self.root.after(0, lambda: self.run_failed(str(e)))

    def run_completed(self):
        self.running = False
        self.run_btn.config(state="normal")
        self.stop_btn.config(state="disabled")
        self.progress.stop()
        self.output_text.insert(tk.END, "KIM run completed successfully!\n")
        self.output_text.see(tk.END)
        messagebox.showinfo("Success", "KIM run completed!")

    def run_failed(self, error_msg):
        self.running = False
        self.run_btn.config(state="normal")
        self.stop_btn.config(state="disabled")
        self.progress.stop()
        self.output_text.insert(tk.END, f"KIM run failed: {error_msg}\n")
        self.output_text.see(tk.END)
        messagebox.showerror("Error", f"KIM run failed: {error_msg}")

    def stop_kim(self):
        # This is a placeholder - actual implementation would need to handle process termination
        self.running = False
        self.run_btn.config(state="normal")
        self.stop_btn.config(state="disabled")
        self.progress.stop()
        self.output_text.insert(tk.END, "KIM run stopped.\n")
        self.output_text.see(tk.END)

    def plot_results(self):
        try:
            runpath = self.runpath.get()
            m_mode = int(self.config_vars["kim_setup"]["m_mode"].get())
            n_mode = int(self.config_vars["kim_setup"]["n_mode"].get())

            # Create a new window for plotting
            plot_window = tk.Toplevel(self.root)
            plot_window.title("KIM Results")
            plot_window.geometry("1000x700")

            # Create matplotlib figure with 3 subplots in 1 column
            fig, axes = plt.subplots(3, 1, figsize=(10, 12))
            fig.suptitle(f"KIM Results - m={m_mode}, n={n_mode}")

            # Try to load and plot results
            try:
                # Load phi_sol
                phi_sol = np.loadtxt(
                    os.path.join(runpath, f"out/m{m_mode}_n{n_mode}/fields/phi_sol.dat")
                )
                r = phi_sol[:, 0]
                phi = phi_sol[:, 1]

                # Plot Phi: real, imaginary, and absolute value
                axes[0].plot(r, phi.real, label="Real", linewidth=2)
                axes[0].plot(r, phi.imag, label="Imaginary", linewidth=2)
                axes[0].plot(r, np.abs(phi), label="Absolute", linewidth=2, linestyle="--")
                axes[0].set_xlabel("r [cm]")
                axes[0].set_ylabel(r"$\Phi$")
                axes[0].legend()

                # Load E_perp components
                E_perp = np.loadtxt(
                    os.path.join(runpath, f"out/m{m_mode}_n{n_mode}/fields/E_perp.dat")
                )
                E_perp_psi = np.loadtxt(
                    os.path.join(runpath, f"out/m{m_mode}_n{n_mode}/fields/E_perp_psi.dat")
                )
                E_perp_MA = np.loadtxt(
                    os.path.join(runpath, f"out/m{m_mode}_n{n_mode}/fields/E_perp_MA.dat")
                )

                # Plot combined E_perp components (absolute values only)
                axes[1].plot(E_perp[:, 0], np.abs(E_perp[:, 1]), label=r"$|E_\perp|$", linewidth=2)
                axes[1].plot(
                    E_perp_psi[:, 0],
                    np.abs(E_perp_psi[:, 1]),
                    label=r"$|E_\perp^{[\psi]}|$",
                    linewidth=2,
                )
                axes[1].plot(
                    E_perp_MA[:, 0], np.abs(E_perp_MA[:, 1]), label=r"$|E_\perp^{MA}|$", linewidth=2
                )
                axes[1].set_xlabel("r [cm]")
                axes[1].set_ylabel("Electric Field")
                axes[1].legend()

                # Load Br_pert
                Br_pert = np.loadtxt(
                    os.path.join(runpath, f"out/m{m_mode}_n{n_mode}/fields/Br_pert.dat")
                )

                # Plot Br_pert (absolute value)
                axes[2].plot(Br_pert[:, 0], np.abs(Br_pert[:, 1]), linewidth=2, color="red")
                axes[2].set_xlabel("r [cm]")
                axes[2].set_ylabel(r"$B_r$ [G]")

                list(map(lambda x: x.tick_params(top=True, right=True, which="both"), axes))
                list(map(lambda x: x.grid(True, ls=":"), axes))

            except FileNotFoundError as e:
                # Show error message in plot
                axes[0].text(
                    0.5,
                    0.5,
                    f"Results not found:\n{str(e)}",
                    ha="center",
                    va="center",
                    transform=axes[0].transAxes,
                    fontsize=12,
                )
                axes[1].axis("off")
                axes[2].axis("off")

            plt.tight_layout()

            # Create frame for matplotlib widgets
            plot_frame = ttk.Frame(plot_window)
            plot_frame.pack(fill="both", expand=True)

            # Embed plot in tkinter window
            canvas = FigureCanvasTkAgg(fig, plot_frame)
            canvas.draw()
            canvas.get_tk_widget().pack(fill="both", expand=True)

            # Add navigation toolbar for zooming and panning
            toolbar = NavigationToolbar2Tk(canvas, plot_frame)
            toolbar.update()

        except Exception as e:
            messagebox.showerror("Error", f"Failed to plot results: {str(e)}")

    def open_output_dir(self):
        runpath = self.runpath.get()
        output_dir = os.path.join(runpath, "out")
        if os.path.exists(output_dir):
            if sys.platform == "win32":
                os.startfile(output_dir)
            elif sys.platform == "darwin":
                subprocess.call(["open", output_dir])
            else:
                subprocess.call(["xdg-open", output_dir])
        else:
            messagebox.showwarning("Warning", "Output directory does not exist yet.")

    def setup_window(self):
        # Update all widgets to ensure proper sizing
        self.root.update_idletasks()

        # Set window size to fit content with some padding
        # Start with the run tab selected (first tab)
        self.notebook.select(self.run_frame)
        self.root.update_idletasks()

        # Set a reasonable window size that fits the two-column layout
        window_width = 900  # Enough for two columns with comfortable spacing
        window_height = 700  # Enough for all configuration sections

        # Center the window on screen
        screen_width = self.root.winfo_screenwidth()
        screen_height = self.root.winfo_screenheight()
        x = (screen_width - window_width) // 2
        y = (screen_height - window_height) // 2

        # Set geometry and bring to front
        self.root.geometry(f"{window_width}x{window_height}+{x}+{y}")
        self.root.lift()
        self.root.attributes("-topmost", True)
        self.root.after(100, lambda: self.root.attributes("-topmost", False))
        self.root.focus_force()


def main():
    root = tk.Tk()
    app = KIMConfigGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
