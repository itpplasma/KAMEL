import numpy as np
import matplotlib.pyplot as plt

class MASTU_config:

    R0 = 80.0
    r_eff_wall = 95.0 # ideal wall effective radius
    r_eff_plasma = 85.0 # approximate effective plasma radius determined by separatrix
    delta_r_antenna = 3.0 # distance of RMP antenna from plasma boundary in cm
    r_eff_antenna = r_eff_plasma + delta_r_antenna # antenna effective radius default
    Btor = 6550.0 # approximate value to give order of magnitudeo
    I0_rmp = 6.0e12

    def __init__(self):
        pass

    def plot_1d_radial_config(self, save=False):

        regions = [
            ("Plasma", 0.0, self.r_eff_plasma, "tab:orange"),
            ("Vacuum", self.r_eff_plasma, self.r_eff_wall, "whitesmoke"),
            #("Vacuum", self.r_eff_antenna, self.r_eff_wall, "whitesmoke"),
        ]
        points = [
            ("RMP Antenna", self.r_eff_antenna, "tab:red", 3),
            ("Wall", self.r_eff_wall, 'tab:gray', 5)
        ]
        major = [
            ("Major radius", 0.0, self.R0, 'k')
        ]

        fig, ax = plt.subplots(figsize=(6, 2))
        height = 1
        for name, pos, color, lw in points:
            ax.plot([pos, pos], [-height/2,height/2], color=color, label=name, lw=lw)
        for name, start, end, color in regions:
            ax.barh(0, width=end-start, height=height, left=start, color=color, edgecolor=color, label=name)

        # Formatting
        ax.set_xlim(0, self.r_eff_wall+0.1)
        ax.set_ylim(-height/2, height/2)
        #ax.set_yticks([])
        ax.set_xlabel(r"Effective radius $r_{\rm{eff}}$ [cm]")
        ax.get_yaxis().set_visible(False)

        xticks = np.arange(0, self.r_eff_wall + 0.1, 10)  # Increase the number of x-ticks with a step of 5
        ax.set_xticks(xticks)
        ax.set_xticklabels([f"{tick}" for tick in xticks])
        
        ax.legend(loc="center", bbox_to_anchor=(0.5, 1.25), ncol=len(regions) + len(points))
        #plt.title("MAST-U 1d radial configuration")
        plt.tight_layout()

        if save:
            plt.savefig('MASTU_radial_config.pdf', bbox_inches='tight', transparent=False)
        
        plt.show()
