#!/usr/bin/env python3
"""
Standalone launcher for KIM GUI that can be run from anywhere
"""

import os
import sys


def main():
    # Get the absolute path of this script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Add the KIMgui directory to Python path if not already there
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)

    # Import and run the GUI
    try:
        from kim_gui import main as gui_main

        gui_main()
    except ImportError as e:
        print(f"Error importing KIM GUI: {e}")
        print("Please ensure you are running this from the KAMEL/python/KIMgui directory")
        print("or that the KIMgui module is properly installed.")
        sys.exit(1)


if __name__ == "__main__":
    main()
