# Author: Kanak Varma
"""
Thin shim — delegates to tools/run_all.py so you can run the pipeline from
the project root without changing directory.

    python run_all.py bengaluru
    python run_all.py chennai --steps core,maps
    python run_all.py --list
"""
import os
import subprocess
import sys

tools_runner = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tools", "run_all.py")
sys.exit(subprocess.run([sys.executable, tools_runner] + sys.argv[1:]).returncode)

