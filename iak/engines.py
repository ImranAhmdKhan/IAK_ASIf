from __future__ import annotations

import os
import shutil
import subprocess
import sys


ENGINE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "iak_engine")

XTB_URLS = [
    "https://github.com/grimme-lab/xtb/releases/download/v6.7.1/xtb-6.7.1-linux-x86_64.tar.xz",
    "https://github.com/grimme-lab/xtb/releases/download/v6.6.1/xtb-6.6.1-linux-x86_64.tar.xz",
]

CREST_URLS = [
    "https://github.com/grimme-lab/crest/releases/download/v3.0.2/crest-x86_64-unknown-linux-gnu.tar.xz",
    "https://github.com/grimme-lab/crest/releases/download/v3.0.2/crest-x86_64-pc-linux-gnu.tar.xz",
    "https://github.com/grimme-lab/crest/releases/download/v3.0.2/crest-linux-x86_64.tar.xz",
    "https://github.com/grimme-lab/crest/releases/download/v3.0.2/crest.tar.xz",
    "https://github.com/grimme-lab/crest/releases/download/v3.0.1/crest-x86_64-unknown-linux-gnu.tar.xz",
    "https://github.com/grimme-lab/crest/releases/download/v2.12/crest.tar.xz",
    "https://github.com/grimme-lab/crest/releases/download/v2.11.2/crest.tar.xz",
]

XTB_DIR = None
CREST_DIR = None
ORCA_DIR = None
ORCA_IS_WINDOWS = False
_WSL_ORCA_EXISTS = None


def is_tool_available(tool_name):
    global _WSL_ORCA_EXISTS
    if tool_name == "xtb":
        return XTB_DIR is not None or shutil.which("xtb") is not None
    if tool_name == "crest":
        return CREST_DIR is not None or shutil.which("crest") is not None
    if tool_name == "orca":
        if ORCA_DIR is not None or shutil.which("orca") is not None or shutil.which("orca.exe") is not None:
            return True
        if sys.platform == "win32":
            if _WSL_ORCA_EXISTS is None:
                try:
                    proc = subprocess.run(
                        "wsl -e bash -c 'export PATH=\"/usr/bin:/bin:/usr/local/bin:$PATH\"; which orca'",
                        shell=True,
                        capture_output=True,
                        timeout=3,
                    )
                    _WSL_ORCA_EXISTS = proc.returncode == 0
                except Exception:
                    _WSL_ORCA_EXISTS = False
            return _WSL_ORCA_EXISTS
        return False
    return False


def get_wsl_path(win_path):
    drive, tail = os.path.splitdrive(win_path)
    if drive:
        return f"/mnt/{drive[0].lower()}{tail.replace(os.sep, '/')}"
    return win_path.replace(os.sep, "/")


def inject_embedded_engines():
    global XTB_DIR, CREST_DIR, ORCA_DIR, ORCA_IS_WINDOWS
    if not os.path.exists(ENGINE_DIR):
        return
    for root, dirs, files in os.walk(ENGINE_DIR):
        try:
            if "xtb" in files and os.path.basename(root) == "bin" and os.path.isfile(os.path.join(root, "xtb")):
                XTB_DIR = root
                if sys.platform == "win32":
                    subprocess.run(f"wsl chmod +x \"{get_wsl_path(os.path.join(root, 'xtb'))}\"", shell=True, capture_output=True)

            if "crest" in files and CREST_DIR is None and os.path.isfile(os.path.join(root, "crest")):
                CREST_DIR = root
                if sys.platform == "win32":
                    subprocess.run(f"wsl chmod +x \"{get_wsl_path(os.path.join(root, 'crest'))}\"", shell=True, capture_output=True)

            if ORCA_DIR is None:
                if "orca.exe" in files and os.path.isfile(os.path.join(root, "orca.exe")):
                    ORCA_DIR = root
                    ORCA_IS_WINDOWS = True
                elif "orca" in files and os.path.isfile(os.path.join(root, "orca")):
                    ORCA_DIR = root
                    ORCA_IS_WINDOWS = False
                    if sys.platform == "win32":
                        subprocess.run(f"wsl chmod +x \"{get_wsl_path(os.path.join(root, 'orca'))}\"", shell=True, capture_output=True)
                        subprocess.run(f"wsl chmod +x \"{get_wsl_path(root)}/orca_\"* 2>/dev/null", shell=True, capture_output=True)
        except OSError:
            continue



inject_embedded_engines()
