#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IAK - Intelligent Automated Kluster-generator Pipeline v11.2
The Job-Report and Visual Analytics Edition

New in v11.2:
- Per-job JSON/CSV/Markdown summary reports.
- High-level input correctness/incorrectness assessment for every job.
- GUI progress monitor with percent completion and ETA.
- JACS-grade plotting controls with custom colors, pan, zoom, and crop.
- Interactive molecular model viewer (Ball-and-stick, CPK, Wireframe, Licorice).
- User-selectable atom, bond, background, graph, axis, and grid colors.

Bug fixes in this version:
- ORCA base_method now uses word-based keyword filtering (fixes TightOpt Freq corruption).
- MPI crash detection uses specific error tokens (fixes false positive from OpenMPI header text).
- Tier 1 ORCA log now correctly reports Serial vs Parallel mode.

Run with: python IAK_ASIf.py
"""

from iak.app import IAKApp

if __name__ == "__main__":
    IAKApp().mainloop()
