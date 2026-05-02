#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IAK - Intelligent Automated Kluster-generator Pipeline
Thin launcher — the application lives in the ``iak`` package.
"""

from iak.app import IAKApp

if __name__ == "__main__":
    IAKApp().mainloop()
