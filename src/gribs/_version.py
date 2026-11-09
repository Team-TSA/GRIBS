"""Program version (A-6 unit 10, U-39 (a)). The single definition of PROGRAM_VERSION.

This module imports nothing, so that `import gribs` (which reads __version__ from here) does not import numpy;
gribs.constants imports PROGRAM_VERSION from here. The pyproject version is a placeholder and is never used in
any output (ADR_A6_PACKAGE_SPLIT).
"""
PROGRAM_VERSION = "0.6.2-alpha"
