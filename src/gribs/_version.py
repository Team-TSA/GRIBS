"""Program version (A-6 unit 10, U-39 (a)). The single definition of PROGRAM_VERSION.

This module imports nothing, so that `import gribs` (which reads __version__ from here) does not import numpy;
gribs.constants imports PROGRAM_VERSION from here. The pyproject version is the PEP 440 form of the release version and is never used in
any output (ADR_A6_PACKAGE_SPLIT).
"""
PROGRAM_VERSION = "0.7.0-alpha"
