# Security Policy

## Supported versions

GRIBS is currently alpha software. Security fixes are applied to the latest development line when practical.

## Reporting a vulnerability

Please do not disclose a suspected security vulnerability in a public Issue.

Report the vulnerability privately to the repository maintainers through GitHub's private vulnerability reporting feature, if enabled. If private reporting is unavailable, contact the maintainers through a private channel listed by the Team-TSA organization.

Include:

- affected version or commit
- operating system and Python version
- reproduction steps
- potential impact
- suggested mitigation, if available

Do not include passwords, personal access tokens, SSH private keys, proprietary propellant data, export-controlled information, or third-party software that cannot be redistributed.

## Scope notes

GRIBS invokes an external executable when the `cea2` backend is selected. Users must obtain the executable and databases from a trusted source, verify file integrity, and understand that executing an untrusted binary can compromise the host system.

Configuration files should also be treated as untrusted input. Review JSON files before use, especially paths selecting executables, data directories, cache directories, and output directories.
