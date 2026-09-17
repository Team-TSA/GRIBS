# Public release checklist

- [ ] Replace `YOUR-ORGANIZATION` in README and CITATION metadata.
- [ ] Confirm the repository name and description.
- [ ] Confirm copyright ownership and the MIT License choice.
- [ ] Keep the internal technical manual outside the public repository.
- [ ] Review the complete Git history for secrets and restricted material.
- [ ] Run `python GRIBS.py --selftest` and `pytest -q` in a clean environment.
- [ ] Check Python versions supported by CI.
- [ ] Confirm that default values are clearly identified as examples.
- [ ] Enable branch protection for `main` and require passing checks.
- [ ] Enable Dependabot alerts, secret scanning, push protection, and private vulnerability reporting where available.
- [ ] Create a signed or annotated `v0.1.0-alpha` tag and a GitHub Release.
- [ ] Recheck public pages while signed out of GitHub.
