# Third-party components

PLIQ bundles third-party code under `pliq/_vendor/`. This repository’s **MIT
license** applies to the PLIQ Python wrapper and integration code only.

## TM-score

- **Path:** `pliq/_vendor/tmscore/TMscore.cpp`
- **Origin:** Zhang group TM-score distribution
- **Usage:** Compiled locally into the user cache; not redistributed as a binary
- **Action:** Comply with the Zhang group’s terms for TM-score when publishing
  or distributing derivatives

## BINANA

- **Path:** `pliq/_vendor/binana_python/`
- **Origin:** BINANA Python implementation (bundled for interaction fingerprinting)
- **Usage:** Invoked via Open Babel–prepared PDBQT inputs

Check the upstream tree for license text if you republish or modify BINANA.

## Other dependencies (PyPI)

Installed via `pip` — see each project’s license:

- RDKit, ProDy, PoseBusters, pandas, numpy, openpyxl, tqdm

## Open Babel & OTMol

Not bundled. Must be installed separately on the host system.
