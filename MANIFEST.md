# Engineered library snapshots

Each folder `<TASK_ID>/<library>-<upstream-short-sha>/` holds exactly the library files a task's Docker
build copies into the agent checkout, at their checkout-relative paths. Docker builds fetch this repository
(`git fetch --depth 1 <repo> main`), check that `git rev-parse HEAD:<folder>` equals the tree hash below,
copy the folder's files, and drop `.git`. Nothing else from this repository is read by a build.

| Task | Folder | Upstream | Engineering delta from upstream | Folder tree hash |
|---|---|---|---|---|
| SWS_009 | `SWS_009/proteoclade-342720b` | https://github.com/HeldLab/ProteoClade @ `342720bce1d4532cad0dca7cae7468e9c17ebe6c` (GPL v3) | Package files `proteoclade/*.py` and `LICENSE` only. Three comments reworded (`pcannotate.py` two, `pcutilities.py` one); line endings normalised to LF and trailing blank lines trimmed in `__init__.py`, `pcconstants.py`, `pcdb.py`. No behaviour change. | `9eb9e8e684bbcaefc6e7ebe032d6ba465147b927` |
| SWS_006 | `SWS_006/gemmi-5cc1c23` | https://github.com/project-gemmi/gemmi @ `5cc1c23c6007e0e6cbd69289c6f7c0bff50e943e` (MPL-2.0) | No upstream file is modified or removed. Adds a task-authored Python package `study/` (`__init__.py`, `pipeline.py`, `model.py`, `structure.py`, `reporting.py`), a methods document `study/METHODS.md`, and two data assets under `study/assets/` (COD 2023162 geometry CIF; synthetic ensemble `ensemble.npz`). The upstream `LICENSE.txt` is included unmodified. Everything else in the checkout is stock gemmi, fetched from upstream. | `fa7ecf663ce85ec7a36ed0d518463f516349d85c` |
| SWS_008 | `SWS_008/study-232b605` | https://github.com/nextflow-io/nextflow @ `232b60569865e9a4577e48c1955409238359d6ca` (Apache-2.0; docs CC-BY-SA-4.0) | The Nextflow tree itself is byte-identical to upstream and is **not** in this folder: fetch it from upstream. This folder holds only the task-authored study material at study-relative paths: `analysis/` (Python package `hicparams`, `pyproject.toml`) and `protocol/` (workflow `main.nf`, `nextflow.config`, `assay_settings.json`, `deliver.sh`, `bin/*.py`, `README.md`). Recorded run outputs (`delivery/`) are task data and stay in the task package. | `36b24a43ea7a31035a28fe80a1db853031c94bf1` |

No row for SWS_000 (installs released versions with pip) or SWS_004 (native defects at a pinned upstream commit, so the build fetches upstream and deletes a fixed path list). SWS_001 to SWS_003, SWS_005 and SWS_007 were not audited for this repository.
