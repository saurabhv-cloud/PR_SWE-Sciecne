# Engineered library snapshots

Each folder `<TASK_ID>/<library>-<upstream-short-sha>/` holds exactly the library files a task's Docker
build copies into the agent checkout, at their checkout-relative paths. Docker builds fetch this repository
(`git fetch --depth 1 <repo> main`), check that `git rev-parse HEAD:<folder>` equals the tree hash below,
copy the folder's files, and drop `.git`. Nothing else from this repository is read by a build.

| Task | Folder | Upstream | Engineering delta from upstream | Folder tree hash |
|---|---|---|---|---|
| SWS_009 | `SWS_009/proteoclade-342720b` | https://github.com/HeldLab/ProteoClade @ `342720bce1d4532cad0dca7cae7468e9c17ebe6c` (GPL v3) | Package files `proteoclade/*.py` and `LICENSE` only. Three comments reworded (`pcannotate.py` two, `pcutilities.py` one); line endings normalised to LF and trailing blank lines trimmed in `__init__.py`, `pcconstants.py`, `pcdb.py`. No behaviour change. | `9eb9e8e684bbcaefc6e7ebe032d6ba465147b927` |
