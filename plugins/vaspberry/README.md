# VASPBERRY for Codex

Codex skill plugin for VASPBERRY topology and response analysis. It guides engine installation, input selection, actual released commands, scientific validation and source/method credit. A standard-library Python helper diagnoses an engine checkout, generates exact-version credit, and runs a small released analytic model demonstration.

The plugin version is **0.1.3**; its validated engine is the **2026-10-09 help amendment of VASPBERRY v1.6.6**, commit `2e7067b8b8b6d3ab444eef5bc025f4ce5720ca4d`. For a new installation the skill checks GitHub for the latest stable release, downloads that release separately, reads its own instructions and records the exact version/commit used. It preserves explicitly requested versions and existing user checkouts. No engine binaries, private research input, VASP source, credentials, MCP server or network service are bundled.

Use a Python 3.10+ interpreter for the following commands. From this plugin directory, try the helper directly:

```bash
python3 skills/vaspberry/scripts/vaspberry_agent.py doctor --source /path/to/VASPBERRY
python3 skills/vaspberry/scripts/vaspberry_agent.py demo --source /path/to/VASPBERRY --output /new/demo-directory
python3 skills/vaspberry/scripts/vaspberry_agent.py cite --source /path/to/VASPBERRY --output /new/citations-directory --method fhs
```

Source discovery uses `--source`, then `VASPBERRY_ROOT`, then the working directory. Demo dependencies and exact-commit reproduction commands are in [setup.md](skills/vaspberry/references/setup.md). Every output directory must be new. The analytic demo writes its engine outputs under `calculation/`, stdout/stderr logs, `receipt.json`, and citation/provenance files. It is a software demonstration, not a material prediction.

Example requests after installing the plugin in Codex:

- “Use VASPBERRY to check whether these VASP outputs support a Chern-number calculation, then calculate it with validation and citations.”
- “Run the VASPBERRY analytic demo and explain the Chern values and the limits of this test.”
- “Use VASPBERRY to plot my charge Hall calculation and report which mesh and band-window checks are still missing.”

## Scientific reporting

The [scientific validation guide](skills/vaspberry/references/scientific-validation.md) defines observables and abbreviations, physical input conditions, numerical checks, convergence and interpretation. The [reporting guide](skills/vaspberry/references/reporting.md) gives Methods/Results and caption guidance with explicit evidence requirements. Real-material results require their own supported inputs and validation; the public analytic demonstration is a bounded software/model test.

## Installation and sharing

This plugin is distributed from `plugins/vaspberry/` in the canonical
[VASPBERRY repository](https://github.com/Infant83/VASPBERRY). The manifest is
`.codex-plugin/plugin.json`; the complete standalone skill is
`skills/vaspberry/`. See the [agent installation and usage guide](https://github.com/Infant83/VASPBERRY/blob/master/docs/AGENT_USAGE.md)
for first-run requests, engine setup and reproducible analysis.

With the Codex CLI, add the repository catalog using a sparse checkout and
install the plugin:

```bash
GIT_LFS_SKIP_SMUDGE=1 codex plugin marketplace add Infant83/VASPBERRY --ref master --sparse .agents/plugins --sparse plugins/vaspberry
codex plugin add vaspberry@vaspberry
```

The sparse checkout fetches the catalog and plugin; the separately installed
engine and large example inputs are not needed to install the plugin. The
catalog name is `vaspberry`. This repository catalog is not an OpenAI public
directory listing.

Without the Codex CLI, use the repository's
[local skill installer](https://github.com/Infant83/VASPBERRY/blob/master/tools/install_agent_skill.py)
from a repository checkout:

```bash
python3 tools/install_agent_skill.py
```

It copies the **entire** `plugins/vaspberry/skills/vaspberry/` directory,
including scripts, references, agents metadata and LICENSE, into
`~/.agents/skills/vaspberry/`. It refuses an existing destination so a previous
installation is preserved. A manual copy must preserve the same complete
folder and must not merge or overwrite an existing skill without choosing an
update.

Start a **new Codex conversation** after either installation, then request:
“Use VASPBERRY to install the latest stable engine, run its public analytic
demonstration, and give me the result and citations.” The plugin and standalone
skill provide the same workflow; choose one installation route to avoid
loading duplicate skills. The engine is installed separately using
[setup.md](skills/vaspberry/references/setup.md).

## Engine updates

VASPBERRY source and releases remain in the canonical GitHub repository. The plugin does not carry a second copy. Ask “Install the latest stable VASPBERRY” or “Check for a VASPBERRY update”; the skill resolves the current release and sets up a separate versioned checkout/environment as needed. Existing analyses keep their recorded version. There is no background polling or unattended replacement of an engine in use.

An engine release does not by itself require repackaging the plugin. Update the plugin when its workflow, helper, supported input/command assumptions or citation handling need to change. New releases are checked against their own documentation and validations; the older plugin baseline is not evidence of their compatibility.

## Credit and licensing

The MIT [LICENSE](LICENSE) covers the independently authored files in this standalone plugin. Preserve the engine's WAVETRANS attribution to R. M. Feenstra and M. Widom and other contributor/third-party notices. The VASPBERRY engine is separately installed; the plugin license does not relicense engine or third-party code.

The citation helper uses the actual VERSION/CITATION.cff and Git identity, and follows the author's README scholarly citation guidance; see [citation notes](skills/vaspberry/references/citations.md). Software credit uses the GitHub source/version, while the author's PRB/PRL papers and applicable numerical-method references have separate roles. Zenodo is optional. The 2018 Zenodo DOI for VASPBERRY 1.0 is excluded from later-version credit. Other scientific methods beyond the bundled FHS/FH/background references must be selected from the engine's technical report and verified for the actual calculation.
