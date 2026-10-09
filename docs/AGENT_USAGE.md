# Use VASPBERRY with an agent

VASPBERRY includes a Codex plugin and a standalone skill for agents that can
read local skill files and run commands. You can ask for installation in
ordinary language:

> Install the VASPBERRY agent plugin from https://github.com/Infant83/VASPBERRY
> and set it up for use in this project.

In Korean:

> https://github.com/Infant83/VASPBERRY 의 VASPBERRY를 내 에이전트에
> 설치하고 사용할 수 있게 설정해줘.

Give the repository URL on the first request so the agent can find the
maintained installation instructions. The agent needs local file/terminal
access and network access for downloads; a text-only chat cannot install
software on your computer. Installing the workflow does not require VASP,
a Fortran compiler or an OpenAI API key. Native calculations have the
dependencies described in [BUILD.md](BUILD.md).

## Codex plugin

The repository provides the `vaspberry` marketplace and its `vaspberry`
plugin. In a Codex environment with the `codex` CLI, the agent can run:

```bash
GIT_LFS_SKIP_SMUDGE=1 codex plugin marketplace add Infant83/VASPBERRY --ref master --sparse .agents/plugins --sparse plugins/vaspberry --json
codex plugin add vaspberry@vaspberry --json
codex plugin list --marketplace vaspberry --json
```

This uses the repository's default branch, `master`, for the plugin
workflow. The calculation engine is selected separately from stable releases.
The sparse checkout downloads the plugin and catalog without checking out
the large optional material example inputs. No external service connection
or separate account authentication is configured by this skills-only plugin.

Start a new Codex conversation after installation. If the CLI is unavailable,
the agent can install the standalone skill below. These GitHub installation
routes do not require public plugin-directory publication. A public-directory
listing is a separate submission and review process.

## Standalone skill

For Codex, local skills are stored in `~/.agents/skills/<skill-name>/` for
user-wide use, or `<project>/.agents/skills/<skill-name>/` for one project.
The maintained source is
[`plugins/vaspberry/skills/vaspberry/`](../plugins/vaspberry/skills/vaspberry/).
Other agents must support the `SKILL.md` format and use their own documented
skill location; installing a Codex skill does not itself establish support
in every agent.

From an existing VASPBERRY checkout, use Python 3.10+:

```bash
python3 tools/install_agent_skill.py
```

The installer copies the complete skill, including its references, helper
and MIT notice, to `~/.agents/skills/vaspberry`. It refuses to overwrite an
existing destination and prints a JSON receipt with the copied file hashes.
To install for one project or into another agent's supported skill directory:

```bash
python3 tools/install_agent_skill.py --destination /path/to/project/.agents/skills/vaspberry
```

You can preview the paths without installing:

```bash
python3 tools/install_agent_skill.py --dry-run
```

If you only need the skill, the agent can fetch a small checkout into a new
directory. These command-specific Git settings let this work without a
Git LFS executable and leave optional large inputs as pointer files:

```bash
git -c filter.lfs.process= -c filter.lfs.clean= -c filter.lfs.smudge= -c filter.lfs.required=false clone --depth 1 --filter=blob:none --sparse --branch master https://github.com/Infant83/VASPBERRY.git VASPBERRY-agent
git -C VASPBERRY-agent -c filter.lfs.process= -c filter.lfs.clean= -c filter.lfs.smudge= -c filter.lfs.required=false sparse-checkout set plugins/vaspberry tools
python3 VASPBERRY-agent/tools/install_agent_skill.py
```

Alternatively, copy the entire `plugins/vaspberry/skills/vaspberry` directory
to the chosen skill location. Do not copy only `SKILL.md`: it links to the
supporting references and helper. Use either the plugin or the standalone
skill in a given agent to avoid duplicate copies. Start a new conversation;
restart Codex if the new skill is not detected.

## First use and scientific results

After installation, try:

- “Install the latest stable VASPBERRY and run its public analytic demo.
  Explain the Chern numbers and include citations.”
- “Check whether the VASP outputs in this directory support a Chern-number
  calculation, then run the supported method with numerical checks.”
- “Use VASPBERRY to plot these charge Hall results. Explain the units and
  report which mesh and band-window convergence checks are still missing.”
- “Write a Methods paragraph for this VASPBERRY calculation using the
  recorded inputs, parameters, software version and method references.”

The skill reads the selected engine's own documentation, diagnoses available
dependencies and checks the input requirements for the requested method.
It preserves explicitly requested versions and existing user checkouts.
When an engine installation is needed, it resolves the latest stable release,
creates a separate versioned source/environment, verifies the source commit
and runs the checks appropriate to that release. Native WAVECAR calculations
need the compatible compiler and numerical libraries in [BUILD.md](BUILD.md).

The included analytic demo is a software/model check. It does not validate
a real material or establish convergence. Real-material reports must describe
input compatibility, numerical diagnostics, convergence evidence, physical
assumptions and limitations. See the skill's
[scientific validation](../plugins/vaspberry/skills/vaspberry/references/scientific-validation.md)
and [reporting](../plugins/vaspberry/skills/vaspberry/references/reporting.md)
guides. Missing input or validation should be reported explicitly.

## Updates and reproducibility

Ask “Check for a VASPBERRY update and install the latest stable release.”
The skill looks up the current engine release when requested and records
the exact version and commit. Earlier results keep their recorded source.
There is no background updater or automatic replacement of an engine in use.

Refresh the GitHub plugin marketplace, then check the installed plugin with:

```bash
codex plugin marketplace upgrade vaspberry --json
codex plugin list --marketplace vaspberry --json
```

Start a new conversation after an update. The engine and plugin have separate
version numbers. A new engine release alone does not require copying the
engine into the plugin; update the workflow when commands, input requirements
or guidance change. For a standalone skill, preserve the existing copy and
install the updated skill to a new destination before replacing the active
copy through your agent's usual installation workflow.

## Credit

The agent follows the author's [README citation guidance](../README.md#citation)
and the skill's [citation notes](../plugins/vaspberry/skills/vaspberry/references/citations.md).
Results record the actual VASPBERRY version/commit and GitHub source, the
author's applicable PRB/PRL papers, and the numerical-method references used.
Citation exports are available through the helper. An old DOI must not be
presented as the identifier of a newer version. The independent plugin and
skill files use MIT; the engine's existing contributor and WAVETRANS notices
remain applicable.

## Distribution files and official host guidance

- [Repository catalog](../.agents/plugins/marketplace.json)
- [Codex manifest](../plugins/vaspberry/.codex-plugin/plugin.json)
- [Skill instructions](../plugins/vaspberry/skills/vaspberry/SKILL.md)
- [Plugin overview](../plugins/vaspberry/README.md)
- [OpenAI plugin packaging and marketplaces](https://developers.openai.com/plugins/build/plugins)
- [Codex plugin CLI](https://learn.chatgpt.com/docs/cli/reference)
- [Codex local skill locations](https://learn.chatgpt.com/docs/build-skills#where-codex-loads-local-skills)
