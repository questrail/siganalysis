# The Python versions listed in the pyproject.toml classifiers, which the
# test-all recipe runs the suite against.
python_versions := "3.12 3.13"

# List the available justfile recipes
[group('general')]
@default:
  just --list --unsorted

# List the lines of code in the project
[group('general')]
loc:
  scc --remap-unknown "-*- Justfile -*-":"justfile"

# Search pydoc for given term
[group('general')]
doc term:
  uv run python -m pydoc {{term}}

# Type check using ty
[group('test')]
check:
  uv run ty check

# Lint and format code without making changes
[group('test')]
lint:
  uv run ruff format --check
  uv run ruff check

# Lint and format code and apply changes
[group('test')]
fix:
  uv run ruff format
  uv run ruff check --fix

# Test code using pytest
[group('test')]
test *args:
  uv run pytest {{args}}

# Test code using pytest on every supported Python
[group('test')]
test-all:
  #!/usr/bin/env bash
  set -euo pipefail
  # There is no CI, so this is the only thing keeping the versions claimed in
  # the pyproject.toml classifiers honest.
  for v in {{python_versions}}; do
    echo "== Python $v =="
    uv run --isolated --python "$v" --extra plotting --with pytest pytest -q
  done

# Add dependency
[group('dependencies')]
add dep:
  uv add {{dep}}

# Add dependency to the development group
[group('dependencies')]
dev dep:
  uv add --dev {{dep}}

# Update dependency in the project dependencies or any group
[group('dependencies')]
up dep:
  uv lock -P {{dep}}
  uv sync

# List the outdated dependencies
[group('dependencies')]
out:
  uv pip list --outdated

# Lock/freeze dependencies
[group('dependencies')]
lock:
  uv lock

# Check, test, build, smoke test, publish to PyPI, and then tag
[group('deploy')]
deploy: release-check lint check test-all build smoke
  #!/usr/bin/env bash
  set -euo pipefail
  version="$(uv version --short)"
  # Prompt for the token rather than reading it from a store, so that nothing
  # extra has to be installed and the token stays out of the shell history.
  # Nothing echoes while it is pasted.
  read -rsp "PyPI token for v${version} (pypi-...): " token
  echo
  if [ -z "${token}" ]; then echo "No token given."; exit 1; fi
  case "${token}" in
    pypi-*) ;;
    *) echo "That does not look like a PyPI API token."; exit 1 ;;
  esac
  # Tag only after PyPI accepts the upload, so that a tag is evidence the
  # release shipped rather than evidence it was attempted.
  UV_PUBLISH_TOKEN="${token}" uv publish
  git tag -a "v${version}" -m "Release v${version}"
  git push origin "v${version}"
  echo "Published and tagged v${version}"

# Confirm the tree is ready to release the version in pyproject.toml
[group('deploy')]
release-check:
  #!/usr/bin/env bash
  set -euo pipefail
  version="$(uv version --short)"
  echo "Preparing v${version}"
  if [ -n "$(git status --porcelain)" ]; then echo "Working tree is dirty."; exit 1; fi
  if [ "$(git rev-parse --abbrev-ref HEAD)" != "master" ]; then echo "Not on master."; exit 1; fi
  git fetch --quiet origin
  if [ -n "$(git log origin/master..HEAD --oneline)" ]; then echo "Unpushed commits."; exit 1; fi
  if git rev-parse -q --verify "refs/tags/v${version}" >/dev/null; then
    echo "Tag v${version} already exists."; exit 1
  fi
  if ! grep -q "## v${version}" CHANGELOG.md; then
    echo "No CHANGELOG entry for v${version}."; exit 1
  fi
  # A clean tree rules out an untracked file, so existing here means committed.
  if [ ! -f "docs/releases/v${version}.md" ]; then
    echo "No release notes at docs/releases/v${version}.md."; exit 1
  fi
  echo "Ready."

# Build the sdist and wheel from a clean dist/
[group('deploy')]
build:
  # Clearing dist/ leaves uv publish exactly the current version to upload,
  # rather than every build the directory has ever collected.
  rm -rf dist
  uv build

# Install the built wheel on its own and check that it works
[group('deploy')]
smoke:
  #!/usr/bin/env bash
  set -euo pipefail
  # The tests all run against the source tree, so only this catches a wheel
  # missing a module or one that makes matplotlib a hard dependency.
  version="$(uv version --short)"
  wheel="./dist/siganalysis-${version}-py3-none-any.whl"
  run() { uv run --isolated --python 3.12 --no-project --with "$wheel" python -c "$1"; }
  run "import siganalysis, importlib.util; assert siganalysis.__version__ == '${version}', siganalysis.__version__; assert callable(siganalysis.stft); assert importlib.util.find_spec('matplotlib') is None, 'matplotlib is a hard dependency'; print('imports without matplotlib: ok')"
  # Reaching a plotting name without matplotlib has to fail, and say so.
  if run "import siganalysis; siganalysis.plot_spectrogram" 2>/dev/null; then
    echo "plot_spectrogram resolved without matplotlib installed."; exit 1
  fi
  echo "wheel ok"
