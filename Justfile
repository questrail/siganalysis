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

# Check, test, build, and publish to PyPI
[group('deploy')]
deploy: lint check test
  rm -rf dist
  uv build
  uv publish
