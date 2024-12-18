# List the available justfile recipes
[group('general')]
@default:
  just --list --unsorted

# List the lines of code in the project
[group('general')]
loc:
  scc --remap-unknown "-*- Justfile -*-":"justfile"

# Lint code using ruff
[group('test')]
lint: 
  ruff check

# Test code using nose2
[group('test')]
test: 
  uv run nose2

# Add/update dependency
[group('dependencies')]
add dep:
  uv add {{dep}}

# Add/update dependency to the development group
[group('dependencies')]
dev dep:
  uv add --dev {{dep}}

# List the outdated dependencies
[group('dependencies')]
out:
  uv pip list --outdated

# Lock/freeze dependencies
[group('dependencies')]
lock:
  uv lock
