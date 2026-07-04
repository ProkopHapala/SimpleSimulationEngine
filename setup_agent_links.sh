#!/usr/bin/env bash
# Setup symlinks so AI coding assistants can find shared skills & workflows.
# Canonical source: doc/AGENTs/skills/ and doc/AGENTs/workflows/
# Run this after cloning the repo on a new machine.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$REPO_ROOT"

CANONICAL_SKILLS="doc/AGENTs/skills"
CANONICAL_WORKFLOWS="doc/AGENTs/workflows"

# Tool config dirs that use skills/ and workflows/ subdirs
TOOL_DIRS=(.devin .windsurf .cursor .claude .codex)

link() {
    local target="$1" link_path="$2"
    if [[ -L "$link_path" ]]; then
        local current
        current="$(readlink "$link_path")"
        if [[ "$current" == "$target" ]]; then
            echo "  OK (exists): $link_path -> $target"
            return 0
        fi
        rm "$link_path"
    elif [[ -e "$link_path" ]]; then
        echo "  SKIP (not a symlink, won't overwrite): $link_path"
        return 0
    fi
    ln -s "$target" "$link_path"
    echo "  LINKED: $link_path -> $target"
}

echo "Setting up agent skills/workflows symlinks..."
echo "Repo root: $REPO_ROOT"
echo

if [[ ! -d "$CANONICAL_SKILLS" ]]; then
    echo "ERROR: $CANONICAL_SKILLS not found. Is this the right repo?" >&2
    exit 1
fi

for dir in "${TOOL_DIRS[@]}"; do
    echo "[$dir]"
    mkdir -p "$dir"
    link "../$CANONICAL_SKILLS" "$dir/skills"
    link "../$CANONICAL_WORKFLOWS" "$dir/workflows"
    echo
done

echo "Done."
