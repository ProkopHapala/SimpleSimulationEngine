#!/usr/bin/env python3
"""
Hybrid heuristic + LLM-assisted classification of .md files.

Purpose:
  Automatically classify markdown documents by type, quality, and confidence,
  and optionally rename git-tracked files to a category-based naming scheme.
  This makes document purpose immediately visible from the filename and reduces
  ambiguity in a large knowledge-base repository.

Category → extension mapping (renaming only applies to selected categories):
  user_guide            → .guide.md     (manuals, how-to, tutorials)
  implementation_report → .report.md    (polished reports with validation)
  technical_note        → .note.md      (math derivations, algorithm notes)
  progress_report       → .progress.md  (debug logs, session notes, dated entries)
  draft_note            → .plan.md      (WIP, exploratory, temporary docs)
  llm_chat              → .chat.md      (LLM conversation transcripts)

Categories never renamed: meta_doc, config, api_doc, design_doc

Rename rules:
  1. Only git-tracked files are renamed (untracked files are skipped).
  2. Excluded paths (never renamed): doc/topical_audit/, doc/Markdown/cpp/
  3. Excluded patterns (never renamed): auto-generated source docs (.h.md, .cpp.md, etc.)
  4. Trailing ext words are stripped before adding the new extension to avoid
     duplication (e.g. Foo_report.md → Foo.report.md, not Foo_report.report.md).
  5. Idempotent: files already having .{ext}.md are never double-renamed.
  6. Dry-run mode is available to preview changes before executing.
  7. --backup copies tracked .md files to a backup directory before any rename.
  8. --commit auto-commits the rename after executing.

Quality levels:
  high    - Polished, comprehensive, well-structured
  medium  - Useful but informal or incomplete
  low     - Very brief, disorganized, temporary, or misclassified

Confidence: 0-100 heuristic score (higher = more certain)
"""

import os, re, sys, math, shutil, subprocess, argparse
from pathlib import Path
from typing import Dict, List, Tuple, Set

REPO = "/home/prokophapala/git/SimpleSimulationEngine"

# ─────────────────────────────────────────────
# EXCLUSIONS (these files are classified but never renamed)
# ─────────────────────────────────────────────

# Folders excluded from renaming (relative to REPO)
EXCLUDED_RENAME_PATHS = [
    "doc/topical_audit/",
    "doc/Markdown/cpp/",
]

# Auto-generated source-to-markdown suffixes to skip renaming
# e.g. Foo.h.md, Foo.cpp.md, Foo.py.md
EXCLUDED_RENAME_SUFFIXES = {".h", ".cpp", ".c", ".py", ".f90", ".js", ".cl", ".hpp"}


def is_excluded_from_rename(filepath: str, repo_path: str = REPO) -> bool:
    """Return True if file should not be renamed (excluded folder or auto-generated source doc)."""
    rel = os.path.relpath(filepath, repo_path)
    # Path exclusions
    for excl in EXCLUDED_RENAME_PATHS:
        if rel.startswith(excl):
            return True
    # Pattern exclusions: stem ends with a source-code suffix
    basename = os.path.basename(filepath)
    stem, ext = os.path.splitext(basename)  # ext = .md
    for src_ext in EXCLUDED_RENAME_SUFFIXES:
        if stem.lower().endswith(src_ext.lower()):
            return True
    return False

# ─────────────────────────────────────────────
# QUALITY THRESHOLDS
# ─────────────────────────────────────────────
QUALITY_HIGH_MIN_SECTIONS = 6
QUALITY_HIGH_MIN_LENGTH = 2000
QUALITY_LOW_MAX_SECTIONS = 2
QUALITY_LOW_MAX_LENGTH = 1500
QUALITY_LOW_MAX_MARKERS = 4

# ─────────────────────────────────────────────
# REGEX PATTERNS
# ─────────────────────────────────────────────

# LLM model name patterns (headings level 1-3)
LLM_HEADING_RE = re.compile(
    r'^#{1,3}\s*('
    r'USER'
    r'|Gemini'
    r'|ChatGPT[^\w]'
    r'|Claude[^\w]'
    r'|DeepSeek[^\w]'
    r'|Kimi[^\w]'
    r'|SWE-[^\w]'
    r'|GPT-[^\w]'
    r'|Devin[^\w]'
    r'|Copilot[^\w]'
    r'|Post-mortem'
    r'|Cody[^\w]'
    r'|Tabby[^\w]'
    r'|OpenAI[^\w]'
    r'|gpt-[^\w]'
    r'|o1[^\w]'
    r'|o3[^\w]'
    r'|llama[^\w]'
    r'|mistral[^\w]'
    r'|codestral[^\w]'
    r'|perplexity[^\w]'
    r'|cursor[^\w]'
    r'|sweep-[^\w]'
    r')',
    re.MULTILINE | re.IGNORECASE
)

# Share links for LLM chat transcripts
SHARE_LINK_RE = re.compile(
    r'https://(chatgpt\.com/share|gemini\.google\.com/share|claude\.ai/share|chat\.deepseek\.com/share|www\.kimi\.com/share)',
    re.IGNORECASE
)

# Code block detection
CODE_BLOCK_RE = re.compile(r'^```\w+', re.MULTILINE)
INLINE_CODE_RE = re.compile(r'`[^`]+`')

# Section detection
SECTION_RE = re.compile(r'^#{2,3}\s+(.+)$', re.MULTILINE)

# Table detection
TABLE_RE = re.compile(r'^\|.*\|', re.MULTILINE)

# Date patterns for dev notes / debug logs
# YYYY-MM-DD, YYYY/MM/DD, or Month DD, YYYY
DATE_RE = re.compile(
    r'(?:^#{1,4}\s*(\d{4}[-/]\d{2}[-/]\d{2})'  # heading with date
    r'|\b(\d{4}[-/]\d{2}[-/]\d{2})\s+(session notes|parity recap|progress update|debugging|testing|check|meeting|review|log)'
    r'|\b(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)[a-z]*\s+\d{1,2},?\s+\d{4}\b'
    r'|\b\d{4}[-/]\d{2}[-/]\d{2}\b.*(?:notes|recap|update|status|log|debug|fix|parity|test))',
    re.MULTILINE | re.IGNORECASE
)

# ─────────────────────────────────────────────
# KEYWORD DICTIONARIES (type → keywords)
# ─────────────────────────────────────────────

USER_GUIDE_KEYWORDS = [
    "user guide", "manual", "how to use", "usage guide", "getting started",
    "prerequisites", "installation", "setup", "configuration",
    "overview", "introduction", "what is", "purpose",
    "tutorial", "quick start", "step-by-step", "walkthrough", "beginner",
    "lesson", "course", "example usage", "try it", "run the following",
    "copy and paste", " hands-on",
]

API_DOC_KEYWORDS = [
    "api", "interface", "function", "method", "class", "parameters",
    "returns", "arguments", "signature", "call", "invoke",
    "public api", "module", "module reference", "api reference",
    "docstring", "documented", "exposed", "wrapper",
]

DESIGN_DOC_KEYWORDS = [
    "design document", "architecture", "system design", "data flow",
    "component", "subsystem", "module architecture", "design pattern",
    "algorithmic strategy", "implementation plan", "blueprint",
    "global memory map", "kernel launch", "memory layout",
]

IMPLEMENTATION_REPORT_KEYWORDS = [
    "implementation report", "what we achieved", "validation results",
    "key takeaways", "test outcome", "results", "benchmark",
    "performance", "speedup", "accuracy", "comparison",
    "reproducible", "verified", "validated", "working",
]

TECHNICAL_NOTE_KEYWORDS = [
    "mathematical", "derivation", "equation", "formula", "proof",
    "algorithm", "complexity", "convergence", "stability",
    "numerical", "analytical", "theoretical", "physics",
    "debugging", "nan", "precision",
]

PROGRESS_REPORT_KEYWORDS = [
    "session notes", "parity recap", "progress update",
    "root cause", "problems encountered", "fixes applied", "next steps",
    "current status", "test outcome", "key findings", "bottlenecks",
    "crashes", "workaround", "regression", "memory leak", "buffer overrun",
    "race condition", "segfault", "mismatch", "discrepancy",
    "double-free", "out-of-bounds", "undefined behavior", "hot path",
    "TODO", "FIXME", "HACK", "WIP", "WARN",
]

DRAFT_NOTE_KEYWORDS = [
    "draft", "sketch", "proposal", "idea", "exploratory",
    "work in progress", "unfinished", "incomplete", "placeholder",
    "future work", "to be done", "not yet", "preliminary",
]

LLM_KEYWORDS = [
    "# USER\n", "# Gemini\n", "# ChatGPT", "# Claude", "# DeepSeek",
    "# Kimi", "# SWE-", "# GPT-", "# Copilot", "# Devin",
    "# Post-mortem", "# gpt-", "# o1", "# o3",
]

# Filename patterns (type → patterns)
FILENAME_PATTERNS = {
    "user_guide": ["readme", "guide", "manual", "overview", "introduction", "tutorial", "_tutorial", "quickstart", "getting_started"],
    "api_doc": ["_api", "api_", "_reference", "reference_", "_doc", "_docs", "_manifest"],
    "design_doc": ["_design", "design_", "_architecture", "architecture_"],
    "implementation_report": ["_report", "report_", "_summary", "summary_", "_progress", "progress_"],
    "technical_note": ["_note", "_notes", "_derivation", "_math", "_algorithm"],
    "progress_report": ["_debug", "debug_", "_dev", "dev_", "_log", "_session", "_progress", "progress_", "_status", "status_"],
    "draft_note": ["_draft", "draft_", "_wip", "_todo", "_plan", "_sketch", "_research"],
}


def count_keywords(text_lower: str, keywords: list) -> int:
    """Count occurrences of keywords in text (case-insensitive)."""
    return sum(1 for kw in keywords if kw.lower() in text_lower)


def score_content(text: str, filename: str) -> Dict[str, int]:
    """
    Score a markdown file across all document-type categories.
    Returns dict: type → score (higher = stronger evidence).
    """
    lower = text.lower()
    fname_lower = filename.lower()
    basename = os.path.basename(fname_lower)

    scores = {
        "llm_chat": 0,
        "user_guide": 0,
        "api_doc": 0,
        "design_doc": 0,
        "implementation_report": 0,
        "technical_note": 0,
        "progress_report": 0,
        "draft_note": 0,
    }

    # ── LLM Chat (highest priority detections) ──
    if SHARE_LINK_RE.search(text):
        scores["llm_chat"] += 15

    llm_headings = LLM_HEADING_RE.findall(text)
    if llm_headings:
        scores["llm_chat"] += len(llm_headings) * 6

    user_turns = len(re.findall(r'^#{1,3}\s*USER\s*$', text, re.MULTILINE | re.IGNORECASE))
    if user_turns > 0:
        scores["llm_chat"] += user_turns * 4

    # LLM conversational patterns
    if re.search(r'^#{1,3}\s*(Gemini|ChatGPT|Claude|DeepSeek)\s*$', text, re.MULTILINE | re.IGNORECASE):
        scores["llm_chat"] += 8

    # ── Content keyword scoring ──
    scores["user_guide"] = count_keywords(lower, USER_GUIDE_KEYWORDS)
    scores["api_doc"] = count_keywords(lower, API_DOC_KEYWORDS)
    scores["design_doc"] = count_keywords(lower, DESIGN_DOC_KEYWORDS)
    scores["implementation_report"] = count_keywords(lower, IMPLEMENTATION_REPORT_KEYWORDS)
    scores["technical_note"] = count_keywords(lower, TECHNICAL_NOTE_KEYWORDS)
    scores["progress_report"] = count_keywords(lower, PROGRESS_REPORT_KEYWORDS)
    scores["draft_note"] = count_keywords(lower, DRAFT_NOTE_KEYWORDS)

    # ── Date-based progress_report boost ──
    date_matches = len(DATE_RE.findall(text))
    if date_matches > 0:
        scores["progress_report"] += date_matches * 3  # stronger date signal

    # ── Structural signals ──
    section_count = len(SECTION_RE.findall(text))
    code_blocks = len(CODE_BLOCK_RE.findall(text))
    tables = len(TABLE_RE.findall(text))

    # User guides / tutorials often have many code blocks
    scores["user_guide"] += min(code_blocks, 5)
    # API docs often have tables
    scores["api_doc"] += min(tables // 2, 3)
    # Reports often have tables
    scores["implementation_report"] += min(tables // 2, 3)

    # ── Filename pattern scoring ──
    for doc_type, patterns in FILENAME_PATTERNS.items():
        for pat in patterns:
            if pat in basename:
                scores[doc_type] += 8
                break

    # Special filename overrides
    if basename == "readme.md" or basename == "index.md":
        scores["user_guide"] += 15
    if "_tutorial" in basename or basename.startswith("tutorial_"):
        scores["user_guide"] += 12
    if "_api" in basename or basename.startswith("api_"):
        scores["api_doc"] += 12
    if "_report" in basename or basename.startswith("report_"):
        scores["implementation_report"] += 10
    if "_design" in basename or basename.startswith("design_"):
        scores["design_doc"] += 10
    if "_debug" in basename or basename.startswith("debug_"):
        scores["progress_report"] += 10
    if "_progress" in basename or basename.startswith("progress_"):
        scores["progress_report"] += 12
    if "_draft" in basename or basename.startswith("draft_"):
        scores["draft_note"] += 12

    # ── Path-based signals ──
    if "/topics/" in fname_lower:
        # Topics are usually reports or design docs
        scores["implementation_report"] += 2
        scores["design_doc"] += 1
    if "/markdown/" in fname_lower or "/doc/" in fname_lower:
        scores["user_guide"] += 1
        scores["api_doc"] += 1
    if "/devnotes/" in fname_lower or "/dev_notes/" in fname_lower:
        scores["progress_report"] += 4
        scores["draft_note"] += 2
    if "/tests/" in fname_lower and ("readme" in basename or "notes" in basename):
        scores["user_guide"] += 3
    if "/tutorials/" in fname_lower or "/examples/" in fname_lower:
        scores["user_guide"] += 5

    # ── Cross-type penalties ──
    if scores["llm_chat"] > 5:
        # LLM chats shouldn't be classified as other types
        for k in scores:
            if k != "llm_chat":
                scores[k] = max(0, scores[k] - scores["llm_chat"] // 2)

    return scores


def is_config(filepath: str) -> bool:
    """Agent/IDE config files (skills, workflows, AGENTS, coding rules, protocols, etc.)"""
    f = filepath.lower()
    if "/.windsurf/" in f or "/.cursor/" in f or "/.devin/" in f or "/.aider/" in f:
        return True
    if "/skills/" in f and f.endswith("skill.md"):
        return True
    if "/workflows/" in f and (f.endswith(".md") or f.endswith("_workflow")):
        return True
    if "/protocols/" in f and f.endswith(".md"):
        return True
    if "/codingrules/" in f and f.endswith(".md"):
        return True
    basename = os.path.basename(f)
    if basename in ("agents.md", "agents_alternative.md", "global_rules.md", "workflow_template.md"):
        return True
    return False


# ShaderToy hash suffix pattern: filenames like Noise_gradient_2D_XdXGW8.md
SHADER_HASH_RE = re.compile(r'_[A-Za-z0-9]{6}\.md$', re.IGNORECASE)


def is_shader_snippet(filepath: str) -> bool:
    """Detect ShaderToy downloaded shader files (GLSL code in markdown wrapper)."""
    f = filepath.lower()
    if "shadertoy_inspiration" in f:
        return True
    # Hash suffix pattern in filename (e.g. _XdXGW8.md, _MslfDf.md)
    if SHADER_HASH_RE.search(os.path.basename(f)):
        return True
    return False


WIP_MARKERS = [
    "TODO", "FIXME", "HACK", "WIP", "DRAFT",
    "TEMPORARY", "PLACEHOLDER", "PENDING",
    "UNFINISHED", "INCOMPLETE", "ROUGH", "SKETCH",
    "FUTURE WORK", "TO BE DONE",
]


def assess_quality(filepath: str, text: str, doc_type: str) -> Tuple[str, int]:
    """
    Assess document quality (high/medium/low) and return confidence 0-100.
    Returns: (quality, confidence)
    """
    lower = text.lower()
    fname = os.path.basename(filepath).lower()

    section_count = len(SECTION_RE.findall(text))
    text_len = len(text)
    marker_count = sum(1 for m in WIP_MARKERS if m.lower() in lower)
    code_blocks = len(CODE_BLOCK_RE.findall(text))
    tables = len(TABLE_RE.findall(text))

    # ── Base quality scoring ──
    quality_score = 0

    # Structure bonus
    if section_count >= QUALITY_HIGH_MIN_SECTIONS:
        quality_score += 3
    elif section_count >= 4:
        quality_score += 1
    elif section_count <= QUALITY_LOW_MAX_SECTIONS:
        quality_score -= 2

    # Length bonus
    if text_len >= QUALITY_HIGH_MIN_LENGTH:
        quality_score += 2
    elif text_len >= 3000:
        quality_score += 1
    elif text_len <= QUALITY_LOW_MAX_LENGTH:
        quality_score -= 1

    # Content richness
    quality_score += min(code_blocks, 2)
    quality_score += min(tables, 2)

    # WIP penalty
    if marker_count >= 5:
        quality_score -= 3
    elif marker_count >= 3:
        quality_score -= 2
    elif marker_count >= 1:
        quality_score -= 1

    # Filename quality signals
    if fname.startswith("readme") or fname.startswith("index"):
        quality_score += 2
    if "_draft" in fname or "_wip" in fname or "_todo" in fname:
        quality_score -= 2
    if "_report" in fname or "_guide" in fname or "_tutorial" in fname:
        quality_score += 1

    # LLM chat is inherently lower quality as documentation
    if doc_type == "llm_chat":
        quality_score -= 2

    # Type-specific adjustments
    if doc_type == "api_doc" and code_blocks < 2:
        quality_score -= 1  # API docs should have examples
    if doc_type == "user_guide" and code_blocks < 2:
        quality_score -= 1  # Guides should have examples
    if doc_type == "implementation_report" and tables < 1:
        quality_score -= 1  # Reports should have data tables

    # Map to quality level
    if quality_score >= 4:
        quality = "high"
    elif quality_score >= 1:
        quality = "medium"
    else:
        quality = "low"

    # Confidence based on strength of classification signals
    confidence = min(100, max(20, 50 + quality_score * 10))

    return quality, confidence


# Generic meta filenames that are not user guides per se
META_NAMES = {"readme.md", "codemap.md", "agents.md", "agents_alternative.md",
              "dox_agents.md", "index.md", "code_of_conduct.md", "contributing.md",
              "license.md", "changelog.md", "history.md", "todo.md", "todo_effi.md"}


def is_meta_doc(filepath: str) -> bool:
    """Detect generic project-meta files (README, CODEMAP, AGENTS, index, etc.)."""
    fname = os.path.basename(filepath).lower()
    if fname in META_NAMES:
        return True
    # CODEMAP copy etc.
    if fname.startswith("codemap") and fname.endswith(".md"):
        return True
    # AGENTS in any directory
    if fname.startswith("agents") and fname.endswith(".md"):
        return True
    # DOX_AGENTS variants
    if fname.startswith("dox_agents") and fname.endswith(".md"):
        return True
    return False


def classify_file(filepath: str, text: str) -> Tuple[str, int]:
    """
    Classify a markdown file and return (doc_type, confidence_score).
    """
    # Highest priority: generic meta files (README, CODEMAP, AGENTS, index)
    if is_meta_doc(filepath):
        return "meta_doc", 95

    # Second priority: config files (skills, workflows, protocols, coding rules)
    if is_config(filepath):
        return "config", 95

    # Third priority: shader snippets (ShaderToy downloads — GLSL code, not docs)
    if is_shader_snippet(filepath):
        return "shader_snippet", 95

    # Design document detection: first line contains "Design Document"
    first_lines = text[:500].lower()
    if "design document" in first_lines or "# design document" in first_lines:
        return "design_doc", 85

    # Principles/guide detection: filename contains "principles" or "guidelines"
    basename = os.path.basename(filepath).lower()
    if "principles" in basename or "guidelines" in basename:
        return "technical_note", 85

    scores = score_content(text, filepath)

    # LLM chat: require share link or explicit LLM headings, not just keyword matches
    if SHARE_LINK_RE.search(text):
        return "llm_chat", min(100, 50 + scores["llm_chat"] * 3)
    llm_headings = LLM_HEADING_RE.findall(text)
    if llm_headings and scores["llm_chat"] >= 15:
        return "llm_chat", min(100, 50 + scores["llm_chat"] * 3)

    # Find best type
    best_type = max((k for k in scores if k != "llm_chat"), key=scores.get)
    best_score = scores[best_type]

    # If no strong signals, default based on content
    if best_score <= 0:
        lower = text.lower()
        if re.search(r'\$\$.+\$\$', text) or re.search(r'\\[a-zA-Z]+\{', text):
            return "technical_note", 60
        return "user_guide", 50

    # If close between two types, pick the more specific one
    sorted_types = sorted(
        ((k, v) for k, v in scores.items() if k != "llm_chat"),
        key=lambda x: x[1], reverse=True
    )
    if len(sorted_types) >= 2:
        top_type, top_score = sorted_types[0]
        second_type, second_score = sorted_types[1]
        if top_score - second_score <= 2:
            # Ambiguous — lower confidence
            confidence = 40 + top_score * 5
        else:
            confidence = 50 + top_score * 5
    else:
        confidence = 50 + best_score * 5

    confidence = min(100, max(20, confidence))

    return best_type, confidence


def get_git_tracked_files(repo_path: str) -> Set[str]:
    """Return set of absolute paths of files tracked by git."""
    try:
        result = subprocess.run(
            ["git", "-C", repo_path, "ls-files"],
            capture_output=True, text=True, check=True
        )
        tracked = set()
        for line in result.stdout.strip().split("\n"):
            if line:
                tracked.add(os.path.abspath(os.path.join(repo_path, line)))
        return tracked
    except Exception as e:
        print(f"WARN: cannot get git tracked files: {e}", file=sys.stderr)
        return set()


# Category → file extension mapping (applied to git-tracked files only)
# Categories NOT listed here are never renamed (meta_doc, config, api_doc, design_doc)
CATEGORY_EXT = {
    "user_guide": "guide",
    "implementation_report": "report",
    "technical_note": "note",
    "progress_report": "progress",
    "draft_note": "plan",
    "llm_chat": "chat",
}


# Extra suffixes to strip when a category switches its canonical extension
STRIP_ALIASES = {
    "chat": ["discussion"],
    "discussion": ["chat"],
}


def clean_stem_for_ext(stem: str, ext: str) -> str:
    """
    Remove trailing words matching `ext` (and aliases) from stem to avoid duplication.
    e.g.  Foo_discussion  + chat  →  Foo
          Foo.report  + report  →  Foo
    """
    stem_lower = stem.lower()
    exts_to_strip = [ext.lower()] + [a.lower() for a in STRIP_ALIASES.get(ext, [])]

    for e in exts_to_strip:
        # Remove trailing .ext
        if stem_lower.endswith(f".{e}"):
            stem = stem[:-(len(e) + 1)]
            stem_lower = stem.lower()

        # Remove trailing _ext or -ext
        for sep in ("_", "-"):
            if stem_lower.endswith(f"{sep}{e}"):
                stem = stem[:-(len(e) + len(sep))]
                stem_lower = stem.lower()
                break

    return stem


def backup_files(
    results: List[Tuple[str, str, str, int, Dict]],
    repo_path: str,
    backup_dir: str = ".backup_md"
) -> int:
    """
    Copy all tracked .md files to a backup directory, preserving relative structure.
    Returns number of files backed up.
    """
    tracked = get_git_tracked_files(repo_path)
    backup_root = os.path.join(repo_path, backup_dir)
    count = 0
    for old_path, doc_type, *_ in results:
        abs_path = os.path.abspath(old_path)
        if abs_path not in tracked:
            continue
        rel = os.path.relpath(old_path, repo_path)
        dest = os.path.join(backup_root, rel)
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        shutil.copy2(abs_path, dest)
        count += 1
    print(f"Backed up {count} tracked .md files to {backup_root}")
    return count


def rename_by_category(
    results: List[Tuple[str, str, str, int, Dict]],
    repo_path: str,
    tracked: Set[str],
    dry_run: bool = True
) -> List[Tuple[str, str]]:
    """
    Rename git-tracked files according to CATEGORY_EXT mapping.
    Returns list of (old_path, new_path) tuples.
    Skips excluded, non-tracked, and already-correctly-named files.
    Idempotent: running twice produces no additional changes.
    """
    renamed = []

    for old_path, doc_type, *_ in results:
        if doc_type not in CATEGORY_EXT:
            continue

        abs_path = os.path.abspath(old_path)
        if abs_path not in tracked:
            print(f"SKIP (not tracked by git): {old_path}")
            continue

        if is_excluded_from_rename(old_path, repo_path):
            print(f"SKIP (excluded path/pattern): {old_path}")
            continue

        ext = CATEGORY_EXT[doc_type]
        dirname = os.path.dirname(old_path)
        basename = os.path.basename(old_path)
        stem, _ = os.path.splitext(basename)  # strip .md

        # Guard: already has correct double extension? (e.g. Foo.chat.md)
        if stem.lower().endswith(f".{ext}"):
            print(f"SKIP (already has .{ext}.md): {old_path}")
            continue

        # Clean stem: remove existing trailing ext words
        new_stem = clean_stem_for_ext(stem, ext)

        new_basename = f"{new_stem}.{ext}.md"
        new_path = os.path.join(dirname, new_basename)

        # No-op / idempotency check
        if new_path == old_path:
            print(f"SKIP (no change): {old_path}")
            continue

        # Collision check
        if os.path.exists(new_path):
            print(f"SKIP (collision): {old_path} -> {new_path}")
            continue

        if dry_run:
            print(f"DRY-RUN: git mv {old_path} {new_path}")
            renamed.append((old_path, new_path))
        else:
            try:
                subprocess.run(
                    ["git", "-C", repo_path, "mv", old_path, new_path],
                    check=True, capture_output=True
                )
                print(f"RENAMED: {old_path} -> {new_path}")
                renamed.append((old_path, new_path))
            except subprocess.CalledProcessError as e:
                print(f"FAIL: {old_path} -> {new_path}: {e.stderr.decode()}", file=sys.stderr)

    return renamed


def main():
    parser = argparse.ArgumentParser(description="Classify markdown files in repo")
    parser.add_argument("--rename", action="store_true", help="Rename git-tracked files to category extensions (.guide.md, .report.md, etc.)")
    parser.add_argument("--rename-dry-run", action="store_true", help="Show what would be renamed without doing it")
    parser.add_argument("--backup", action="store_true", help="Copy all tracked .md files to backup dir before any rename")
    parser.add_argument("--backup-dir", default=".backup_md", help="Directory for markdown backups (default: .backup_md)")
    parser.add_argument("--commit", action="store_true", help="After renaming, auto-commit the changes to git")
    parser.add_argument("--file-list", default="/tmp/all_md_files.txt", help="Path to file containing list of .md files")
    args = parser.parse_args()

    # Read the list of all .md files
    with open(args.file_list) as f:
        paths = [p.strip() for p in f if p.strip()]

    # Results: list of (path, type, quality, confidence, scores)
    results = []

    for p in paths:
        try:
            with open(p, "r", encoding="utf-8", errors="ignore") as f:
                text = f.read(15_000)
        except Exception as e:
            print(f"WARN: cannot read {p}: {e}", file=sys.stderr)
            continue

        # 1. Meta-doc filter (highest priority: README, CODEMAP, AGENTS, index)
        if is_meta_doc(p):
            doc_type = "meta_doc"
            confidence = 95
            scores = {}
        # 2. Config filter (SKILL.md, workflows, protocols, coding rules, etc.)
        elif is_config(p):
            doc_type = "config"
            confidence = 95
            scores = {}
        # 3. Shader snippet filter (ShaderToy downloads — GLSL code, not docs)
        elif is_shader_snippet(p):
            doc_type = "shader_snippet"
            confidence = 95
            scores = {}
        else:
            doc_type, confidence = classify_file(p, text)
            scores = score_content(text, p)

        # 2. Quality assessment
        quality, _ = assess_quality(p, text, doc_type)

        results.append((p, doc_type, quality, confidence, scores))

    # Sort by type, then quality, then path
    TYPE_ORDER = [
        "meta_doc", "config", "shader_snippet", "user_guide", "api_doc", "design_doc",
        "implementation_report", "technical_note", "progress_report", "draft_note", "llm_chat",
    ]
    QUALITY_ORDER = {"high": 0, "medium": 1, "low": 2}

    def sort_key(item):
        p, doc_type, quality, conf, scores = item
        type_idx = TYPE_ORDER.index(doc_type) if doc_type in TYPE_ORDER else 99
        qual_idx = QUALITY_ORDER.get(quality, 1)
        return (type_idx, qual_idx, p)

    results.sort(key=sort_key)

    # Group by type
    categories: Dict[str, List[Tuple]] = {}
    for p, doc_type, quality, confidence, scores in results:
        categories.setdefault(doc_type, []).append((p, quality, confidence, scores))

    # Print summary to stdout
    print("# Markdown File Classification (v2)")
    print("")
    total = len(results)
    print(f"Total files: {total}")
    print("")

    for cat in TYPE_ORDER:
        if cat not in categories:
            continue
        items = categories[cat]
        print(f"## {cat} ({len(items)})")
        for p, quality, confidence, scores in items:
            rel = os.path.relpath(p, REPO)
            conf_str = f"{confidence}%"
            print(f"- [{quality}] {rel}  (conf: {conf_str})")
        print("")

    # ── Backup if requested ──
    if args.backup:
        print("\n" + "=" * 50)
        print("Backing up tracked .md files")
        print("=" * 50)
        backup_files(results, REPO, backup_dir=args.backup_dir)

    # ── Rename files by category if requested ──
    if args.rename or args.rename_dry_run:
        dry_run = args.rename_dry_run or not args.rename
        print("\n" + "=" * 50)
        print(f"{'DRY-RUN' if dry_run else 'LIVE'}: Renaming files by category")
        print("=" * 50)

        tracked = get_git_tracked_files(REPO)
        renamed = rename_by_category(results, REPO, tracked, dry_run=dry_run)

        if dry_run:
            print(f"\nWould rename {len(renamed)} files. Use --rename to execute.")
        else:
            print(f"\nRenamed {len(renamed)} files.")
            if renamed and args.commit:
                try:
                    subprocess.run(
                        ["git", "-C", REPO, "commit", "-m",
                         f"Rename {len(renamed)} markdown files to category extensions"],
                        check=True, capture_output=True
                    )
                    print(f"Committed rename of {len(renamed)} files.")
                except subprocess.CalledProcessError as e:
                    print(f"FAIL: git commit failed: {e.stderr.decode()}", file=sys.stderr)

    # Write structured markdown output
    out_path = os.path.join(REPO, "doc", "MarkdownFileClassification.md")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("# Markdown File Classification\n\n")
        f.write("This document lists all `.md` files in the repository and classifies them by type, quality, and confidence.\n")
        f.write("Generated automatically by `scripts/classify_md.py` (hybrid heuristic approach).\n\n")
        f.write("## Legend\n\n")
        f.write("- **Type**: document classification based on content and purpose\n")
        f.write("- **Quality**: `high` (polished), `medium` (useful but informal), `low` (temporary/misclassified)\n")
        f.write("- **Confidence**: heuristic certainty score (0-100%, higher = more certain)\n\n")
        f.write("## Categories\n\n")
        f.write("- **meta_doc** — Project meta files: README, CODEMAP, AGENTS, index (AI/agentic navigation)\n")
        f.write("- **config** — Agent/IDE configuration: SKILL.md, workflows, protocols, coding rules, .windsurf, .cursor, .devin\n")
        f.write("- **shader_snippet** — Downloaded ShaderToy GLSL shaders (not documentation)\n")
        f.write("- **user_guide** — Manuals, how-to guides, and tutorials\n")
        f.write("- **api_doc** — API reference documentation (functions, classes, parameters)\n")
        f.write("- **design_doc** — Architecture and design documents\n")
        f.write("- **implementation_report** — Polished reports with validation/results\n")
        f.write("- **technical_note** — Mathematical derivations, algorithm explanations\n")
        f.write("- **progress_report** — Debug logs, session notes, dated progress entries\n")
        f.write("- **draft_note** — Temporary / exploratory / WIP / unfinished docs\n")
        f.write("- **llm_chat** — Conversation transcripts with LLMs\n\n")

        for cat in TYPE_ORDER:
            if cat not in categories:
                continue
            items = categories[cat]
            f.write(f"## {cat} ({len(items)})\n\n")
            for p, quality, confidence, scores in items:
                rel = os.path.relpath(p, REPO)
                f.write(f"- `[{quality}]` `{rel}`  (conf: {confidence}%)\n")
            f.write("\n")

        # Add quality summary table
        f.write("## Quality Summary\n\n")
        quality_counts: Dict[str, int] = {}
        for _, _, quality, _, _ in results:
            quality_counts[quality] = quality_counts.get(quality, 0) + 1
        f.write("| Quality | Count | Description |\n")
        f.write("|---------|-------|-------------|\n")
        for q in ["high", "medium", "low"]:
            if q in quality_counts:
                desc = {"high": "Polished, comprehensive, well-structured",
                        "medium": "Useful but informal or incomplete",
                        "low": "Temporary, misclassified, or very rough"}[q]
                f.write(f"| {q} | {quality_counts[q]} | {desc} |\n")
        f.write("\n")

    print(f"\nWrote classification to: {out_path}")


if __name__ == "__main__":
    main()
