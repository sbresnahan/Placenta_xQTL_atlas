#!/usr/bin/env python3
"""
config_get.py — Read config.yml and emit shell export statements.

Called once per shell script via:
    eval "$(python3 "${SCRIPTS_DIR}/config_get.py" "${CONFIG}" --cohort "${COHORT}")"

Exports all config values as UPPERCASE shell variables:
  - Top-level scalar keys -> UPPERCASE var (e.g. output_base -> OUTPUT_BASE)
  - Nested dict keys flattened with underscores
    (e.g. rna_editing.edit_sites_min_coverage -> RNA_EDITING_EDIT_SITES_MIN_COVERAGE)
  - If --cohort given: cohort-specific keys exported as
    STRANDEDNESS, FASTQ_DIR, FASTQ_MAP, SAMPLES_FILE

Booleans are emitted as lowercase 'true'/'false' so shell scripts can test
them with [ "$VAR" = "true" ].

Does NOT require PyYAML — uses a built-in line-based parser that handles
arbitrary nesting depth via indentation.

Usage:
    python3 config_get.py <config.yml> [--cohort <name>]
"""

import sys


def parse_config(text):
    """Parse config.yml into a nested dict.

    Handles arbitrary nesting depth via indentation. Each increase in
    indentation (by any number of spaces) creates a nested dict.

    Does NOT handle lists, multi-line strings, or flow-style YAML.
    """
    root = {}
    # Stack of (indent, dict) pairs — tracks the current nesting path.
    stack = [(-1, root)]

    for raw_line in text.splitlines():
        # Strip inline comments.
        stripped = raw_line.split("#", 1)[0].rstrip()
        if not stripped.strip():
            continue

        indent = len(stripped) - len(stripped.lstrip(" "))
        line = stripped.strip()

        if ":" not in line:
            continue

        key, _, value = line.partition(":")
        key = key.strip()
        value = value.strip()

        # Pop the stack until we find the parent (indent < current indent).
        while stack and stack[-1][0] >= indent:
            stack.pop()

        parent = stack[-1][1] if stack else root

        if value == "":
            # Nested block — push onto stack so subsequent lines go inside.
            new_dict = {}
            parent[key] = new_dict
            stack.append((indent, new_dict))
        else:
            parent[key] = coerce_scalar(value)

    return root


def coerce_scalar(value):
    """Coerce a YAML scalar string to int/float/bool/str."""
    if value in ("true", "True", "TRUE"):
        return True
    if value in ("false", "False", "FALSE"):
        return False
    if value in ("null", "Null", "~", ""):
        return None
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        pass
    if (value.startswith('"') and value.endswith('"')) or \
       (value.startswith("'") and value.endswith("'")):
        return value[1:-1]
    return value


def shell_quote(value):
    """Single-quote a value for safe shell export.

    Booleans are emitted as lowercase 'true'/'false' so shell scripts can
    test them with [ "$VAR" = "true" ].
    """
    if value is None:
        return "''"
    if isinstance(value, bool):
        s = "true" if value else "false"
    else:
        s = str(value)
    return "'" + s.replace("'", "'\"'\"'") + "'"


def flatten(prefix, obj, out):
    """Flatten nested dict into UPPERCASE_KEY -> value pairs."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            key = (prefix + "_" + k).upper() if prefix else k.upper()
            flatten(key, v, out)
    else:
        out[prefix] = obj


def main():
    if len(sys.argv) < 2:
        sys.stderr.write("Usage: config_get.py <config.yml> [--cohort <name>]\n")
        sys.exit(1)

    config_path = sys.argv[1]
    cohort = None
    if "--cohort" in sys.argv:
        idx = sys.argv.index("--cohort")
        if idx + 1 < len(sys.argv):
            cohort = sys.argv[idx + 1]

    try:
        with open(config_path) as f:
            config = parse_config(f.read())
    except Exception as e:
        sys.stderr.write("ERROR loading config %s: %s\n" % (config_path, e))
        sys.exit(1)

    if not isinstance(config, dict):
        sys.stderr.write("ERROR: config is not a dict\n")
        sys.exit(1)

    # Flatten all top-level keys (including nested dicts like rna_editing,
    # intron_retention, lsf) into UPPERCASE vars.
    flat = {}
    for key, value in config.items():
        if key == "cohorts":
            continue  # handled separately below
        flatten("", {key: value}, flat)

    # Emit exports for all flattened keys.
    for key in sorted(flat.keys()):
        value = flat[key]
        if value is None:
            continue
        print("export %s=%s" % (key, shell_quote(value)))

    # Cohort-specific keys.
    if cohort:
        cohorts = config.get("cohorts", {})
        if cohort not in cohorts:
            sys.stderr.write(
                "ERROR: cohort '%s' not found in config. Available: %s\n"
                % (cohort, list(cohorts.keys()))
            )
            sys.exit(1)
        cohort_cfg = cohorts[cohort]
        if not isinstance(cohort_cfg, dict):
            sys.stderr.write(
                "ERROR: cohort '%s' is not a dict in config.\n" % cohort
            )
            sys.exit(1)
        for key, value in cohort_cfg.items():
            if value is None:
                continue
            print("export %s=%s" % (key.upper(), shell_quote(value)))


if __name__ == "__main__":
    main()
