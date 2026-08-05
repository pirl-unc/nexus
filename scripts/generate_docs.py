#!/usr/bin/env python3
"""
Auto-generate Quarto documentation site from existing README.md and .nf files.
Outputs to docs/ directory.
"""


import ast
import os
import re
import shutil
import yaml


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SUBWORKFLOWS_DIR = os.path.join(REPO_ROOT, "src", "nexuslib", "subworkflows")
WORKFLOWS_DIR = os.path.join(REPO_ROOT, "src", "nexuslib", "workflows")
UTILITIES_DIR = os.path.join(REPO_ROOT, "src", "nexuslib", "scripts")
PYPROJECT_PATH = os.path.join(REPO_ROOT, "pyproject.toml")
DOCS_DIR = os.path.join(REPO_ROOT, "docs")

CATEGORY_LABELS = {
    "alignment": "Alignment",
    "antigen_prediction": "Antigen Prediction",
    "assembly": "Assembly",
    "haplotagging": "Haplotagging",
    "hla_typing": "HLA Typing",
    "isoform_characterization": "Isoform Characterization",
    "peptide_prediction": "Peptide Prediction",
    "quantification": "Quantification",
    "read_clustering": "Read Clustering",
    "read_error_correction": "Read Error Correction",
    "sequencing_simulation": "Sequencing Simulation",
    "utilities": "Utilities",
    "variant_annotation": "Variant Annotation",
    "variant_calling": "Variant Calling",
}

WORKFLOW_LABELS = {
    "assembly_long-read-rna": "Assembly (Long-read RNA)",
    "assembly_short-read-rna": "Assembly (Short-read RNA)",
    "haplotagging_long-read-dna": "Haplotagging (Long-read DNA)",
    "haplotagging_long-read-rna": "Haplotagging (Long-read RNA)",
    "haplotagging_short-read-dna": "Haplotagging (Short-read DNA)",
    "hla_typing_long-read": "HLA typing (Long-read)",
    "hla_typing_long-read-dna": "HLA typing (Long-read DNA)",
    "hla_typing_long-read-rna": "HLA typing (Long-read RNA)",
    "hla_typing_short-read": "HLA typing (Short-read)",
    "hla_typing_short-read-dna": "HLA typing (Short-read DNA)",
    "hla_typing_short-read-rna": "HLA typing (Short-read RNA)",
    "isoform_characterization_long-read": "Isoform Characterization (Long-read)",
    "isoform_characterization_short-read": "Isoform Characterization (Short-read)",
    "quantification_long-read-rna": "Quantification (Long-read RNA)",
    "quantification_short-read-rna": "Quantification (Short-read RNA)",
    "variant_calling_long-read-dna-germline": "Variant Calling (Long-read DNA Germline)",
    "variant_calling_long-read-dna-somatic": "Variant Calling (Long-read DNA Somatic)",
    "variant_calling_long-read-rna": "Variant Calling (Long-read RNA)",
    "variant_calling_short-read-dna-germline": "Variant Calling (Short-read DNA Germline)",
    "variant_calling_short-read-dna-somatic": "Variant Calling (Short-read DNA Somatic)",
}


def is_source_dir(path, dirname):
    """Check whether a directory holds a subworkflow/workflow category or tool.

    Skips dotfiles and dunder directories so Python bytecode caches
    (`__pycache__`, which any test run leaves behind under `src/`) are
    never mistaken for a category or a tool. Without this guard the
    generator emits a bogus `__pycache__` sidebar section plus a link to
    a page it never writes, which breaks `quarto render`.
    """
    return os.path.isdir(path) and not dirname.startswith(('.', '__'))


def extract_tool_name(dirname):
    """Extract the tool name from a subworkflow directory name.

    Strips the canonical `<category>_` prefix if present. If no known
    category matches (e.g. the directory name has a typo like
    `sequenicng_simulation_beers2`), falls back to returning the last
    underscore-separated token — this is consistent with the convention
    that subworkflow directories are named `<category>_<tool>` and tool
    names use dashes (never underscores), so the trailing token is
    always the tool name.
    """
    for category in CATEGORY_LABELS:
        prefix = category + "_"
        if dirname.startswith(prefix):
            return dirname[len(prefix):]
    if "_" in dirname:
        return dirname.rsplit("_", 1)[-1]
    return dirname


def resolve_nf_path(tool_path, tool_dir):
    """Resolve the .nf file inside a subworkflow/workflow directory.

    Tries the canonical `<tool_dir>.nf` first. If that file doesn't
    exist (e.g. because the directory name has a typo but the .nf
    inside is correctly spelled — as is the case for
    `sequenicng_simulation_beers2/sequencing_simulation_beers2.nf`),
    falls back to the sole .nf file in the directory.

    Returns the tuple (nf_filename, nf_path). When no .nf is found,
    returns the canonical (non-existent) names so the caller's
    `read_nf_file` path still degrades gracefully to empty content.
    """
    canonical = f"{tool_dir}.nf"
    canonical_path = os.path.join(tool_path, canonical)
    if os.path.exists(canonical_path):
        return canonical, canonical_path

    if not os.path.isdir(tool_path):
        return canonical, canonical_path

    nf_files = sorted(f for f in os.listdir(tool_path) if f.endswith('.nf'))
    if len(nf_files) == 1:
        fname = nf_files[0]
        return fname, os.path.join(tool_path, fname)

    # Zero matches or ambiguous (multiple .nf files in one dir): keep
    # the canonical path so the caller's existing missing-file behavior
    # applies. Multiple .nf in a single subworkflow dir is unusual
    # enough that picking arbitrarily would silently mislabel pages.
    return canonical, canonical_path


def read_extras_block(tool_path):
    """Return hand-written extra docs for a tool, or '' if none exist.

    `generate_docs.py` overwrites the per-tool .qmd on every run, so any
    long-form content (deep parameter references, gotchas, examples)
    written directly into the generated .qmd is lost on regen. To keep
    such content version-controlled and regen-safe, drop a Markdown file
    named `docs_extras.md` next to the subworkflow's .nf file:

        src/nexuslib/subworkflows/<category>/<tool_dir>/
            ├── <tool_dir>.nf
            └── docs_extras.md    ← appended verbatim to the generated qmd

    The content is appended after the auto-generated sections (banner,
    usage, callout, parameters table), so the extras sit at the bottom
    of the page. Use level-2 (`##`) headings inside the extras file to
    match the surrounding sections.
    """
    extras_path = os.path.join(tool_path, "docs_extras.md")
    if not os.path.exists(extras_path):
        return ""
    with open(extras_path, "r") as f:
        return f.read().strip()


def read_nf_file(nf_path):
    """Read an .nf file and extract banner, help text, known_methods, and params with defaults."""
    if not os.path.exists(nf_path):
        return "", "", [], []
    with open(nf_path, 'r') as f:
        content = f.read()

    # Extract banner (description line between === bars)
    banner = ""
    banner_match = re.search(r'log\.info\s*"""\\\s*\n\s*(=+\n.*?\n\s*=+)', content, re.DOTALL)
    if banner_match:
        for line in banner_match.group(1).strip().split('\n'):
            line = line.strip()
            if line and not line.startswith('='):
                banner = line
                break

    # Extract help text block
    help_text = ""
    help_match = re.search(
        r'if\s*\(params\.help\)\s*\{\s*log\.info\s*"""\\\s*\n(.*?)"""\s*\.stripIndent',
        content, re.DOTALL
    )
    if help_match:
        help_text = help_match.group(1)

    # Extract known_methods list
    methods = []
    methods_match = re.search(r'def known_methods\s*=\s*\[(.*?)\]', content, re.DOTALL)
    if methods_match:
        methods = [m.strip().strip("'\"") for m in methods_match.group(1).split(',')]
        methods = [m for m in methods if m and m != 'all']

    # Extract params with defaults: params.xxx = 'value'
    param_defaults = []
    skip_params = {'help', 'samples_tsv_file', 'output_dir'}
    for match in re.finditer(r"params\.(\w+)\s*=\s*'([^']*)'", content):
        name = match.group(1)
        default = match.group(2)
        if name not in skip_params:
            param_defaults.append((name, default))
    # Also match params.xxx = "" (double quotes)
    for match in re.finditer(r'params\.(\w+)\s*=\s*"([^"]*)"', content):
        name = match.group(1)
        default = match.group(2)
        if name not in skip_params:
            # Avoid duplicates
            if not any(p[0] == name for p in param_defaults):
                param_defaults.append((name, default))

    return banner, help_text, methods, param_defaults


def parse_help_to_params_table(help_text):
    """Parse help text required/optional arguments into a markdown table."""
    if not help_text:
        return ""

    rows = []
    current_param = None
    current_desc = []

    for line in help_text.split('\n'):
        stripped = line.strip()
        # Match parameter lines: --param_name  :  description
        param_match = re.match(r'(--\S+)\s+:\s+(.*)', stripped)
        if param_match:
            # Save previous param
            if current_param:
                rows.append((current_param, ' '.join(current_desc).strip()))
            current_param = param_match.group(1)
            current_desc = [param_match.group(2).strip()]
        elif current_param and stripped and not stripped.startswith('usage:') and not stripped.startswith('required') and not stripped.startswith('optional') and not stripped.startswith('workflow:') and not stripped.startswith('-c ') and not stripped.startswith('-w '):
            # Continuation of previous description
            current_desc.append(stripped)

    # Save last param
    if current_param:
        rows.append((current_param, ' '.join(current_desc).strip()))

    if not rows:
        return ""

    lines = ["| parameter | description |", "|:----------|:------------|"]
    for param, desc in rows:
        # Escape pipes in description
        desc = desc.replace('|', '\\|')
        # Wrap CLI-style --flag tokens in backticks so Quarto's smart-typography
        # doesn't convert the leading "--" into an en-dash. Only match flags
        # that aren't already inside backticks.
        desc = re.sub(r"(?<![`\w])(--[a-zA-Z][\w-]*)(?![`\w])", r"`\1`", desc)
        lines.append(f"| `{param}` | {desc} |")
    lines.append("")
    lines.append(": {.striped .hover}")
    return '\n'.join(lines)


def _placeholder_for_param(name, help_text=""):
    """Return a sensible placeholder based on the parameter name.

    Resolution order:
      1. An explicit per-parameter override, for semantic values no generic
         heuristic can infer (e.g. a tumor purity must lie in (0, 1], so the
         generic float "0.0" is a nonsensical example -> use "0.5").
      2. A suffix match on the param name (e.g. `_fastq_file` ->
         `/path/to/file.fastq`). For gzip-capable formats (fastq/fasta/vcf), if
         `help_text` states the file is compressed (mentions e.g. "fastq.gz"),
         a ".gz" is appended (`/path/to/file.fastq.gz`) so the example matches
         what the tool actually reads/writes.
      3. A help-text scan for a known file-format token, when no suffix
         matches. Here a trailing digit in the param name (e.g. `out_hap1`,
         `out_hap2`) is appended to the example filename so paired outputs
         render as `/path/to/file1.fastq.gz` / `/path/to/file2.fastq.gz`.
      4. Otherwise `""`.
    """
    # (1) Per-parameter overrides for values a heuristic can't infer.
    example_overrides = {
        'tumor_purity': '0.5',
        'sex': 'female',
    }
    if name in example_overrides:
        return example_overrides[name]

    file_extensions = [
        ('_fasta_file', '/path/to/file.fasta'),
        ('_fa_file', '/path/to/file.fa'),
        ('_fastq_file', '/path/to/file.fastq'),
        ('_gtf_file', '/path/to/file.gtf'),
        ('_gff_file', '/path/to/file.gff'),
        ('_bed_file', '/path/to/file.bed'),
        ('_bam_file', '/path/to/file.bam'),
        ('_bai_file', '/path/to/file.bai'),
        ('_vcf_file', '/path/to/file.vcf'),
        ('_tsv_file', '/path/to/file.tsv'),
        ('_csv_file', '/path/to/file.csv'),
        ('_txt_file', '/path/to/file.txt'),
        ('_json_file', '/path/to/file.json'),
        ('_yaml_file', '/path/to/file.yaml'),
        ('_xml_file', '/path/to/file.xml'),
        ('_png_file', '/path/to/file.png'),
        ('_sif_file', '/path/to/file.sif'),
        ('_file', '/path/to/file'),
        ('_dir', '/path/to/dir/'),
        ('_path', '/path/to/dir/'),
    ]
    # (2) Suffix match, refined to ".gz" when the help text says the file is
    # gzip-compressed (only for formats that are commonly gzipped).
    gzippable = ('fastq', 'fasta', 'vcf')
    for suffix, placeholder in file_extensions:
        if name.endswith(suffix):
            ext = placeholder.rsplit('.', 1)[-1]
            if help_text and ext in gzippable and ('%s.gz' % ext) in help_text.lower():
                return placeholder + '.gz'
            return placeholder

    # Help-text fallback. Order matters: longer tokens first so "fastq.gz"
    # wins over the substring "fastq".
    if help_text:
        digit_match = re.search(r'(\d+)$', name)
        digit_suffix = digit_match.group(1) if digit_match else ''
        help_lower = help_text.lower()
        help_extensions = [
            ('fastq.gz', f'/path/to/file{digit_suffix}.fastq.gz'),
            ('fasta.gz', f'/path/to/file{digit_suffix}.fasta.gz'),
            ('vcf.gz',   f'/path/to/file{digit_suffix}.vcf.gz'),
            ('fastq',    f'/path/to/file{digit_suffix}.fastq'),
            ('fasta',    f'/path/to/file{digit_suffix}.fasta'),
            ('.bam',     f'/path/to/file{digit_suffix}.bam'),
            ('.vcf',     f'/path/to/file{digit_suffix}.vcf'),
            ('.tsv',     f'/path/to/file{digit_suffix}.tsv'),
            ('.gtf',     f'/path/to/file{digit_suffix}.gtf'),
            ('.bed',     f'/path/to/file{digit_suffix}.bed'),
        ]
        for keyword, placeholder in help_extensions:
            if keyword in help_lower:
                return placeholder

    return '""'


def generate_subworkflow_page(tool_name, nf_filename, nf_path, tool_dir, tool_path, cat_docs_dir):
    """Generate a .qmd page for a subworkflow from its .nf file."""
    banner, help_text, _, param_defaults = read_nf_file(nf_path)
    params_path = os.path.join(tool_path, "params.yaml")
    has_params_yaml = os.path.exists(params_path)

    lines = [
        f'---',
        f'title: "{tool_name}"',
        f'---',
        '',
    ]

    if banner:
        lines.append(banner)
        lines.append('')

    # Usage section
    lines.append('## usage')
    lines.append('')
    if has_params_yaml:
        lines.append('```bash')
        lines.append(f'nexus run --nf-workflow {nf_filename} -params-file params.yaml')
        lines.append('```')
    else:
        lines.append('```bash')
        lines.append(f'nexus run --nf-workflow {nf_filename} \\')
        lines.append(f'    -c nextflow.config \\')
        lines.append(f'    -w work/ \\')
        lines.append(f'    --samples_tsv_file samples.tsv \\')
        lines.append(f'    --output_dir results/ \\')

        # Add all parameters with their defaults
        for i, (name, default) in enumerate(param_defaults):
            if default:
                val = f'"{default}"'
            else:
                val = _placeholder_for_param(name)
            is_last = (i == len(param_defaults) - 1)
            if is_last:
                lines.append(f'    --{name} {val}')
            else:
                lines.append(f'    --{name} {val} \\')

        if not param_defaults:
            # Remove trailing backslash from --output_dir line
            lines[-1] = lines[-1].rstrip(' \\')

        lines.append('```')
    lines.append('')

    # Config file note
    lines.append('::: {.callout-note}')
    lines.append('Nextflow config files are available')
    lines.append('[here](https://github.com/pirl-unc/nexus/tree/main/nextflow).')
    lines.append('Use the config file that matches your installed Nexus version (e.g. `nexus_v0.2.0_nextflow_slurm.config`).')
    lines.append(':::')
    lines.append('')

    # Parameters section
    if has_params_yaml:
        params_dest_name = f"{tool_dir}_params.yaml"
        shutil.copy2(params_path, os.path.join(cat_docs_dir, params_dest_name))
        with open(params_path, 'r') as f:
            params_content = f.read()
        lines.append('## parameters')
        lines.append('')
        lines.append(f'[Download params.yaml]({params_dest_name}){{.btn .btn-primary download="{params_dest_name}"}}')
        lines.append('')
        lines.append('```yaml')
        lines.append(params_content.strip())
        lines.append('```')
        lines.append('')
    else:
        params_table = parse_help_to_params_table(help_text)
        if params_table:
            lines.append('## parameters')
            lines.append('')
            lines.append(params_table)
            lines.append('')

    # Hand-written extras (see read_extras_block docstring for usage).
    # Appended verbatim — author is responsible for headings/formatting.
    extras = read_extras_block(tool_path)
    if extras:
        lines.append(extras)
        lines.append('')

    return '\n'.join(lines)


def generate_subworkflow_pages():
    """Generate .qmd pages for all subworkflows."""
    sidebar_items = {}

    for category in sorted(os.listdir(SUBWORKFLOWS_DIR)):
        category_path = os.path.join(SUBWORKFLOWS_DIR, category)
        if not is_source_dir(category_path, category):
            continue

        label = CATEGORY_LABELS.get(category, category.replace('_', ' '))
        sidebar_items[category] = {"label": label, "tools": []}

        cat_docs_dir = os.path.join(DOCS_DIR, "subworkflows", category)
        os.makedirs(cat_docs_dir, exist_ok=True)

        for tool_dir in sorted(os.listdir(category_path)):
            tool_path = os.path.join(category_path, tool_dir)
            if not is_source_dir(tool_path, tool_dir):
                continue

            tool_name = extract_tool_name(tool_dir)
            nf_filename, nf_path = resolve_nf_path(tool_path, tool_dir)

            qmd_content = generate_subworkflow_page(tool_name, nf_filename, nf_path, tool_dir, tool_path, cat_docs_dir)
            qmd_path = os.path.join(cat_docs_dir, f"{tool_dir}.qmd")

            with open(qmd_path, 'w') as f:
                f.write(qmd_content)

            sidebar_items[category]["tools"].append({
                "name": tool_name,
                "file": f"subworkflows/{category}/{tool_dir}.qmd"
            })

    return sidebar_items


def generate_workflow_pages():
    """Generate .qmd pages for all workflows."""
    sidebar_items = {}

    for category in sorted(os.listdir(WORKFLOWS_DIR)):
        category_path = os.path.join(WORKFLOWS_DIR, category)
        if not is_source_dir(category_path, category):
            continue

        label = CATEGORY_LABELS.get(category, category.replace('_', ' '))
        sidebar_items[category] = {"label": label, "workflows": []}

        cat_docs_dir = os.path.join(DOCS_DIR, "workflows", category)
        os.makedirs(cat_docs_dir, exist_ok=True)

        for wf_dir in sorted(os.listdir(category_path)):
            wf_path = os.path.join(category_path, wf_dir)
            if not is_source_dir(wf_path, wf_dir):
                continue

            wf_label = WORKFLOW_LABELS.get(wf_dir, wf_dir.replace('_', ' ').replace('-', ' '))
            nf_filename, nf_path = resolve_nf_path(wf_path, wf_dir)
            params_path = os.path.join(wf_path, "params.yaml")

            banner, help_text, methods, _ = read_nf_file(nf_path)

            params_content = ""
            if os.path.exists(params_path):
                with open(params_path, 'r') as f:
                    params_content = f.read()

            qmd_lines = [
                f'---',
                f'title: "{wf_label}"',
                f'---',
                f'',
            ]

            if banner:
                qmd_lines.append(banner)
                qmd_lines.append('')

            # Tools section
            if methods:
                qmd_lines.append('## tools')
                qmd_lines.append('')
                qmd_lines.append('The following tools run by default (`methods: "all"`):')
                qmd_lines.append('')
                qmd_lines.append(', '.join(f'`{m}`' for m in methods))
                qmd_lines.append('')

            # Usage section
            qmd_lines.append('## usage')
            qmd_lines.append('')
            qmd_lines.append('```bash')
            qmd_lines.append(f'nexus run --nf-workflow {nf_filename} -params-file params.yaml')
            qmd_lines.append('```')
            qmd_lines.append('')

            # Config file note
            qmd_lines.append('::: {.callout-note}')
            qmd_lines.append('Nextflow config files are available')
            qmd_lines.append('[here](https://github.com/pirl-unc/nexus/tree/main/nextflow).')
            qmd_lines.append('Use the config file that matches your installed nexus version (e.g. `nexus_v0.2.0_nextflow_slurm.config`).')
            qmd_lines.append(':::')
            qmd_lines.append('')

            # Parameters section
            if params_content:
                # Copy params.yaml into docs dir for download
                params_dest_name = f"{wf_dir}_params.yaml"
                shutil.copy2(params_path, os.path.join(cat_docs_dir, params_dest_name))

                qmd_lines.append('## parameters')
                qmd_lines.append('')
                qmd_lines.append(f'[Download params.yaml]({params_dest_name}){{.btn .btn-primary download="{params_dest_name}"}}')
                qmd_lines.append('')
                qmd_lines.append('```yaml')
                qmd_lines.append(params_content.strip())
                qmd_lines.append('```')
                qmd_lines.append('')

            # Hand-written extras (same mechanism as subworkflows).
            extras = read_extras_block(wf_path)
            if extras:
                qmd_lines.append(extras)
                qmd_lines.append('')

            qmd_path = os.path.join(cat_docs_dir, f"{wf_dir}.qmd")
            with open(qmd_path, 'w') as f:
                f.write('\n'.join(qmd_lines) + '\n')

            sidebar_items[category]["workflows"].append({
                "name": wf_label,
                "file": f"workflows/{category}/{wf_dir}.qmd"
            })

    return sidebar_items


def read_cli_commands():
    """Read pyproject.toml to get CLI command -> module mappings."""
    cli_map = {}  # module_name -> cli_command
    if not os.path.exists(PYPROJECT_PATH):
        return cli_map
    with open(PYPROJECT_PATH, 'r') as f:
        content = f.read()
    in_section = False
    for line in content.split('\n'):
        stripped = line.strip()
        if stripped == '[project.scripts]':
            in_section = True
            continue
        if in_section:
            if stripped.startswith('['):
                break
            match = re.match(r'(\S+)\s*=\s*"nexuslib\.scripts\.(\w+):run"', stripped)
            if match:
                cli_map[match.group(2)] = match.group(1)
    return cli_map


def parse_utility_argparse(py_path):
    """Parse a utility script to extract argparse description and arguments using AST."""
    with open(py_path, 'r') as f:
        source = f.read()

    tree = ast.parse(source)
    description = ""
    arguments = []

    for node in ast.walk(tree):
        # Find ArgumentParser(description=...)
        if isinstance(node, ast.Call):
            func = node.func
            is_argparser = False
            if isinstance(func, ast.Attribute) and func.attr == 'ArgumentParser':
                is_argparser = True
            elif isinstance(func, ast.Name) and func.id == 'ArgumentParser':
                is_argparser = True
            if is_argparser:
                for kw in node.keywords:
                    if kw.arg == 'description' and isinstance(kw.value, ast.Constant):
                        description = ' '.join(kw.value.value.split())

        # Find add_argument(...) calls
        if isinstance(node, ast.Call):
            func = node.func
            if isinstance(func, ast.Attribute) and func.attr == 'add_argument':
                arg_name = ""
                help_text = ""
                default = None
                arg_type = ""
                required = False
                choices = []

                # Get positional arg name(s)
                for arg in node.args:
                    if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                        if arg.value.startswith('-'):
                            arg_name = arg.value

                # Get keyword arguments
                for kw in node.keywords:
                    if kw.arg == 'help' and isinstance(kw.value, ast.Constant):
                        help_text = ' '.join(kw.value.value.split())
                    elif kw.arg == 'default' and isinstance(kw.value, ast.Constant):
                        default = kw.value.value
                    elif kw.arg == 'type' and isinstance(kw.value, ast.Name):
                        arg_type = kw.value.id
                    elif kw.arg == 'required' and isinstance(kw.value, ast.Constant):
                        required = kw.value.value
                    elif kw.arg == 'choices' and isinstance(kw.value, (ast.List, ast.Tuple)):
                        choices = [e.value for e in kw.value.elts
                                   if isinstance(e, ast.Constant)]

                if arg_name:
                    arguments.append({
                        'name': arg_name,
                        'help': help_text,
                        'default': default,
                        'type': arg_type,
                        'required': required,
                        'choices': choices,
                    })

    return description, arguments


# Per-script Usage examples, keyed by module name then by the exact flag string.
# Use this only when the generic placeholder heuristics can't (or shouldn't)
# produce the value: e.g. short-flag tools whose flag (`-i`/`-o`) carries no
# descriptive name to drive a suffix match, or a bespoke output path that must
# not leak into other tools that share the same flag/dest name (e.g.
# `output_tsv_file` is also used by merge_rna_assemblies and the CCF tool).
USAGE_EXAMPLE_OVERRIDES = {
    'convert_netmhcpan_txt2tsv': {
        '-i': '/path/to/file.txt',
        '-o': '/path/to/file/output.tsv',
    },
    'create_abra2_targets_bed_file': {
        '--bedtools': '/path/to/bedtools.static.binary',
    },
    'create_beers2_input_data': {
        '--sample-id': 'sample001',
    },
}


def generate_utility_pages():
    """Generate .qmd pages for all utilities."""
    utilities_docs_dir = os.path.join(DOCS_DIR, "utilities")
    os.makedirs(utilities_docs_dir, exist_ok=True)

    cli_map = read_cli_commands()
    sidebar_items = []

    for filename in sorted(os.listdir(UTILITIES_DIR)):
        if not filename.endswith('.py') or filename == '__init__.py':
            continue

        module_name = filename[:-3]  # strip .py
        py_path = os.path.join(UTILITIES_DIR, filename)
        description, arguments = parse_utility_argparse(py_path)

        cli_command = cli_map.get(module_name, module_name.replace('_', '-'))
        title = cli_command

        lines = [
            '---',
            f'title: "{title}"',
            '---',
            '',
        ]

        if description:
            lines.append(description)
            lines.append('')

        # Usage
        lines.append('## Usage')
        lines.append('')
        lines.append('```bash')
        usage_parts = [cli_command]
        script_examples = USAGE_EXAMPLE_OVERRIDES.get(module_name, {})
        for arg in arguments:
            if arg['name'] in script_examples:
                # Per-(script, flag) example for cases no generic heuristic can
                # infer (short-flag tools, or a bespoke path that must not leak
                # into other tools sharing the same flag/dest name).
                usage_parts.append(f"{arg['name']} {script_examples[arg['name']]}")
            elif arg['default'] is not None:
                usage_parts.append(f"[{arg['name']} {arg['default']}]")
            else:
                # Convert --arg-name to arg_name for placeholder lookup. Pass
                # the help text in so paired outputs like --out-hap1/--out-hap2
                # can be inferred from their help strings ("Output fastq.gz
                # for haplotype 1." → /path/to/file1.fastq.gz).
                param_key = arg['name'].lstrip('-').replace('-', '_')
                placeholder = _placeholder_for_param(param_key, arg['help'])
                # For non-file params, use type-aware placeholder
                if placeholder == '""' and arg['type']:
                    type_placeholders = {'int': '0', 'float': '0.0', 'str': '""'}
                    placeholder = type_placeholders.get(arg['type'], '""')
                usage_parts.append(f"{arg['name']} {placeholder}")
        lines.append(' \\\n    '.join(usage_parts))
        lines.append('```')
        lines.append('')

        # Parameters table
        if arguments:
            lines.append('## Parameters')
            lines.append('')
            lines.append('| Parameter | Type | Default | Description |')
            lines.append('|:----------|:-----|:--------|:------------|')
            for arg in arguments:
                default_str = f'`{arg["default"]}`' if arg['default'] is not None else 'required'
                type_str = f'`{arg["type"]}`' if arg['type'] else ''
                help_str = arg['help'].replace('|', '\\|') if arg['help'] else ''
                # Surface argparse choices so enum-style params (e.g. --sex,
                # --missing-from-bam) document their allowed values in the table.
                if arg.get('choices'):
                    choice_str = ' or '.join(f'`{c}`' for c in arg['choices'])
                    help_str = (f"{help_str} " if help_str else '') + f"Choices: {choice_str}."
                lines.append(f"| `{arg['name']}` | {type_str} | {default_str} | {help_str} |")
            lines.append('')
            lines.append(': {.striped .hover}')
            lines.append('')

        qmd_path = os.path.join(utilities_docs_dir, f"{module_name}.qmd")
        with open(qmd_path, 'w') as f:
            f.write('\n'.join(lines))

        sidebar_items.append({
            "name": title,
            "file": f"utilities/{module_name}.qmd"
        })

    # Generate utilities/index.qmd
    index_lines = [
        '---',
        'title: "Utilities"',
        '---',
        '',
        'Standalone command-line utilities installed with Nexus.',
        '',
        '| Command | Description |',
        '|:--------|:------------|',
    ]
    for filename in sorted(os.listdir(UTILITIES_DIR)):
        if not filename.endswith('.py') or filename == '__init__.py':
            continue
        module_name = filename[:-3]
        py_path = os.path.join(UTILITIES_DIR, filename)
        description, _ = parse_utility_argparse(py_path)
        cli_command = cli_map.get(module_name, module_name.replace('_', '-'))
        desc_short = description.split('.')[0] + '.' if description else ''
        index_lines.append(f"| [{cli_command}]({module_name}.qmd) | {desc_short} |")
    index_lines.append('')
    index_lines.append(': {.striped .hover tbl-colwidths="[40, 60]"}')

    index_path = os.path.join(DOCS_DIR, "utilities", "index.qmd")
    with open(index_path, 'w') as f:
        f.write('\n'.join(index_lines))

    return sidebar_items


def generate_quarto_yml(subworkflow_sidebar, workflow_sidebar, utility_sidebar):
    """Generate _quarto.yml with full navigation."""
    # Sort each category's contents by displayed tool/workflow name so the
    # sidebar order matches the user-facing label, not the (possibly
    # typo'd or category-prefixed) source directory name. Without this, a
    # directory like `sequenicng_simulation_beers2/` sorts after
    # `sequencing_simulation_*/` and `beers2` ends up at the bottom of the
    # Sequencing Simulation section even though `b` < `n` alphabetically.
    sw_contents = []
    for category in sorted(subworkflow_sidebar.keys()):
        info = subworkflow_sidebar[category]
        tools_sorted = sorted(info["tools"], key=lambda t: t["name"].lower())
        sw_contents.append({
            "section": info["label"],
            "contents": [t["file"] for t in tools_sorted]
        })

    wf_contents = []
    for category in sorted(workflow_sidebar.keys()):
        info = workflow_sidebar[category]
        workflows_sorted = sorted(info["workflows"], key=lambda w: w["name"].lower())
        wf_contents.append({
            "section": info["label"],
            "contents": [w["file"] for w in workflows_sorted]
        })

    sc_contents = [s["file"] for s in utility_sidebar]

    quarto_config = {
        "project": {"type": "website", "output-dir": "_site"},
        "website": {
            "title": "Nexus",
            "navbar": {
                "left": [
                    {"text": "Home", "file": "index.qmd"},
                    {"text": "Subworkflows", "file": "subworkflows/index.qmd"},
                    {"text": "Workflows", "file": "workflows/index.qmd"},
                    {"text": "Utilities", "file": "utilities/index.qmd"},
                    {"text": "FAQ", "file": "faq.qmd"},
                ],
                "right": [
                    {"icon": "github", "href": "https://github.com/pirl-unc/nexus"}
                ]
            },
            "sidebar": [
                {
                    "id": "subworkflows", "title": "Subworkflows", "style": "floating",
                    "contents": [{"section": "Overview", "contents": ["subworkflows/index.qmd"]}] + sw_contents
                },
                {
                    "id": "workflows", "title": "Workflows", "style": "floating",
                    "contents": [{"section": "Overview", "contents": ["workflows/index.qmd"]}] + wf_contents
                },
                {
                    "id": "utilities", "title": "Utilities", "style": "floating",
                    "contents": [{"section": "Overview", "contents": ["utilities/index.qmd"]}] + [{"section": "Commands", "contents": sc_contents}]
                }
            ]
        },
        # `custom.scss` lives next to _quarto.yml under docs/. Listing it here
        # (as a list, not a string) tells Quarto to layer it on top of the
        # cosmo Bootstrap theme.
        "format": {"html": {"theme": ["cosmo", "custom.scss"], "toc": True}}
    }

    with open(os.path.join(DOCS_DIR, "_quarto.yml"), 'w') as f:
        yaml.dump(quarto_config, f, default_flow_style=False, sort_keys=False)


def main():
    print("Generating Quarto documentation site...")

    print("  Generating subworkflow pages...")
    sw_sidebar = generate_subworkflow_pages()
    print(f"    Generated {sum(len(v['tools']) for v in sw_sidebar.values())} subworkflow pages")

    print("  Generating workflow pages...")
    wf_sidebar = generate_workflow_pages()
    print(f"    Generated {sum(len(v['workflows']) for v in wf_sidebar.values())} workflow pages")

    print("  Generating utility pages...")
    sc_sidebar = generate_utility_pages()
    print(f"    Generated {len(sc_sidebar)} utility pages")

    print("  Generating _quarto.yml...")
    generate_quarto_yml(sw_sidebar, wf_sidebar, sc_sidebar)

    print(f"Done! Site generated in {DOCS_DIR}/")
    print(f"  NOTE: index.qmd, subworkflows/index.qmd, workflows/index.qmd are NOT overwritten.")
    print(f"  Preview with: cd docs && quarto preview")


if __name__ == "__main__":
    main()
