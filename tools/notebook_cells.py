#!/usr/bin/env python3
"""Edit notebook source one cell at a time without exposing saved outputs."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import os


def digest(data):
    return hashlib.sha256(data).hexdigest()


def source(cell):
    value = cell['source']
    return value if isinstance(value, str) else ''.join(value)


def filename(index, cell):
    extension = {'code': 'py', 'markdown': 'md', 'raw': 'txt'}[cell['cell_type']]
    return f'{index:04d}.{extension}'


def read_notebook(data):
    notebook = json.loads(data)
    if notebook.get('nbformat') != 4 or not isinstance(notebook.get('cells'), list):
        raise ValueError('Expected a version 4 notebook with a cells list')
    for cell in notebook['cells']:
        filename(0, cell)
        source(cell)
    return notebook


def split(notebook_path, workspace):
    data = notebook_path.read_bytes()
    notebook = read_notebook(data)
    workspace.mkdir(parents=True, exist_ok=False)
    (workspace / 'original.ipynb').write_bytes(data)
    manifest = {'notebook': str(notebook_path.resolve()), 'sha256': digest(data)}
    (workspace / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    cells_dir = workspace / 'cells'
    cells_dir.mkdir()
    index = []
    for i, cell in enumerate(notebook['cells']):
        name = filename(i, cell)
        text = source(cell)
        (cells_dir / name).write_bytes(text.encode('utf-8'))
        index.append(f'{name}\t{len(text)} characters\t{len(text.splitlines())} lines')
    (workspace / 'index.txt').write_text('\n'.join(index) + '\n')
    print(f'Extracted {len(index)} cells; original notebook unchanged. Index: {workspace / "index.txt"}')


def build(workspace, allow_code=False):
    manifest = json.loads((workspace / 'manifest.json').read_text())
    data = (workspace / 'original.ipynb').read_bytes()
    if digest(data) != manifest['sha256']:
        raise ValueError('Backup changed; refusing to rebuild')
    notebook_path = Path(manifest['notebook'])
    if digest(notebook_path.read_bytes()) != manifest['sha256']:
        raise ValueError('Original notebook changed since split; start a fresh workspace')
    original = read_notebook(data)
    result = copy.deepcopy(original)
    expected = {filename(i, c) for i, c in enumerate(original['cells'])}
    actual = {p.name for p in (workspace / 'cells').iterdir()}
    if actual != expected:
        raise ValueError('Cell files missing or added; refusing to rebuild')
    changes = []
    for i, cell in enumerate(result['cells']):
        name = filename(i, cell)
        text = (workspace / 'cells' / name).read_bytes().decode('utf-8')
        if text == source(cell):
            continue
        if cell['cell_type'] == 'code' and not allow_code:
            raise ValueError(f'{name}: code changed; inspect it and explicitly use --allow-code')
        cell['source'] = text if isinstance(cell['source'], str) else text.splitlines(keepends=True)
        if cell['cell_type'] == 'code':
            cell['outputs'] = []
            cell['execution_count'] = None
        changes.append(name)
    # All fields other than intentional source edits and stale code outputs survive.
    encoded = data if not changes else (json.dumps(result, ensure_ascii=False, indent=1) + '\n').encode('utf-8')
    if read_notebook(encoded) != result:
        raise ValueError('Round-trip verification failed')
    return notebook_path, encoded, changes


def apply(workspace, allow_code=False):
    path, data, changes = build(workspace, allow_code)
    if not changes:
        print('No changes; original bytes preserved')
        return
    # Write alongside the original and replace only after the complete file exists.
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        temporary.chmod(path.stat().st_mode)
        # Recheck immediately before replacement; use only one worker per notebook.
        manifest = json.loads((workspace / 'manifest.json').read_text())
        if digest(path.read_bytes()) != manifest['sha256']:
            raise ValueError('Original changed during apply; refusing to replace it')
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    print('Applied: ' + ', '.join(changes))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    extract = sub.add_parser('split')
    extract.add_argument('notebook', type=Path)
    extract.add_argument('workspace', type=Path)
    for command in ('check', 'apply'):
        item = sub.add_parser(command)
        item.add_argument('workspace', type=Path)
        item.add_argument('--allow-code', action='store_true')
    args = parser.parse_args()
    try:
        if args.command == 'split':
            split(args.notebook, args.workspace)
        elif args.command == 'check':
            _, _, changes = build(args.workspace, args.allow_code)
            print('Changed cells: ' + (', '.join(changes) or 'none'))
        else:
            apply(args.workspace, args.allow_code)
    except (ValueError, KeyError, OSError) as error:
        parser.exit(1, f'{error}\n')


if __name__ == '__main__':
    main()
