"""Stage the WCRP PNW code mirror in a clean destination checkout."""

import os
import subprocess
import sys
from pathlib import Path

SOURCE_FOLDER = "Code"
TARGET_FOLDER = "WCRP_PNW_heatwave_2021"


def git(repository, *args):
    return subprocess.check_output(["git", "--literal-pathspecs", "-C", repository, *args])


def read_tree(repository, folder):
    entries = {}
    for entry in git(repository, "ls-tree", "-rz", "HEAD", "--", folder).split(b"\0"):
        if not entry:
            continue
        metadata, path = entry.decode("utf-8").split("\t", 1)
        mode, kind, object_id = metadata.split()
        if kind != "blob" or mode not in {"100644", "100755"}:
            raise ValueError(f"Unsupported Git entry: {path} ({mode} {kind})")
        entries[path.removeprefix(f"{folder}/")] = (mode, object_id)
    return entries


def stage_mirror(source, target):
    if source == target or git(target, "status", "--porcelain", "--untracked-files=all"):
        raise ValueError("Use a separate, clean target checkout.")
    if git(target, "ls-files", "--others", "--ignored", "--exclude-standard", "--", TARGET_FOLDER):
        raise ValueError(f"The target checkout contains ignored files under {TARGET_FOLDER}/.")

    source_files = {
        path: entry
        for path, entry in read_tree(source, SOURCE_FOLDER).items()
        if not path.lower().endswith(".bak")
    }
    if not source_files:
        raise ValueError(f"No source files found under {SOURCE_FOLDER}/.")

    target_files = read_tree(target, TARGET_FOLDER)
    for path in sorted(target_files.keys() - source_files.keys()):
        git(target, "rm", "--", f"{TARGET_FOLDER}/{path}")

    for path, (mode, object_id) in source_files.items():
        if target_files.get(path) == (mode, object_id):
            continue
        target_path = f"{TARGET_FOLDER}/{path}"
        output = target / target_path
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(git(source, "cat-file", "blob", object_id))
        executable = "+x" if mode == "100755" else "-x"
        if os.name != "nt":
            output.chmod(0o755 if mode == "100755" else 0o644)
        git(target, "add", "--chmod=" + executable, "--", target_path)

    staged = {}
    for entry in git(target, "ls-files", "--stage", "-z", "--", TARGET_FOLDER).split(b"\0"):
        if entry:
            metadata, path = entry.decode("utf-8").split("\t", 1)
            mode, object_id, stage = metadata.split()
            if stage != "0":
                raise ValueError(f"Unmerged index entry: {path}")
            staged[path.removeprefix(f"{TARGET_FOLDER}/")] = (mode, object_id)
    if staged != source_files:
        raise ValueError("The staged target tree differs from the source Git blobs or modes.")


if __name__ == "__main__":
    source_repo, target_repo = (Path(arg).resolve(strict=True) for arg in sys.argv[1:])
    stage_mirror(source_repo, target_repo)
