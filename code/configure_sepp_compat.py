#!/usr/bin/env python3
"""Point SEPP's pplacer interface to the EPA-ng compatibility adapter."""

from __future__ import annotations

import argparse
import configparser
from pathlib import Path
import shutil
import tempfile


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--adapter", required=True, type=Path)
    parser.add_argument("--record", required=True, type=Path)
    args = parser.parse_args()

    config_path = args.config.resolve()
    adapter = args.adapter.resolve()
    record = args.record.resolve()
    if not config_path.is_file():
        raise FileNotFoundError(f"SEPP main.config not found: {config_path}")
    if not adapter.is_file():
        raise FileNotFoundError(f"Compatibility adapter not found: {adapter}")

    parser_obj = configparser.ConfigParser()
    parser_obj.optionxform = str
    parser_obj.read(config_path)
    if "pplacer" not in parser_obj or "path" not in parser_obj["pplacer"]:
        raise ValueError(f"Unexpected SEPP configuration structure: {config_path}")

    backup = config_path.with_name(config_path.name + ".before_epa_ng_compat")
    if not backup.exists():
        shutil.copy2(config_path, backup)

    parser_obj["pplacer"]["path"] = str(adapter)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=config_path.parent,
        prefix=config_path.name + ".",
        suffix=".tmp",
        delete=False,
    ) as handle:
        parser_obj.write(handle)
        temp_path = Path(handle.name)
    temp_path.replace(config_path)

    record.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(config_path, record)
    print(f"Configured SEPP pplacer interface: {adapter}")


if __name__ == "__main__":
    main()
