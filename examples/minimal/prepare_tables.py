"""Unpack the CCGeant4 example data using the Python standard library."""
import argparse
import shutil
import json
import lzma
import os
from pathlib import Path
import tempfile

EXAMPLE = Path(__file__).resolve().parent


def unpack(destination):
    manifest = json.loads((EXAMPLE / "SOURCE_DATA.json").read_text(encoding="utf-8"))
    for entry in manifest["files"]:
        destination_file = destination / entry["file"]
        if destination_file.is_file() and destination_file.stat().st_size == entry["bytes"]:
            continue
        archive = EXAMPLE / entry["archive"]
        if not archive.is_file():
            raise FileNotFoundError(f"Missing example data: {archive}")
        destination_file.parent.mkdir(parents=True, exist_ok=True)
        temporary = None
        try:
            with lzma.open(archive, "rb") as source, tempfile.NamedTemporaryFile(
                dir=destination_file.parent, delete=False
            ) as target:
                temporary = Path(target.name)
                shutil.copyfileobj(source, target, 1024 * 1024)
            if temporary.stat().st_size != entry["bytes"]:
                raise ValueError(f"Unexpected decompressed size for {entry['archive']}")
            os.replace(temporary, destination_file)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
    return destination


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    print(unpack(parser.parse_args().destination.resolve()))
