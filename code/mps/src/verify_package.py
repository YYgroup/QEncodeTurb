"""Check frozen sources and supplied data, without numerical recomputation."""
import hashlib
import json
from pathlib import Path
from paths import PACKAGE, REFERENCE

def main():
    manifest = json.loads((PACKAGE/'source_manifest.json').read_text(encoding='utf-8'))
    for relative, expected in manifest.items():
        if hashlib.sha256((PACKAGE/relative).read_bytes()).hexdigest() != expected:
            raise RuntimeError('Source checksum mismatch: '+relative)
    saved = json.loads((REFERENCE/'recheck/independent_validation.json').read_text(encoding='utf-8'))
    assert saved['passed'] and saved['summary']['checked_saved_mps'] == 45
    plot_manifest = json.loads((REFERENCE/'figure/data_manifest.json').read_text(encoding='utf-8'))
    for relative, expected in plot_manifest.items():
        assert hashlib.sha256((REFERENCE/'figure'/relative).read_bytes()).hexdigest() == expected
    package_manifest = PACKAGE/'package_manifest.json'
    if package_manifest.exists():
        for relative, expected in json.loads(package_manifest.read_text(encoding='utf-8')).items():
            assert hashlib.sha256((PACKAGE/relative).read_bytes()).hexdigest()==expected, relative
    print(json.dumps(dict(passed=True,frozen_source_files=len(manifest),
                          plotting_data_files=len(plot_manifest),
                          supplied_reference_approximations=45),indent=2))

if __name__ == '__main__':
    main()
