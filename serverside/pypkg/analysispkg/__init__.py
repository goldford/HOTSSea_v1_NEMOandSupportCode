from pathlib import Path

pyfiles = Path.glob(Path(__file__).parent, 'pkg_*.py')
modules = sorted([p.stem for p in pyfiles])

__all__ = modules
