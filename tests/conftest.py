import sys
from pathlib import Path

# make hippunfold/workflow/lib importable as `lib`, as it is for workflow scripts
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "hippunfold" / "workflow"))
