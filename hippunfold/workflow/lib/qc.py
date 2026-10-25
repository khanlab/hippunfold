"""
QC metrics and checks.

Each QC step writes one or more QC result JSON files, typically from the rule
that computes the thing being checked, so the metrics are recorded without
failing that rule. A report rule then visualizes the results, and a separate
verify rule (which also requires the report, so it always exists) raises an
error if any error-severity check failed.

A QC result file has the form:

    {
      "qc": "srcsink",                       # name of the QC step
      "entities": {"subject": "001", ...},   # what was checked
      "metrics": {"name": {"value": ..., "description": ...}, ...},
      "checks": [
        {"name": ..., "passed": bool, "severity": "error" | "warning",
         "message": ..., "description": ...},
        ...
      ],
      "data": {...}                          # arrays etc for plotting
    }
"""

import json

import numpy as np

ERROR = "error"
WARNING = "warning"
SEVERITIES = (ERROR, WARNING)


class QCCheckError(Exception):
    """Raised when error-severity QC checks fail."""


def _to_builtin(obj):
    """Convert numpy types (recursively) to JSON-serializable builtins."""
    if isinstance(obj, dict):
        return {str(k): _to_builtin(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_builtin(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _to_builtin(obj.tolist())
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


class QCResult:
    """Metrics, checks and plotting data for one QC step on one input."""

    def __init__(self, qc, entities=None):
        self.qc = qc
        self.entities = dict(entities or {})
        self.metrics = {}
        self.checks = []
        self.data = {}

    def add_metric(self, name, value, description=""):
        self.metrics[name] = {"value": _to_builtin(value), "description": description}
        return value

    def add_check(self, name, passed, message="", severity=ERROR, description=""):
        if severity not in SEVERITIES:
            raise ValueError(f"severity must be one of {SEVERITIES}, not {severity}")
        self.checks.append(
            {
                "name": name,
                "passed": bool(passed),
                "severity": severity,
                "message": message,
                "description": description,
            }
        )
        return bool(passed)

    def add_data(self, name, value):
        self.data[name] = _to_builtin(value)

    def failed(self, severity=None):
        return [
            c
            for c in self.checks
            if not c["passed"] and (severity is None or c["severity"] == severity)
        ]

    @property
    def status(self):
        """'fail' if an error check failed, 'warn' if a warning did, else 'pass'."""
        if self.failed(ERROR):
            return "fail"
        if self.failed(WARNING):
            return "warn"
        return "pass"

    @property
    def label(self):
        """Short human-readable name for the entities, e.g. 'hemi-L label-hipp'."""
        return " ".join(f"{k}-{v}" for k, v in self.entities.items())

    def to_dict(self):
        return {
            "qc": self.qc,
            "entities": self.entities,
            "metrics": self.metrics,
            "checks": self.checks,
            "data": self.data,
        }

    def save(self, path):
        with open(path, "w") as f:
            json.dump(_to_builtin(self.to_dict()), f, indent=1)

    @classmethod
    def load(cls, path):
        with open(path) as f:
            d = json.load(f)
        result = cls(d["qc"], d.get("entities"))
        result.metrics = d.get("metrics", {})
        result.checks = d.get("checks", [])
        result.data = d.get("data", {})
        return result


def overall_status(results):
    statuses = {r.status for r in results}
    for status in ("fail", "warn"):
        if status in statuses:
            return status
    return "pass"


def verify(results, warn_only=False, logger=None, report=None):
    """
    Raise QCCheckError if any error-severity check failed (unless warn_only).
    Failed warnings, and failed errors when warn_only, are logged.

    Returns rows (dicts) summarizing every check, e.g. for writing a TSV.
    """
    rows = []
    errors = []
    for result in results:
        for check in result.checks:
            rows.append(
                {
                    "qc": result.qc,
                    **result.entities,
                    "check": check["name"],
                    "severity": check["severity"],
                    "passed": check["passed"],
                    "message": check["message"],
                }
            )
            if check["passed"]:
                continue
            msg = f"[{result.qc}] {result.label}: {check['name']}: {check['message']}"
            if check["severity"] == ERROR:
                errors.append(msg)
            if logger:
                logger.warning(f"QC {check['severity']}: {msg}")

    if errors and not warn_only:
        where = f"\nSee the QC report for details: {report}" if report else ""
        raise QCCheckError(
            f"{len(errors)} QC check(s) failed:\n  "
            + "\n  ".join(errors)
            + where
            + "\nTo continue anyway, re-run with --qc_warn_only."
        )
    return rows
