import csv

from lib.qc import QCResult, verify
from lib.utils import setup_logger

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)

results = [
    QCResult.load(path, warn_checks=snakemake.params.warn_checks)
    for path in snakemake.input.metrics
]
rows = verify(
    results,
    warn_only=snakemake.params.warn_only,
    logger=logger,
    report=snakemake.input.report,
)

columns = list(dict.fromkeys(k for row in rows for k in row))
with open(snakemake.output.tsv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

logger.info(f"{len(rows)} QC checks written to {snakemake.output.tsv}")
