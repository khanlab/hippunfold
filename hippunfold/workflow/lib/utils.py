import logging
import sys


def setup_logger(log_file=None):
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    handler = (
        logging.FileHandler(log_file) if log_file else logging.StreamHandler(sys.stderr)
    )
    handler.setFormatter(formatter)

    logger.addHandler(handler)
    return logger
