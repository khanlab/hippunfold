import logging
import sys
import os
from appdirs import AppDirs


def setup_logger(log_file=None):
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    if log_file:
        handler = logging.FileHandler(log_file)
    else:
        handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    if log_file:
        # Redirect stdout and stderr to the same log file
        sys.stdout = open(log_file, "a")
        sys.stderr = open(log_file, "a")

    return logger


def get_download_dir():
    if "HIPPUNFOLD_CACHE_DIR" in os.environ.keys():
        download_dir = os.environ["HIPPUNFOLD_CACHE_DIR"]
    else:
        # create local download dir if it doesn't exist
        dirs = AppDirs("hippunfold", "khanlab")
        download_dir = dirs.user_cache_dir
    return download_dir
