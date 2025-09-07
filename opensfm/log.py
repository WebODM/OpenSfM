import logging
import os
from typing import Optional
import vmem

def setup() -> None:
    logging.basicConfig(
        format="[%(levelname)s %(asctime)s] %(message)s", level=logging.INFO, force=True, datefmt='%H:%M:%S'
    )


def memory_usage() -> float:
    return vmem.virtual_memory().used / 1024 / 1024 / 1024


def memory_available() -> Optional[int]:
    """Available memory in MB.
    """
    return vmem.virtual_memory().available / 1024 / 1024
