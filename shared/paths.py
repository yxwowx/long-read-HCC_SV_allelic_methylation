"""shared/paths.py ==========
Repo-root path configuration, loaded from .env (see .env.example for the template).
Import from pipeline scripts as:
    from shared.paths import HCC_DATA_DIR, REFERENCE_DIR, PATIENT_CODE_MAP
"""
import os
from pathlib import Path

from dotenv import load_dotenv

REPO_ROOT = Path(__file__).resolve().parent.parent
load_dotenv(REPO_ROOT / ".env")

HCC_DATA_DIR = Path(os.environ["HCC_DATA_DIR"])
REFERENCE_DIR = Path(os.environ["REFERENCE_DIR"])
PATIENT_CODE_MAP = Path(os.environ["PATIENT_CODE_MAP"])
