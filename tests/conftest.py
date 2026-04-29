import os
import sys
from pathlib import Path


os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

project_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(project_root))

mpl_config_dir = Path(
    os.environ.get("MPLCONFIGDIR")
    or project_root / ".pytest-cache" / "matplotlib"
)
mpl_config_dir.mkdir(parents=True, exist_ok=True)
os.environ["MPLCONFIGDIR"] = str(mpl_config_dir)
