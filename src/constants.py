import os
import re
from pathlib import Path
import json

from dotenv import load_dotenv

__version__ = "1.0.1"
NAME = "SmileBot"

# load sensitive data .env
load_dotenv()

# set sensitive data from .env as variables
SUPPORT_GUILD_ID:int = int(os.getenv("SUPPORT_GUILD_ID"))
OWNERS_ID:list = [int(user) for user in os.getenv("OWNERS").split(",")]
TOKEN:str = os.getenv("TOKEN")
WEBHOOK_URL:str = os.getenv("WEBHOOK_URL")

# set pattern for auto_detect
AUTO_DETECT_PATTERN = re.compile(r"^&[^&]+&$")

# smiles
SMILE_BG:tuple[int, int, int] = (55, 56, 61)

base_dir = Path(__file__).resolve().parent
DISCORD_DARK_PATH = base_dir / "smiles" / "palette.json"
with DISCORD_DARK_PATH.open("r", encoding="utf-8") as f:
    _raw = json.load(f)["discord dark"]

DISCORD_DARK_JSON = {int(k): tuple(v) for k, v in _raw.items()}
