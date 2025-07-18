import os
import re

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
AUTO_DETECT_PATTERN = re.compile(r"&(.+)&")

# smile
SMILE_BG:tuple[int, int, int] = (55, 56, 61)
