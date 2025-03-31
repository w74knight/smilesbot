import os
import re

from dotenv import load_dotenv

__version__ = "1.0.0"
NAME = "SmileBot"

# load sensitive data .env

# set sensitive data from .env as variables
OWNER_ID = os.getenv("OWNER")
TOKEN = os.getenv("TOKEN")
WEBHOOK_URL = os.getenv("WEBHOOK_URL")

# set pattern for auto_detect
AUTO_DETECT_PATTERN = re.compile(r"&(.+)&")

# smile
SMILE_BG = (44, 45, 49)
