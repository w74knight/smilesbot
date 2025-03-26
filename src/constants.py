import re
import os
from dotenv import load_dotenv
from util import smile_rgb

# load sensitive data .env
load_dotenv("env")

# set sensitive data from .env as variables
OWNER_ID = os.getenv("OWNER")
TOKEN = os.getenv("TOKEN")

# set pattern for auto_detect
AUTO_DETECT_PATTERN = re.compile(r"&([^&]+)&")

# smile
SMILE_BG = smile_rgb(44, 45, 49)
