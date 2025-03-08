import os
import json

CONFIG = "config.json"
global_server_config = {}

# Load Configurations
def load_configs():
    global global_server_config

    if os.path.exists(CONFIG):
        with open(CONFIG, "r") as file:
            global_server_config = json.load(file)
            return global_server_config
    return global_server_config

# Save Configurations
def save_configs():
    with open(CONFIG, "w") as file:
        json.dump(global_server_config, file, indent=4)

def get_server_config(guild_id):
    global global_server_config

    # Check if the server has settings, else create default ones
    if str(guild_id) not in global_server_config:
        global_server_config[str(guild_id)] = {"prefix": "/", "auto_detect": False}  # Default prefix and auto_detect off
        save_configs()
    return global_server_config[str(guild_id)]

def set_server_config(guild_id, setting, value):
    global global_server_config
    server_config = get_server_config(guild_id)

    if setting == "prefix":
        server_config["prefix"] = value
    elif setting == "auto_detect":
        if value.lower() == "true":
            server_config["auto_detect"] = True
        elif value.lower() == "false":
            server_config["auto_detect"] = False
        else:
            return False

    save_configs()
    return True