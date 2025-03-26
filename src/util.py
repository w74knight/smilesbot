def get_prefix(bot, message):
    if not message.guild:
        return "/"
    settings = bot.db_handler.get_server_setting(str(message.guild.id))
    prefix = settings['prefix'] if settings else "/"
    return prefix

def smile_rgb(r: int, g: int, b: int):
    return (r/255, g/255, b/255)
