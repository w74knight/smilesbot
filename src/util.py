from discord.ext import commands

def get_prefix(bot, message):
    if not message.guild:
        return "/"
    settings = bot.db_handler.get_server_setting(str(message.guild.id))
    prefix = settings['prefix'] if settings else "/"
    return prefix

def smile_rgb(r: int, g: int, b: int):
    return (r/255, g/255, b/255)

def rgb_to_hex(rgb):
    return '{:02x}{:02x}{:02x}'.format(*rgb)

def admin_only():
    async def predicate(ctx):
        if not ctx.author.guild_permissions.administrator:
            await ctx.send("You do not have permission to use this command.")
        return ctx.author.guild_permissions.administrator
    return commands.check(predicate)

def transform_rgb_to_smile(atoms: dict):
    return {k: smile_rgb(*v) for k, v in atoms.items()}

def complement_color(rgb):
    # took from chatgpt
    return tuple(1.0 - (c / 255.0) for c in rgb)