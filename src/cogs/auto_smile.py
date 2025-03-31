from discord.ext import commands

from db.db import DatabaseHandler
from util import admin_only


class AutoSmileCommand(commands.Cog):
    name = "/auto_smile"
    description = "Enable or disable automatic smile[...] message detection."
    setting_key = 'auto_smile'

    def __init__(self, bot):
        self.bot = bot
        self.database_handler:DatabaseHandler = self.bot.db_handler

    @admin_only()
    @commands.hybrid_command(name="auto_smile", description="Enable or disable automatic smile[...] message detection.")
    async def auto_smile(self, ctx):
        guild_id = str(ctx.guild.id)

        server_config = self.database_handler.get_server_setting(guild_id)
        auto_smile = server_config.get(self.setting_key, False)

        self.database_handler.server_settings.set_server_setting(guild_id, self.setting_key, not auto_smile)

        await ctx.send(f"Automatic smile detection is now {'enabled' if not auto_smile else 'disabled'}.")

async def setup(bot):
    await bot.add_cog(AutoSmileCommand(bot))
