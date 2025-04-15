import os
from discord.ext import commands
from discord import app_commands, Interaction
from logging import getLogger, Logger

from constants import NAME, OWNERS_ID, SUPPORT_GUILD_ID

class SysRestartCommand(commands.Cog):
    name = "/sysrestart"
    description = "Restart the bot"

    def __init__(self, bot):
        self.bot = bot
        self.logger: Logger = getLogger(NAME)
        self.logger.info("SysRestart initialized.")

    def cog_check(self, ctx):
        return ctx.guild.id == SUPPORT_GUILD_ID

    @app_commands.command(name="sysrestart", description="Restart the bot")
    async def restart(self, interaction: Interaction):
        if interaction.user.id not in OWNERS_ID:
            await interaction.response.send_message("🚫 You are not authorized.")
            return

        await interaction.response.send_message("♻️ Restarting bot...")
        await self.bot.close()
        os._exit(1) # exit code 1 indicates a restart


async def setup(bot):
    await bot.add_cog(SysRestartCommand(bot))
