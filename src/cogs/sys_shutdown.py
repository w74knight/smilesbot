import os
from discord.ext import commands
from discord import app_commands, Interaction
from logging import getLogger, Logger

from constants import NAME, OWNERS_ID, SUPPORT_GUILD_ID

class SysShutdownCommand(commands.Cog):
    name = "/sysshutdown"
    description = "Shutdown the bot"

    def __init__(self, bot):
        self.bot = bot
        self.logger: Logger = getLogger(NAME)
        self.logger.info("SysShutdown initialized.")

    def cog_check(self, ctx):
        return ctx.guild.id == SUPPORT_GUILD_ID

    @app_commands.command(name="sysshutdown", description="Shutdown the bot")
    async def shutdown(self, interaction: Interaction):
        if interaction.user.id not in OWNERS_ID:
            await interaction.response.send_message("🚫 You are not authorized.")
            return

        await interaction.response.send_message("Shutting down...")
        await self.bot.close()
        os._exit(0)  # exit gracefully


async def setup(bot):
    await bot.add_cog(SysShutdownCommand(bot))
