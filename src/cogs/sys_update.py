from discord import File
from discord.ext import commands
from discord import app_commands, Interaction
from logging import getLogger, Logger
import subprocess
import os
import io

from constants import NAME, OWNERS_ID, SUPPORT_GUILD_ID

PROJECT_DIR = "/home/smile.bot/smilesbot"
PYTHON_PATH = "/home/smile.bot/smilesbot/env/bin/python3.13"

class SysUpdateCommand(commands.Cog):
    name = "/sysupdate"
    description = "Update the bot"

    def __init__(self, bot):
        self.bot = bot
        self.logger: Logger = getLogger(NAME)
        self.logger.info("SysUpdate initialized.")

    def cog_check(self, ctx):
        return ctx.guild.id == SUPPORT_GUILD_ID
    
    @app_commands.command(name="sysupdate", description="pull the latest changes from git repo")
    async def update(self, interaction: Interaction):
        if interaction.user.id not in OWNERS_ID:
            await interaction.response.send_message("🚫 You are not authorized.")
            return

        await interaction.response.send_message("Updating from Git...")

        try:
            git_process = subprocess.run(["git", "pull"], cwd=PROJECT_DIR, check=True, capture_output=True, text=True)
            git_output = git_process.stdout + git_process.stderr

            if "Already up to date." in git_output:
                await interaction.followup.send("Already up to date.")
                return

            pip_process = subprocess.run(
                [PYTHON_PATH, "-m", "pip", "install", "--upgrade", "-r", "requirements.txt"],
                cwd=PROJECT_DIR,
                check=True,
                capture_output=True,
                text=True
            )
            pip_output = pip_process.stdout + pip_process.stderr

            full_log = (
                f"**Git Output:**\n```\n{git_output}\n```\n"
                f"**Pip Output:**\n```\n{pip_output}\n```\n"
                f"**Restarting bot...**"
            )
            log_file = io.BytesIO(full_log.encode("utf-8"))
            log_file.seek(0)
            log_file.name = "log.txt"

            await interaction.followup.send(content="✅ Update complete. Here's the full log:", file=File(log_file))
            await self.bot.close()
            os._exit(1)  # exit code 1 indicates a restart
        except subprocess.CalledProcessError as e:
            self.logger.error(
                f"An error occurred while updating:\n```\n{e.stderr}\n```\n"
                f"Please check the logs for more details."
            )
        except Exception as e:
            self.logger.error(
                f"An unexpected error occurred:\n```\n{str(e)}\n```\n"
                f"Please check the logs for more details."
            )
        finally:
            self.logger.info("Update process completed.")
            self.logger.info(f"Git Output: {git_output}")


async def setup(bot):
    await bot.add_cog(SysUpdateCommand(bot))
