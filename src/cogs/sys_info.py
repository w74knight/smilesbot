import platform
import psutil
import time
from datetime import timedelta
from logging import getLogger, Logger

import discord
from discord.ext import commands
from discord import app_commands, Interaction

from constants import NAME, OWNERS_ID, SUPPORT_GUILD_ID

try:
    import GPUtil
    GPU_ENABLED = True
except ImportError:
    GPU_ENABLED = False

class SysInfoCommand(commands.Cog):
    name = "/sysinfo"
    description = "Display detailed system information."

    def __init__(self, bot):
        self.bot: commands.Bot = bot
        self.logger: Logger = getLogger(NAME)
        self.logger.info("SystemInfoCommand initialized.")

    def cog_check(self, ctx):
        return ctx.guild.id == SUPPORT_GUILD_ID

    def user_allowed(self, interaction: Interaction) -> bool:
        return (
            interaction.guild
            and interaction.guild.id == SUPPORT_GUILD_ID
            and interaction.user.id in OWNERS_ID
        )

    @app_commands.command(name="sysinfo", description="Show detailed system info (restricted)")
    async def sysinfo(self, interaction: Interaction):
        if not self.user_allowed(interaction):
            await interaction.response.send_message("🚫 You are not authorized to use this command.")
            return

        uname = platform.uname()
        cpu_percent = psutil.cpu_percent(interval=1)
        virtual_mem = psutil.virtual_memory()
        total_mem = virtual_mem.total / (1024 ** 3)
        used_mem = virtual_mem.used / (1024 ** 3)
        disk = psutil.disk_usage('/')
        total_disk = disk.total / (1024 ** 3)
        used_disk = disk.used / (1024 ** 3)
        uptime_seconds = time.time() - psutil.boot_time()
        uptime_str = str(timedelta(seconds=int(uptime_seconds)))

        msg = (
            f"**🖥️ System Info**\n"
            f"**OS:** {uname.system} {uname.release} ({uname.version})\n"
            f"**Machine:** {uname.machine}\n"
            f"**Temperature**: {psutil.sensors_temperatures()['cpu_thermal'][0].current}°C\n"
            f"**Processor:** {uname.processor}\n"
            f"**CPU Usage:** {cpu_percent}%\n"
            f"**Memory Usage:** {used_mem:.2f}GB / {total_mem:.2f}GB\n"
            f"**Disk Usage:** {used_disk:.2f}GB / {total_disk:.2f}GB\n"
            f"**Uptime:** {uptime_str}\n"
        )

        if GPU_ENABLED:
            gpus = GPUtil.getGPUs()
            if gpus:
                for i, gpu in enumerate(gpus):
                    msg += (
                        f"\n**🎮 GPU {i}**\n"
                        f"Name: {gpu.name}\n"
                        f"Load: {gpu.load * 100:.1f}%\n"
                        f"Memory Usage: {gpu.memoryUsed}MB / {gpu.memoryTotal}MB\n"
                        f"Temperature: {gpu.temperature}°C\n"
                    )
            else:
                msg += "\n**🎮 GPU:** No GPU found.\n"

        await interaction.response.send_message(msg)


async def setup(bot):
    await bot.add_cog(SysInfoCommand(bot))
