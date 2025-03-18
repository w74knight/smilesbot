import discord
from discord.ext import commands

class SettingsCommand(commands.Cog):
    name = "/settings"
    description = "display the current settings for the server."

    def __init__(self, bot):
        self.bot = bot

    @commands.hybrid_command(name="settings", description="display the current settings for the server.")
    async def settings(self, ctx):
        guild_id = str(ctx.guild.id)

        server_config = self.bot.db_handler.get_server_setting(guild_id)

        if not server_config:
            self.bot.db_handler.create_server_table(guild_id)
            server_config = self.bot.db_handler.get_server_setting(guild_id)

        message = "**Settings**\n"
        print(server_config)
        for key, value in server_config.items():
            message += f"- {key}: {value}\n"
        
        await ctx.send(message)

async def setup(bot):
    await bot.add_cog(SettingsCommand(bot))
