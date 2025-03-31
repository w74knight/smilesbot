from logging import Logger, getLogger

import discord
from discord.ext import commands

from db.db import DatabaseHandler
from constants import NAME


class ClearCommand(commands.Cog):
    name = "/clear"
    description = "Clear: bot settings, render colors, all settings."

    def __init__(self, bot):
        self.bot:commands.Bot = bot
        self.database_handler:DatabaseHandler = self.bot.db_handler
        self.logger:Logger = getLogger(NAME)

        self.logger.info("ClearCommand initialized.")

    @commands.hybrid_group()
    async def clear(self, name):
        pass

    @clear.command(name="settings", description="Clear stored bot settings.")
    async def settings(self, ctx):
        embed = discord.Embed(title="Clear complete!", description="Cleared stored bot settings.")
        self.database_handler.server_settings.clear(str(ctx.guild.id))
        await ctx.send(embed = embed)

    @clear.command(name="render", description="Clear stored render settings.")
    async def render(self, ctx):
        embed = discord.Embed(title="Clear complete!", description="Cleared stored render settings.")
        self.database_handler.render_options.clear(str(ctx.guild.id))
        await ctx.send(embed = embed)

    @clear.command(name="colors", description="Clear stored render color settings.")
    async def colors(self, ctx):
        embed = discord.Embed(title="Clear complete!", description="Cleared stored render color settings.")
        self.database_handler.element_colors.clear(str(ctx.guild.id))
        await ctx.send(embed = embed)

    @clear.command(name="all", description="Clear all stored settings.")
    async def all(self, ctx):
        embed = discord.Embed(title="Clear complete!", description="Cleared all stored bot settings.")

        # this will clear render colors and server settings
        self.database_handler.clear(str(ctx.guild.id))
        await ctx.send(embed = embed)

async def setup(bot):
    await bot.add_cog(ClearCommand(bot))
