import discord
from discord.ext import commands

class ClearCommand(commands.Cog):
    name = "/clear"
    description = "Clear: bot settings, render colors, all settings."

    def __init__(self, bot):
        self.bot = bot

    @commands.hybrid_group()
    async def clear(self, name):
        pass

    @clear.command(name="settings", description="Clear stored bot settings.")
    async def settings(self, ctx):
        embed = discord.Embed(title="Clear complete!", description="Cleared stored bot settings.")
        self.bot.db_handler.clear_settings(str(ctx.guild.id))
        await ctx.send(embed = embed)

    @clear.command(name="colors", description="Clear stored render color settings.")
    async def colors(self, ctx):
        embed = discord.Embed(title="Clear complete!", description="Cleared stored render color settings.")
        self.bot.db_handler.clear_color(str(ctx.guild.id))
        await ctx.send(embed = embed)

    @clear.command(name="all", description="Clear all stored settings.")
    async def all(self, ctx):
        embed = discord.Embed(title="Clear complete!", description="Cleared all stored bot settings.")
        self.bot.db_handler.clear_color(str(ctx.guild.id))
        self.bot.db_handler.clear_settings(str(ctx.guild.id))
        await ctx.send(embed = embed)

async def setup(bot):
    await bot.add_cog(ClearCommand(bot))
