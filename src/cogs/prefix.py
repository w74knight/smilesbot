import discord
from discord.ext import commands

class PrefixCommand(commands.Cog):
    name = "/setprefix /getprefix"
    description = "Set and get the command prefix for this server."

    def __init__(self, bot):
        self.bot = bot

    @commands.hybrid_command(name="setprefix", description="Set a custom prefix for the server.")
    async def setprefix(self, ctx, prefix: str):
        self.bot.db_handler.set_server_setting(str(ctx.guild.id), prefix)

        self.bot.command_prefix = prefix
        await ctx.send(f"Prefix set to: {prefix}")

    @commands.hybrid_command(name="getprefix", description="Get the command prefix for this server.")
    async def getprefix(self, ctx):
        prefix = self.bot.db_handler.get_server_setting(str(ctx.guild.id)).get("prefix", None)

        if not prefix:
            self.bot.db_handler.create_server_table(str(ctx.guild.id))
            prefix = "/"

        await ctx.send(f"The current prefix is: {prefix}")

async def setup(bot):
    await bot.add_cog(PrefixCommand(bot))