import discord
from discord.ext import commands

class HelpCommand(commands.Cog):
    name = "/help"
    description = "Show this help menu."

    def __init__(self, bot):
        self.bot = bot

    @commands.hybrid_command(name="help", description="Displays help menu")
    async def help(self, ctx):
        embed = discord.Embed(title="Help Menu", description="List of available commands:")
        
        for module_name in self.bot.extensions:
            module = self.bot.extensions[module_name]
            class_name = ''.join(word.capitalize() for word in module_name.split('.')[-1].split('_')) + "Command"
            command_class = getattr(module, class_name)
            name = getattr(command_class, "name", None)
            description = getattr(command_class, "description", "No description provided.")
            if name:
                embed.add_field(name=f"`{name}`", value=description, inline=False)
        
        await ctx.send(embed=embed)

async def setup(bot):
    await bot.add_cog(HelpCommand(bot))
