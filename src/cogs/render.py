from logging import Logger, getLogger
from typing import override

from discord.ext import commands
from constants import NAME

def not_bot():
    async def predicate(ctx):
        return ctx.author.id != ctx.bot.user.id
    return commands.check(predicate)

class RenderCommand(commands.Cog):
    name = "/render"
    description = "Render a molecule."

    def __init__(self, bot):
        self.bot:commands.Bot = bot
        self.logger:Logger = getLogger(NAME)

        self.logger.info("RenderCommand initialized.")

    @commands.hybrid_group()
    async def render(self, name):
        pass
    
    @render.command(name="mlcl", description="Render a molecule.")
    @not_bot()
    async def mlcl(self, ctx, mlcl: str, legends:str = "", highlightatoms: str = ""):
        highlightatoms = tuple(map(int, highlightatoms.split(","))) if highlightatoms else ()

        await self.bot.smile.render_molecule(ctx, mlcl, str(ctx.guild.id), legends=legends, highlightAtoms=highlightatoms)

    @render.command(name="rxn", description="Render a reaction equation.")
    async def rxn(self, ctx, rxn: str):
        await self.bot.smile.render_reaction(ctx, rxn, str(ctx.guild.id))

async def setup(bot):
    await bot.add_cog(RenderCommand(bot))

