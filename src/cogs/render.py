from discord.ext import commands

class RenderCommand(commands.Cog):
    name = "/render"
    description = "Render a molecule."

    def __init__(self, bot):
        self.bot = bot

    @commands.hybrid_group()
    async def render(self, name):
        pass

    @render.command(name="mlcl", description="Render a molecule.")
    async def mlcl(self, ctx, mlcl: str):
        palette = self.bot.db_handler.get_element_colors(str(ctx.guild.id))
        await self.bot.smile.render_molecule(ctx, mlcl, palette)

    @render.command(name="rxn", description="Render a reaction equation.")
    async def rxn(self, ctx, rxn: str):
        palette = self.bot.db_handler.get_element_colors(str(ctx.guild.id))
        await self.bot.smile.render_reaction(ctx, rxn, palette)

async def setup(bot):
    await bot.add_cog(RenderCommand(bot))

