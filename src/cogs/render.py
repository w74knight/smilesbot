from discord.ext import commands

class RenderCommand(commands.Cog):
    name = "/render"
    description = "Render a molecule."

    def __init__(self, bot):
        self.bot = bot

    @commands.hybrid_command(name="render", description="Render a molecule.")
    async def render(self, ctx, mlcl: str):
        pallette = self.bot.db_handler.get_element_colors(str(ctx.guild.id))
        await self.bot.smile.render_molecule(ctx, mlcl, pallette)

async def setup(bot):
    await bot.add_cog(RenderCommand(bot))
