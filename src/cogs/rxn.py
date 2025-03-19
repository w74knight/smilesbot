from discord.ext import commands

class RxnCommand(commands.Cog):
    name = "/rxn"
    description = "Render a reaction equation."

    def __init__(self, bot):
        self.bot = bot

    @commands.hybrid_command(name="rxn", description="Render a reaction equation.")
    async def render(self, ctx, rxn: str):
        palette = self.bot.db_handler.get_element_colors(str(ctx.guild.id))
        await self.bot.smile.render_reaction(ctx, rxn, palette)

async def setup(bot):
    await bot.add_cog(RxnCommand(bot))