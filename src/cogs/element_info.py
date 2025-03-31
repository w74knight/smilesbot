import discord
from discord.ext import commands
from rdkit import Chem


class ElementInfoCommand(commands.Cog):
    name = "/info"
    description = "Get information about an element."

    def __init__(self, bot):
        self.bot = bot
        self.periodic_table = Chem.GetPeriodicTable()

    @commands.command(name='info')
    async def info(self, ctx, element: int):
        try:
            symbol = self.periodic_table.GetElementSymbol(element)
            name = self.periodic_table.GetElementName(element)
            atomic_weight = self.periodic_table.GetAtomicWeight(element)
            covalent_radius = self.periodic_table.GetRcovalent(element)
            vdw_radius = self.periodic_table.GetRvdw(element)

            embed = discord.Embed(title=f"Element: {name} ({symbol})", color=0x00ff00)
            embed.add_field(name="Atomic Number", value=element)
            embed.add_field(name="Symbol", value=symbol)
            embed.add_field(name="Name", value=name)
            embed.add_field(name="Atomic Weight", value=atomic_weight)
            embed.add_field(name="Covalent Radius", value=covalent_radius)
            embed.add_field(name="Van der Waals Radius", value=vdw_radius)

            await ctx.send(embed=embed)
        except:
            await ctx.send(f"Element '{element}' not found. Please provide a valid element symbol.")

async def setup(bot):
    await bot.add_cog(ElementInfoCommand(bot))