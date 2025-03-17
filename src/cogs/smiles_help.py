import discord
from discord.ext import commands

class SmilesHelpCommand(commands.Cog):
    name = "/smileshelp"
    description = "Quick guide to SMILES."

    def __init__(self, bot):
        self.bot = bot

    @commands.hybrid_command(name="smileshelp", description="Quick guide to SMILES.")
    async def smileshelp(self, ctx):
        embed = discord.Embed(title="SMILES Syntax", description="Quick guide to SMILES")
        embed.add_field(name="Bonds",
                        value="` . ` Disconnected structures such as ionic bonds and multiple compounds in an image.", inline=False)
        embed.add_field(name="",
                        value="` - ` Single bonds (optional/usually omitted).", inline=False)
        embed.add_field(name="",
                        value="` = ` Double bonds.", inline=False)
        embed.add_field(name="",
                        value="` # ` Tripple bonds.", inline=False)
        embed.add_field(name="",
                        value="` $ ` Quadruple bonds.", inline=False)
        embed.add_field(name="",
                        value="` : ` Aromatic 'one and a half' bonds.", inline=False)
        embed.add_field(name="",
                        value="` / ` Single bond (directional) use for cis-trans stereochemistry.", inline=False)
        embed.add_field(name="",
                        value="` \ ` Single bond (directional) use for cis-trans stereochemistry.", inline=False)
        embed.add_field(name="Branching",
                        value="`( )` are used to denote branches, for example Isopropyl Alcohol is `CC(O)C`. Bond type is put __inside__ the parentheses like `CCC(=O)O`.", inline=False)
        embed.add_field(name="Stereochemistry",
                        value="` / and \` are used for cis-trans stereochemistry, for example (E)-1,2-Dichloroethene is `Cl/C=C/Cl`.", inline=False)
        embed.add_field(name="",
                        value="` @ and @@ ` are used to denote S (@) and R (@@) stereocenters.", inline=False)
        embed.add_field(name="Isotopes and charges",
                        value="`[]` are used for charges and isotopes, for example Carbon-14 is `[14C]` and Sodium Chloride is `[Na+].[Cl-]`")
        await ctx.send(embed=embed)

async def setup(bot):
    await bot.add_cog(SmilesHelpCommand(bot))
