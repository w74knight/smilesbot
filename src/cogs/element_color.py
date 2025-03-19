import discord
from discord.ext import commands
from .util import admin_only
from rdkit import Chem

def rgb_to_hex(rgb):
    return '{:02x}{:02x}{:02x}'.format(*rgb)

class ElementColorCommand(commands.Cog):
    name = "/setcolor /getcolor"
    description = "Set and get the rgb color for an element when rendering."

    def __init__(self, bot):
        self.bot = bot
        self.pt = Chem.GetPeriodicTable()

    @admin_only()
    @commands.hybrid_command(name="setcolor", description="Set the rgb color for an element when rendering.")
    async def setcolor(self, ctx, element: int, color: str):
        if ctx.author.guild_permissions.administrator:
            if element < 0 or element > 118:
                await ctx.send("Invalid element number. Must be between 0 and 118.")
                return
            # Convert color from '255 255 255' to '#FFFFFF'
            rgb = tuple(map(int, color.split(', ')))
            hex_color = rgb_to_hex(rgb)

            self.bot.db_handler.set_element_color(str(ctx.guild.id), element, color)

            # Create an embed to display the color
            embed = discord.Embed(
                title=f"Color for element {self.pt.GetElementName(element)} set to #{hex_color}",
                color=discord.Color.from_rgb(*rgb)
            )
            embed.add_field(name="Element", value=self.pt.GetElementName(element))
            embed.add_field(name="Color", value=hex_color)
            embed.set_thumbnail(url=f"https://dummyimage.com/100x100/{hex_color}/ffffff")

            return await ctx.send(embed=embed)
        await ctx.send("You don't have the permissions for this command!")

    @commands.hybrid_command(name="getcolor", description="Get rgb color for an element when rendering.")
    async def getcolor(self, ctx, element: int):
        color = self.bot.db_handler.get_element_colors(str(ctx.guild.id)).get(element)
        element_name = self.pt.GetElementName(element)
        if color:
            color_str = ' '.join(map(str, color))
            hex_color = rgb_to_hex(color)
            embed = discord.Embed(
                title=f"Color for element {element_name}",
                color=discord.Color.from_rgb(*color)
            )
            embed.add_field(name="Element", value=element_name)
            embed.add_field(name="Color", value=hex_color)
            embed.set_thumbnail(url=f"https://dummyimage.com/{element_name}/{hex_color[1:]}/ffffff")
            await ctx.send(embed=embed)
        else:
            await ctx.send(f"No color set for element {element_name}")

async def setup(bot):
    await bot.add_cog(ElementColorCommand(bot))