import discord
from discord.ext import commands
from rdkit import Chem
from constants import SMILE_BG
from smile.pallette import DISCORD_DARK
from util import rgb_to_hex

class GetCommand(commands.Cog):
    name = "/get"
    description = "Get the command prefix for this server."

    def __init__(self, bot):
        self.bot = bot
        self.periodic_table = Chem.GetPeriodicTable()

    @commands.hybrid_group()
    async def get(self, name):
        pass

    @get.command(name="prefix", description="Get the command prefix for this server.")
    async def prefix(self, ctx):
        prefix = self.bot.db_handler.get_server_setting(str(ctx.guild.id)).get("prefix", None)

        if not prefix:
            self.bot.db_handler.create_server_table(str(ctx.guild.id))
            prefix = "/"

        await ctx.send(f"The current prefix is: {prefix}")

    @get.command(name="color", description="Get rgb color for an element when rendering.")
    async def getcolor(self, ctx, element: int):
        if element < 0 or element > 118:
            return await ctx.send("Invalid element number. Must be between 0 and 118.")

        color = self.bot.db_handler.get_element_colors(str(ctx.guild.id)).get(element)
        if not color:
            color = DISCORD_DARK.get(element, (1,1,1))
            color = tuple(int(c * 255) for c in color)

        element_name = self.periodic_table.GetElementName(element)
        hex_color = rgb_to_hex(color)
        embed = discord.Embed(
            title=f"Color for element {element_name}",
            color=discord.Color.from_rgb(*color)
        )
        embed.add_field(name="Element", value=element_name)
        embed.add_field(name="Color", value=hex_color)
        embed.set_thumbnail(url=f"https://dummyimage.com/250/{hex_color}/ffffff&text={element_name}")
        await ctx.send(embed=embed)

    @get.command(name="bgcolor", description="Get the background color for rendering.")
    async def bgcolor(self, ctx):
        bg_color = self.bot.db_handler.get_render_option(str(ctx.guild.id)).get("background_color", None)
        if not bg_color:
            bg_color = tuple(int(c * 255) for c in SMILE_BG)
        else:
            bg_color = tuple(map(int, bg_color.split(',')))

        hex_color = rgb_to_hex(bg_color)
        embed = discord.Embed(
            title="Background Color",
            color=discord.Color.from_rgb(*bg_color)
        )
        embed.add_field(name="Color", value=hex_color)
        embed.set_thumbnail(url=f"https://dummyimage.com/250/{hex_color}/ffffff")
        await ctx.send(embed=embed)

async def setup(bot):
    await bot.add_cog(GetCommand(bot))