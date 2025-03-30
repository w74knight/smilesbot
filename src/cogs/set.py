import discord
from discord.ext import commands
from rdkit import Chem
from util import rgb_to_hex
from db.db import DatabaseHandler
from util import admin_only

class SetCommand(commands.Cog):
    name = "/set"
    description = "Set server settings."

    def __init__(self, bot):
        self.bot = bot
        self.pt = Chem.GetPeriodicTable()
        self.db_handler:DatabaseHandler = self.bot.db_handler


    @commands.hybrid_group()
    async def set(self, name):
        pass
    
    @admin_only()
    @set.command(name="prefix", description="Set a custom prefix for the server.")
    async def prefix(self, ctx, prefix: str):
        self.db_handler.server_settings.set_server_setting(str(ctx.guild.id), "prefix", prefix)
        self.bot.command_prefix = prefix
        return await ctx.send(f"Prefix set to: {prefix}")

    @admin_only()
    @set.command(name="color", description="Set the rgb color for an element when rendering.")
    async def color(self, ctx, element: int, color: str):
        if element < 0 or element > 118:
            return await ctx.send("Invalid element number. Must be between 0 and 118.")

        # Convert color from '255 255 255' to '#FFFFFF'
        rgb = tuple(map(int, color.split(',')))
        hex_color = rgb_to_hex(rgb)

        self.db_handler.element_colors.set_element_color(str(ctx.guild.id), element, color)

        # Create an embed to display the color
        embed = discord.Embed(
            title=f"Color for element {self.pt.GetElementName(element)} set to {hex_color}",
            color=discord.Color.from_rgb(*rgb)
        )
        embed.add_field(name="Element", value=self.pt.GetElementName(element))
        embed.add_field(name="Color", value=hex_color)
        embed.set_thumbnail(url=f"https://dummyimage.com/250/{hex_color}/ffffff")

        await ctx.send(embed=embed)

    @admin_only()
    @set.command(name="bgcolor", description="Set the background color for rendering.")
    async def bgcolor(self, ctx, color: str):
        # Convert color from '255 255 255' to '#FFFFFF'
        rgb = tuple(map(int, color.split(',')))
        hex_color = rgb_to_hex(rgb)

        self.db_handler.render_options.set_bgcolor(str(ctx.guild.id), color)

        # Create an embed to display the color
        embed = discord.Embed(
            title=f"Background color set to {hex_color}",
            color=discord.Color.from_rgb(*rgb)
        )
        embed.add_field(name="Background Color", value=hex_color)
        print(hex_color)
        embed.set_thumbnail(url=f"https://dummyimage.com/250/{hex_color}/ffffff&text=Background")

        await ctx.send(embed=embed)

    @admin_only()
    @set.command(name="colorbonds", description="Set whether to color bonds in the render.")
    async def colorbonds(self, ctx, option: bool):
        self.db_handler.render_options.set(str(ctx.guild.id), "colorBonds", option)
        await ctx.send(f"Color bonds set to: {option}")

    @admin_only()
    @set.command(name="atomnumbers", description="Set whether to include atom numbers in the render.")
    async def atomnumbers(self, ctx, option: bool):
        self.db_handler.render_options.set(str(ctx.guild.id), "includeAtomNumbers", option)
        await ctx.send(f"Include atom numbers set to: {option}")

    @admin_only()
    @set.command(name="nocarbon", description="Set whether to show carbon symbols in the render.")
    async def nocarbon(self, ctx, option: bool):
        self.db_handler.render_options.set(str(ctx.guild.id), "noCarbonSymbols", option)
        await ctx.send(f"No carbon symbols set to: {option}")

    @admin_only()
    @set.command(name="wedgedashed", description="Set whether to use wedge dashed bonds in the render.")
    async def wedgedashed(self, ctx, option: bool):
        self.db_handler.render_options.set(str(ctx.guild.id), "wedgeDashedBonds", option)
        await ctx.send(f"Wedge dashed bonds set to: {option}")


async def setup(bot):
    await bot.add_cog(SetCommand(bot))
