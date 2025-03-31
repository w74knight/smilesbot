from logging import Logger, getLogger

import discord
from discord.ext import commands
from rdkit import Chem

from constants import NAME, SMILE_BG
from db.db import DatabaseHandler
from smile.pallette import DISCORD_DARK
from util import admin_only, rgb_to_hex


class GetCommand(commands.Cog):
    name = "/get"
    description = "Get the command prefix for this server."

    def __init__(self, bot):
        self.bot:commands.Bot = bot

        self.db_handler:DatabaseHandler = self.bot.db_handler
        self.logger:Logger = getLogger(NAME)

        self.periodic_table:Chem.PeriodicTable = Chem.GetPeriodicTable()

        self.logger.debug("GetCommand initialized.")

    @commands.hybrid_group()
    async def get(self, name):
        pass
    
    @admin_only()
    @get.command(name="prefix", description="Get the command prefix for this server.")
    async def prefix(self, ctx):
        prefix = self.db_handler.get_server_setting(str(ctx.guild.id)).get("prefix", None)

        await ctx.send(f"The current prefix is: {"/" if not prefix else prefix}")

    @admin_only()
    @get.command(name="color", description="Get rgb color for an element when rendering.")
    async def getcolor(self, ctx, element: int):
        if element < 0 or element > 118:
            return await ctx.send("Invalid element number. Must be between 0 and 118.")

        color = self.db_handler.get_element_colors(str(ctx.guild.id)).get(element)

        if not color: # no color set, display default color
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

    @admin_only()
    @get.command(name="bgcolor", description="Get the background color for rendering.")
    async def bgcolor(self, ctx):
        bg_color = self.db_handler.get_render_option(str(ctx.guild.id)).get("background_color", None)
        if not bg_color:
            bg_color = SMILE_BG
        else:
            bg_color = tuple(int(c) for c in bg_color.split(","))

        hex_color = rgb_to_hex(bg_color)
        embed = discord.Embed(
            title="Background Color",
            color=discord.Color.from_rgb(*bg_color)
        )
        embed.add_field(name="Color", value=hex_color)
        embed.set_thumbnail(url=f"https://dummyimage.com/250/{hex_color}/ffffff")
        await ctx.send(embed=embed)

    @admin_only()
    @get.command(name="render", description="Get the current render options.")
    async def render(self, ctx):
        render_config = self.db_handler.get_render_option(str(ctx.guild.id))

        embed = discord.Embed(title="Render Options", color=discord.Color.blue())
        embed.add_field(name="Background Color", value=render_config.get("background_color", "Default"))
        embed.add_field(name="Color Bonds", value=bool(render_config.get("colorBonds")))
        embed.add_field(name="Include Atom Numbers", value=bool(render_config.get("includeAtomNumbers")), inline=False)
        embed.add_field(name="No Carbon Symbols", value=bool(render_config.get("noCarbonSymbols")))
        embed.add_field(name="Wedge Dashed Bonds", value=bool(render_config.get("wedgeDashedBonds")), inline=False)

        await ctx.send(embed=embed)

async def setup(bot):
    await bot.add_cog(GetCommand(bot))