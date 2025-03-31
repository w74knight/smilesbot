import discord
from discord.ext import commands
from rdkit import Chem

from db.db import DatabaseHandler
from util import admin_only, rgb_to_hex


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

    def __boolOptionHelper(self, ctx, option_name):
        prev_value = self.db_handler.render_options.get(str(ctx.guild.id), option_name)
        self.db_handler.render_options.set(str(ctx.guild.id), option_name, not prev_value)

        return not prev_value

    @admin_only()
    @set.command(name="include_atom_numbers", description="Set whether to include atom numbers in the rendering.")
    async def include_atom_numbers(self, ctx):
        value = self.__boolOptionHelper(ctx, "includeAtomNumbers")
        await ctx.send(f"Include atom numbers set to: {value}")

    @admin_only()
    @set.command(name="add_stereo_annotations", description="Set whether to add stereo annotations in the rendering.")
    async def add_stereo_annotations(self, ctx):
        value = self.__boolOptionHelper(ctx, "addStereoAnnotations")
        await ctx.send(f"Add stereo annotations set to: {value}")

    @admin_only()
    @set.command(name="explicit_methyl", description="Set whether to include explicit methyl groups in the rendering.")
    async def explicit_methyl(self, ctx):
        value = self.__boolOptionHelper(ctx, "explicitMethyl")
        await ctx.send(f"Explicit methyl groups set to: {value}")

    @admin_only()
    @set.command(name="atom_label_deuterium_tritium", description="Set whether to label deuterium and tritium atoms.")
    async def atom_label_deuterium_tritium(self, ctx):
        value = self.__boolOptionHelper(ctx, "atomLabelDeuteriumTritium")
        await ctx.send(f"Label deuterium and tritium atoms set to: {value}")

    @admin_only()
    @set.command(name="dummies_are_attachments", description="Set whether dummies are attachments.")
    async def dummies_are_attachments(self, ctx):
        value = self.__boolOptionHelper(ctx, "dummiesAreAttachments")
        await ctx.send(f"Dummies are attachments set to: {value}")

async def setup(bot):
    await bot.add_cog(SetCommand(bot))
