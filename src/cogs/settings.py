import discord
from discord.ext import commands
from rdkit import Chem

from constants import SMILE_BG
from db.db import DatabaseHandler


class SettingsCommand(commands.Cog):
    name = "/settings"
    description = "display the current settings for the server."

    def __init__(self, bot):
        self.bot = bot
        self.periodic_table = Chem.GetPeriodicTable()
        self.db_handler:DatabaseHandler = self.bot.db_handler

    @commands.hybrid_command(name="settings", description="display the current settings for the server.")
    async def settings(self, ctx):
        guild_id = str(ctx.guild.id)

        server_config = self.db_handler.get_server_setting(guild_id)
        render_config = self.db_handler.get_render_option(guild_id)
        element_colors = self.db_handler.get_element_colors(guild_id)

        embed = discord.Embed(title="Settings", color=discord.Color.blue())

        # Server settings
        embed.add_field(name="Server Settings", value=" ", inline=False)
        embed.add_field(name="Prefix", value=server_config.get("prefix", "/"))
        embed.add_field(name="Auto Smile", value=bool(server_config.get("auto_smile")))

        # Render options
        bg_color = render_config.get("background_color") or SMILE_BG

        embed.add_field(name="Render Options", value=" ", inline=False)
        embed.add_field(name="Background Color", value=bg_color, inline=False)
        embed.add_field(name="Include Atom Numbers", value=bool(render_config.get("includeAtomNumbers")))
        embed.add_field(name="Add Stereo Annotations", value=bool(render_config.get("addStereoAnnotations")))
        embed.add_field(name="Explicit Methyl", value=bool(render_config.get("explicitMethyl")), inline=False)
        embed.add_field(name="Atom Label Deuterium Tritium", value=bool(render_config.get("atomLabelDeuteriumTritium")))
        embed.add_field(name="dummiesAreAttachments", value=bool(render_config.get("dummiesAreAttachments")), inline=False)

        atom_colors = ""
        for element, color in element_colors.items():
            element_name = self.periodic_table.GetElementName(element)
            atom_colors += f"{element_name}: {color}\n"

        if atom_colors:
            embed.add_field(name="Element Colors", value=atom_colors, inline=False)
        
        await ctx.send(embed=embed)

async def setup(bot):
    await bot.add_cog(SettingsCommand(bot))
