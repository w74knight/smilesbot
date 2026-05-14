from logging import Logger, getLogger

import discord
from discord.ext import commands
from rdkit import Chem

from constants import NAME, SMILE_BG
from db.db import DatabaseHandler


class SettingsCommand(commands.Cog):
    name = "/settings"
    description = "display the current settings for the server."

    def __init__(self, bot):
        self.bot:commands.Bot = bot
        self.periodic_table:Chem.PeriodicTable = Chem.GetPeriodicTable()
        self.db_handler:DatabaseHandler = self.bot.db_handler
        self.logger:Logger = getLogger(NAME)

        self.logger.info("SettingsCommand initialized.")

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

        embed = discord.Embed(title="Render Options", color=discord.Color.blue())
        embed.add_field(name="Background Color", value=render_config.get("background_color", "Default"))
        embed.add_field(name="Color Bonds", value=bool(render_config.get("colorBonds")))
        embed.add_field(name="Add Stereo Annotations", value=bool(render_config.get("addStereoAnnotations")), inline=False)
        embed.add_field(name="Explicit Methyl", value=bool(render_config.get("explicitMethyl")))
        embed.add_field(name="Atom Label Deuterium Tritium", value=bool(render_config.get("atomLabelDeuteriumTritium")), inline=False)
        embed.add_field(name="dummiesAreAttachments", value=bool(render_config.get("dummiesAreAttachments")), inline=False)

        atom_colors_lines = []
        for element, color in element_colors.items():
            element_name = self.periodic_table.GetElementName(element)
            atom_colors_lines.append(f"{element_name}: {color}")

        chunk = ""
        field_count = 0
        for element, color in sorted(element_colors.items(), key=lambda item: int(item[0])):
            element_name = self.periodic_table.GetElementName(int(element))
            line = f"{element_name}: {color}"
            if len(chunk) + len(line) + 1 > 1024:
                embed.add_field(
                    name="Element Colors" if field_count == 0 else "Element Colors (cont.)",
                    value=chunk,
                    inline=False
                )
                chunk = ""
                field_count += 1
            chunk += line + "\n"
        if chunk:
            embed.add_field(
                name="Element Colors" if field_count == 0 else "Element Colors (cont.)",
                value=chunk,
                inline=False
            )
        
        await ctx.send(embed=embed)

async def setup(bot):
    await bot.add_cog(SettingsCommand(bot))
