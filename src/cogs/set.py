import discord
from discord.ext import commands
from rdkit import Chem

def rgb_to_hex(rgb):
    return '{:02x}{:02x}{:02x}'.format(*rgb)

class SetCommand(commands.Cog):
    name = "/set"
    description = "Set server settings."

    def __init__(self, bot):
        self.bot = bot
        self.pt = Chem.GetPeriodicTable()


    @commands.hybrid_group()
    async def set(self, name):
        pass

    @set.command(name="prefix", description="Set a custom prefix for the server.")
    async def prefix(self, ctx, prefix: str):
        if ctx.author.guild_permissions.administrator:
            self.bot.db_handler.set_server_setting(str(ctx.guild.id), prefix)
            self.bot.command_prefix = prefix
        return await ctx.send(f"Prefix set to: {prefix}")
        await ctx.send("You don't have the permissions for this command!")

    @set.command(name="color", description="Set the rgb color for an element when rendering.")
    async def color(self, ctx, element: int, color: str):
        if element < 0 or element > 118:
            await ctx.send("Invalid element number. Must be between 0 and 118.")
            return

        # Convert color from '255 255 255' to '#FFFFFF'
        rgb = tuple(map(int, color.split(',')))
        hex_color = rgb_to_hex(rgb)

        self.bot.db_handler.set_element_color(str(ctx.guild.id), element, color)

        # Create an embed to display the color
        embed = discord.Embed(
            title=f"Color for element {self.pt.GetElementName(element)} set to {hex_color}",
            color=discord.Color.from_rgb(*rgb)
        )
        embed.add_field(name="Element", value=self.pt.GetElementName(element))
        embed.add_field(name="Color", value=hex_color)
        embed.set_thumbnail(url=f"https://dummyimage.com/100x100/{hex_color[1:]}/ffffff")

        await ctx.send(embed=embed)


async def setup(bot):
    await bot.add_cog(SetCommand(bot))
