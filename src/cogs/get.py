import discord
from discord.ext import commands

def rgb_to_hex(rgb):
    return '{:02x}{:02x}{:02x}'.format(*rgb)

class GetCommand(commands.Cog):
    name = "/get"
    description = "Get the command prefix for this server."

    def __init__(self, bot):
        self.bot = bot

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
    await bot.add_cog(GetCommand(bot))