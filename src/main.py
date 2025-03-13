import discord
from discord.ext import commands
import asyncio
import re
import io
from util import load_configs, save_configs, get_server_config, set_server_config
from rdkit.Chem import AllChem as Chem, rdChemReactions
from rdkit.Chem import Draw
import cirpy
from stdnum import casrn

from dotenv import load_dotenv
import os

load_dotenv()

TOKEN = os.getenv("TOKEN")
print(TOKEN)

pattern = re.compile(r"&([^&]+)&")

# RDKit molecule drawing setup
d2d = Draw.MolDraw2DCairo(350, 300)
opts = d2d.drawOptions()

# Function to validate a SMILES string
def is_valid_smiles(smiles: str):
    try:
        return Chem.MolFromSmiles(smiles) is not None
    except:
        return False

def is_valid_smarts(smarts: str):
    try:
        return Chem.ReactionFromSmarts(smarts, useSMILES = True) is not None
    except:
        return False


# Discord Bot Setup
intents = discord.Intents.default()
intents.message_content = True  # Required for message detection
bot = commands.Bot(command_prefix=lambda bot, message: get_prefix(bot, message), intents=intents, help_command=None)
tree = bot.tree  # Slash command tree

def get_prefix(bot, message):
    guild_id = str(message.guild.id)
    print(get_server_config(guild_id))
    return get_server_config(guild_id).get("prefix", "/")

@bot.hybrid_command(name="help", description="Displays help menu")
async def help(ctx):
    embed = discord.Embed(title="Help Menu", description="List of available commands:", color=0x2f3136)
    embed.add_field(name="`/help`", value="Show this help menu.", inline=False)
    embed.add_field(name="`/smileshelp`", value="Quick guide to SMILES.", inline=False)
    embed.add_field(name="`/render <render_string>`", value="Generate a molecular structure image from a molecume name, ID, or SMILES string.", inline=False)
    embed.add_field(name="`/rxn <rxn_string>`", value=":construction: UNDER CONSTRUCTION :construction:", inline=False)
    embed.add_field(name="`/auto_detect <True|False>`", value=":construction: UNDER CONSTRUCTION :construction:", inline=False)
    await ctx.send(embed=embed)

@bot.hybrid_command(name="smileshelp", description="Quick guide to SMILES.")
async def smileshelp(ctx):
    embed = discord.Embed(title="SMILES Syntax", description="Quick guide to SMILES", color=0x2f3136)
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

@bot.hybrid_command(name="setprefix", description="Set a custom prefix for the server.")
async def setprefix(ctx, new_prefix: str):
    guild_id = str(ctx.guild.id)
    
    # if not interaction.user.guild_permissions.administrator:
    #     await interaction.response.send_message("You need to be an administrator to use this command!", ephemeral=False)
    #     return

    if set_server_config(guild_id, "prefix", new_prefix):
        await ctx.send(f"The command prefix has been updated to `{new_prefix}`.", ephemeral=False)
    else:
        await ctx.send("An error occurred while setting the prefix.", ephemeral=True)


@bot.hybrid_command(name="settings", description="Settings")
async def settings(ctx: discord.Interaction):
    guild_id = str(ctx.guild.id)

    # Fetch server config for the given guild
    server_config = get_server_config(guild_id)

    message = "**Settings**\n"
    for key, value in server_config.items():
        message += f"- {key}: {value}\n"
    
    # Send the settings message to the server
    await ctx.send(message)

@bot.hybrid_command(name="setwidth", description="Set a custom prefix for the server.")
async def setwidth(ctx, new_width: str):
    guild_id = str(ctx.guild.id)

    # if not interaction.user.guild_permissions.administrator:
    #     await interaction.response.send_message("You need to be an administrator to use this command!", ephemeral=False)
    #     return

    if set_server_config(guild_id, "width", new_width):
        await ctx.send(f"The render width has been changed to `{new_width}`.", ephemeral=False)
    else:
        await ctx.send("An error occurred while setting the prefix.", ephemeral=True)

# Function to create a molecule image
def create_molecule_image(mol):
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except:
        print("Kekulization failed, skipping.")
    opts.bondLineWidth = setprefix
    opts.updateAtomPalette({6: (1, 1, 1)})  # Changed Carbon to White
    opts.updateAtomPalette({6: (1, 1, 1)})
    opts.setBackgroundColour((44 / 255, 45 / 255, 49 / 255))  # Sets render background to Discord Dark theme
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    bio = io.BytesIO(d2d.GetDrawingText())
    bio.seek(0)
    return bio

@bot.hybrid_command(name="render", description="Render a molecule.")
async def render(ctx, molecule_ID: str):
    if is_valid_smiles(molecule_ID):
        mol = Chem.MolFromSmiles(molecule_ID)
        loop = asyncio.get_running_loop()
        img = await loop.run_in_executor(None, create_molecule_image, mol)

        embed = discord.Embed(title=f"`{molecule_ID}`", color=0x2f3136)
        embed.set_image(url="attachment://molecule.png")
        await ctx.send(embed=embed, file=discord.File(img, filename="molecule.png"))
        d2d.ClearDrawing()
    else:
        try:
            resolve = cirpy.resolve(molecule_ID, 'smiles')
            mol_resolve = Chem.MolFromSmiles(resolve)
            loop = asyncio.get_running_loop()
            img = await loop.run_in_executor(None, create_molecule_image, mol_resolve)

            embed = discord.Embed(title=f"`{molecule_ID}`", color=0x2f3136)
            embed.set_image(url="attachment://molecule.png")
            await ctx.send(embed=embed, file=discord.File(img, filename="molecule.png"))
            d2d.ClearDrawing()
        except:
            await ctx.send(f"{molecule_ID} is invalid, please try with a different compound ID or check for typos/erros!")

@bot.hybrid_command(name="auto_detect", description="Enable or disable automatic smile[...] message detection.")
async def auto_detect(ctx: discord.Interaction, option: str):
    global server_configs
    if not ctx.user.guild_permissions.administrator:
        await ctx.response.send_message("You need to be an administrator to use this command!", ephemeral=False)
        return
    
    guild_id = str(ctx.guild.id)

    if set_server_config(guild_id, "auto_detect", option):
        await ctx.send(f"Auto-detection has been set to `{option}`.", ephemeral=False)
    else:
        await ctx.send("Invalid option. Please use `true` or `false` for auto_detect.", ephemeral=True)

@bot.event
async def on_message(message):
    if message.author.bot:
        return  # Ignore bot messages
    guild_id = str(message.guild.id)
    server_config = get_server_config(guild_id)
    # Check if auto-detection is enabled for the server
    if server_config.get("auto_detect", False):
        match = pattern.search(message.content)
        print(match)
        if match:
            smiles_str = match.group(1)
            if is_valid_smiles(smiles_str):
                mol = Chem.MolFromSmiles(smiles_str)
                loop = asyncio.get_running_loop()
                img = await loop.run_in_executor(None, create_molecule_image, mol)
                embed = discord.Embed(title=f"SMILES: `{smiles_str}`", color=0x2f3136)
                embed.set_image(url="attachment://molecule.png")
                await message.channel.send(embed=embed, file=discord.File(img, filename="molecule.png"))
                d2d.ClearDrawing()
            else:
                try:
                    resolve = cirpy.resolve(smiles_str, 'smiles', ['cas_number', 'name_by_cir', 'name_by_opsin'])
                    mol_resolve = Chem.MolFromSmiles(resolve)
                    loop = asyncio.get_running_loop()
                    img = await loop.run_in_executor(None, create_molecule_image, mol_resolve)
                    embed = discord.Embed(title=f"SMILES: `{smiles_str}`", color=0x2f3136)
                    embed.set_image(url="attachment://molecule.png")
                    await message.channel.send(embed=embed, file=discord.File(img, filename="molecule.png"))
                    d2d.ClearDrawing()
                except:
                    await message.channel.send(f"Error rendering {smiles_str}")
    await bot.process_commands(message)
### **Bot Ready Event** ###
#set presence as "watching /smiles"
@bot.event
async def on_ready():
    await tree.sync()
    await bot.change_presence(status=discord.Status.online, activity=discord.Activity(type=discord.ActivityType.watching, name="/smiles"))
# Run Bot
bot.run(TOKEN)