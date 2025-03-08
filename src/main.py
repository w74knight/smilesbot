import discord
from discord.ext import commands
import asyncio
import re
import io
from util import load_configs, save_configs, get_server_config, set_server_config
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw

from dotenv import load_dotenv
import os

load_dotenv("env")

TOKEN = os.getenv("TOKEN")
print(TOKEN)

pattern = re.compile(r"&([^&]+)&")

# RDKit molecule drawing setup
d2d = Draw.MolDraw2DCairo(350, 300)
opts = d2d.drawOptions()
opts.updateAtomPalette({6: (1, 1, 1)})  # Color customization
opts.setBackgroundColour((44 / 255, 45 / 255, 49 / 255))  # Dark theme

# Function to validate a SMILES string
def is_valid_smiles(smiles: str):
    try:
        return Chem.MolFromSmiles(smiles) is not None
    except:
        return False

# Function to create a molecule image
def create_molecule_image(mol):
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except:
        print("Kekulization failed, skipping.")

    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    bio = io.BytesIO(d2d.GetDrawingText())
    bio.seek(0)
    return bio

# Discord Bot Setup
intents = discord.Intents.default()
intents.message_content = True  # Required for message detection
bot = commands.Bot(command_prefix=lambda bot, message: get_prefix(bot, message), intents=intents)
tree = bot.tree  # Slash command tree

def get_prefix(bot, message):
    guild_id = str(message.guild.id)
    print(get_server_config(guild_id))
    return get_server_config(guild_id).get("prefix", "/")

@tree.command(name="help", description="Displays help menu")
async def help(interaction: discord.Interaction):
    embed = discord.Embed(title="Help Menu", description="List of available commands:", color=0x2f3136)
    embed.add_field(name="`/help`", value="Show this help menu.", inline=False)
    embed.add_field(name="`/smileshelp`", value="Quick guide to SMILES.", inline=False)
    embed.add_field(name="`/smiles <SMILES string>`", value="Generate a molecular structure image from a SMILES string.", inline=False)
    embed.add_field(name="`/auto_detect <enable|disable>`", value="Enable or disable auto-detection of `smile[...]` messages.", inline=False)
    await interaction.response.send_message(embed=embed)

@tree.command(name="smileshelp", description="Quick guide to SMILES.")
async def smileshelp(interaction: discord.Interaction):
    await interaction.response.send_message("SMILES (Simplified Molecular Input Line Entry System) is a line notation for representing molecules. Example: `CCO` for ethanol.")

@tree.command(name="setprefix", description="Set a custom prefix for the server.")
async def setprefix(interaction: discord.Interaction, new_prefix: str):
    guild_id = str(interaction.guild.id)
    
    # if not interaction.user.guild_permissions.administrator:
    #     await interaction.response.send_message("You need to be an administrator to use this command!", ephemeral=False)
    #     return

    if set_server_config(guild_id, "prefix", new_prefix):
        await interaction.response.send_message(f"The command prefix has been updated to `{new_prefix}`.", ephemeral=False)
    else:
        await interaction.response.send_message("An error occurred while setting the prefix.", ephemeral=True)


@tree.command(name="settings", description="Settings")
async def settings(interaction: discord.Interaction):
    guild_id = str(interaction.guild.id)

    # Fetch server config for the given guild
    server_config = get_server_config(guild_id)

    message = "**Settings**\n"
    for key, value in server_config.items():
        message += f"- {key}: {value}\n"
    
    # Send the settings message to the server
    await interaction.response.send_message(message)


@tree.command(name="smiles", description="Generate a molecular structure image from a SMILES string.")
async def smiles(interaction: discord.Interaction, smiles_str: str):
    if is_valid_smiles(smiles_str):
        mol = Chem.MolFromSmiles(smiles_str)
        loop = asyncio.get_running_loop()
        img = await loop.run_in_executor(None, create_molecule_image, mol)
        
        embed = discord.Embed(title=f"SMILES: `{smiles_str}`", color=0x2f3136)
        embed.set_image(url="attachment://molecule.png")
        await interaction.response.send_message(embed=embed, file=discord.File(img, filename="molecule.png"))
    else:
        await interaction.response.send_message("Invalid SMILES string. Please provide a valid format.")

@tree.command(name="auto_detect", description="Enable or disable automatic smile[...] message detection.")
async def auto_detect(interaction: discord.Interaction, option: str):
    global server_configs
    if not interaction.user.guild_permissions.administrator:
        await interaction.response.send_message("You need to be an administrator to use this command!", ephemeral=False)
        return
    
    guild_id = str(interaction.guild.id)

    if set_server_config(guild_id, "auto_detect", option):
        await interaction.response.send_message(f"Auto-detection has been set to `{option}`.", ephemeral=False)
    else:
        await interaction.response.send_message("Invalid option. Please use `true` or `false` for auto_detect.", ephemeral=True)

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
                
                # Prepare the embed and send it to the same channel where the message was detected
                embed = discord.Embed(title=f"SMILES: `{smiles_str}`", color=0x2f3136)
                embed.set_image(url="attachment://molecule.png")
                await message.channel.send(embed=embed, file=discord.File(img, filename="molecule.png"))
            else:
                await message.channel.send("Invalid SMILES string. Please provide a valid format.")

    await bot.process_commands(message)

### **Bot Ready Event** ###
@bot.event
async def on_ready():
    await tree.sync()

# Run Bot
bot.run(TOKEN)
