import discord
from discord.ext import commands
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
import io
import os
from dotenv import load_dotenv

load_dotenv()

# Set up bot intents
intents = discord.Intents.default()
intents.message_content = True  # Enable message content intent
bot = commands.Bot(command_prefix="!", intents=intents, help_command=None)


# Function to check if a string is a valid SMILES
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

    img = Draw.MolToImage(mol, size=(350, 350), kekulize=True)
    return img


@bot.command()
async def smiles(ctx, *, smiles_str: str):
    if is_valid_smiles(smiles_str):
        try:
            mol = Chem.MolFromSmiles(smiles_str)
            img = create_molecule_image(mol)

            with io.BytesIO() as img_binary:
                img.save(img_binary, format="PNG")
                img_binary.seek(0)

                embed = discord.Embed(title=f"{ctx.author.name}", color=0x2f3136)
                embed.set_image(url="attachment://molecule.png")

                await ctx.send(embed=embed, file=discord.File(img_binary, filename="molecule.png"))
        except Exception as e:
            await ctx.send(f"An error occurred while processing the SMILES string: {str(e)}")
    else:
        await ctx.send("Invalid SMILES string. Please provide a valid format.")


@bot.command()
async def help(ctx):
    embed = discord.Embed(title="Help Menu", description="List of available commands:", color=0x2f3136)
    embed.add_field(name="`!smiles <SMILES string>`",
                    value="Generate a molecular structure image from a SMILES string.", inline=False)
    embed.add_field(name="`!help`",
                    value="Show this help menu.", inline=False)
    await ctx.send(embed=embed)

#token
bot.run(os.getenv("TOKEN"))
