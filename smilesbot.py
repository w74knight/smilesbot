import discord
from discord.ext import commands
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
import asyncio
import io
import os
from dotenv import load_dotenv

load_dotenv()

d2d = Draw.MolDraw2DCairo(350, 300)
opts = d2d.drawOptions()
opts.updateAtomPalette({6: (1, 1, 1)})
opts.setBackgroundColour((44 / 255, 45 / 255, 49 / 255))

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

    img = d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    bio = io.BytesIO(d2d.GetDrawingText())
    return bio

@bot.command()
async def smiles(ctx, *, smiles_str: str, none=None):
    if is_valid_smiles(smiles_str):
        try:

            mol = Chem.MolFromSmiles(smiles_str)
            loop = asyncio.get_running_loop()
            img = await loop.run_in_executor(None, create_molecule_image, mol)
            img.seek(0)

            embed = discord.Embed(title=f"{ctx.author.name}", color=0x2f3136)
            embed.set_image(url="attachment://molecule.png")
            await ctx.send(embed=embed, file=discord.File(create_molecule_image(mol), filename="molecule.png"))

        except Exception as e:
            await ctx.send(f"An error occurred while processing the SMILES string: {str(e)}")
        finally:
            d2d.ClearDrawing()
    else:
        await ctx.send("Invalid SMILES string. Please provide a valid format.")


@bot.command()
async def help(ctx):
    embed = discord.Embed(title="Help Menu", description="List of available commands:", color=0x2f3136)
    embed.add_field(name="`!help`",
                    value="Show this help menu.", inline=False)
    embed.add_field(name="`!smileshelp`",
                    value="Quick guide to SMILES.", inline=False)
    embed.add_field(name="`!smiles <SMILES string>`",
                    value="Generate a molecular structure image from a SMILES string.", inline=False)
    await ctx.send(embed=embed)

@bot.command()
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

#token
bot.run(os.getenv("TOKEN"))
