import discord
from discord.ext import commands
from rdkit import Chem
from rdkit.Chem import Draw
import io
from PIL import Image

# Set up the bot with the necessary intents
intents = discord.Intents.default()
intents.message_content = True  # Enable the message content intent

bot = commands.Bot(command_prefix="!", intents=intents)


# Function to check if a string is a valid SMILES
def is_valid_smiles(smiles: str):
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except:
        return False


# Function to handle kekulization and image creation
def create_molecule_image(mol):
    try:
        Chem.Kekulize(mol)
    except:
        print("Kekulization failed, skipping.")

        # Create an image from the molecule
    img = Draw.MolToImage(mol, size=(900, 900))

    # Convert the image to RGBA (which supports transparency)
    img = img.convert("RGBA")
    datas = img.getdata()

    # Replace white pixels with transparent pixels (optional)
    new_data = []
    for item in datas:
        if item[:3] == (255, 255, 255):  # white color
            new_data.append((255, 255, 255, 0))  # transparent
        else:
            new_data.append((255, 255, 255, item[3]))  # Keep original alpha
    img.putdata(new_data)

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
        await ctx.send("Invalid SMILES string. Please provide a valid SMILES format.")

# Run the bot with your token
bot.run('MTM0NzI3NTQxNjIyNTk3NjMzMw.GG4lPZ.voqxNW2Qw-MmiKdFAVkawlx_dPs2teUOKEgoB4')
