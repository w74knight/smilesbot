from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw, rdChemReactions
import asyncio
import cirpy
import discord
import io

from constants import SMILE_BG
from .pallette import DISCORD_DARK

class Smile(object):
    def __init__(self):
        self.d2d = Draw.MolDraw2DCairo(500, 500)
        self.opts = self.d2d.drawOptions()
        self.opts.setBackgroundColour(SMILE_BG)

    def __is_valid_smiles(self, smiles: str):
        try:
            return Chem.MolFromSmiles(smiles) is not None
        except:
            return False
    
    def __is_valid_smarts(self, smarts: str):
        try:
            return rdChemReactions.ReactionFromSmarts(smarts, useSMILES = True) is not None
        except:
            return False


    def loadAtomPalette(self, pallette):
        self.opts.setAtomPalette(pallette)

    async def __render(self, ctx, title, img):
        embed = discord.Embed(title=f"`{title}`")
        embed.set_image(url="attachment://molecule.png")
        await ctx.send(embed=embed, file=discord.File(img, filename="molecule.png"))
        self.d2d.ClearDrawing()

    def create_molecule_image(self, mol):
        try:
            Chem.Kekulize(mol, clearAromaticFlags=True)
        except:
            print("Kekulization failed, skipping.")
        self.opts.bondLineWidth = 2.
        self.opts.setBackgroundColour(SMILE_BG)
        self.d2d.DrawMolecule(mol)
        self.d2d.FinishDrawing()
        bio = io.BytesIO(self.d2d.GetDrawingText())
        bio.seek(0)
        return bio

    def create_rxn_image(self, mol):
        try:
            Chem.Kekulize(mol, clearAromaticFlags=True)
        except:
            print("Kekulization failed, skipping.")
        self.opts.bondLineWidth = 2.
        self.opts.setBackgroundColour(SMILE_BG)
        self.d2d.DrawMolecule(mol)
        self.d2d.FinishDrawing()
        bio = io.BytesIO(self.d2d.GetDrawingText())
        bio.seek(0)
        return bio

    async def render_molecule(self, ctx, molecule, palette):
        molecule = molecule.strip()

        if not self.__is_valid_smiles(molecule):
            # check if molecule is identified by name
            try:
                molecule = cirpy.resolve(molecule, 'smiles')
            except:
                await ctx.send(f"{molecule} is invalid, please try with a different compound ID or check for typos/erros!")

        if palette:
            palette = DISCORD_DARK | palette
        else:
            palette = DISCORD_DARK

        self.loadAtomPalette(palette)

        mol = Chem.MolFromSmiles(molecule)
        loop = asyncio.get_running_loop()
        img = await loop.run_in_executor(None, self.create_molecule_image, mol)

        await self.__render(ctx, molecule, img)

    async def render_reaction(self, ctx, reaction, palette):
        if not self.__is_valid_smarts(reaction):
            await f"{reaction} is invalid, please try with a different compound ID or check for typos/erros!"
            return

        rxn = Chem.ReactionFromSmarts(reaction, useSMILES = True)
        loop = asyncio.get_running_loop()
        img = await loop.run_in_executor(None, self.create, rxn)

        await self.__render(ctx, reaction, img)