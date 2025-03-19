from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw, rdChemReactions
import asyncio
import cirpy
import discord
import io

import src.constants
from src.constants import SMILE_BG, discord_dark

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
        
    def loadAtomPalette(self, ctx):
        if ("server_id") in self.bot.db_handler.check_color_settings(str(ctx.guild.id, "server_id")):
            p = self.bot.db_handler.get_element_colors(str(ctx.guild.id))
            p2 = discord_dark | p
        else:
            p2 = src.constants.discord_dark
        palette = self.opts.setAtomPalette(p2)
        return palette

        
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

    async def render_molecule(self, ctx, molecule, palette=None):
        molecule = molecule.strip()

        if not self.__is_valid_smiles(molecule):
            # check if molecule is identified by name
            try:
                molecule = cirpy.resolve(molecule, 'smiles')
            except:
                await ctx.send(f"{molecule} is invalid, please try with a different compound ID or check for typos/erros!")
        
        if palette:
            print(palette)
            self.loadAtomPalette(palette)

        mol = Chem.MolFromSmiles(molecule)
        loop = asyncio.get_running_loop()
        img = await loop.run_in_executor(None, self.create_molecule_image, mol)

        await self.__render(ctx, molecule, img)

    async def render_reaction(self, ctx, reaction, palette = None):
        if not self.__is_valid_smarts(reaction):
            await f"{reaction} is invalid, please try with a different compound ID or check for typos/erros!"
            return

        if palette:
            print(palette)
            self.loadAtomPalette(palette)

        rxn = Chem.ReactionFromSmarts(reaction, useSMILES = True)
        loop = asyncio.get_running_loop()
        img = await loop.run_in_executor(None, self.create, rxn)

        await self.__render(ctx, reaction, img)