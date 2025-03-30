from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdChemReactions as Reactions
from rdkit.Chem import Draw
import asyncio
import cirpy
import discord
import io
from constants import SMILE_BG, smile_rgb
from .pallette import DISCORD_DARK
import functools
from util import transform_rgb_to_smile

from db.db import DatabaseHandler

class Smile(object):
    def __init__(self, database_handler: DatabaseHandler):
        self.db_handler:DatabaseHandler = database_handler

        self.d2d = Draw.MolDraw2DCairo(-1, -1)
        self.opts = self.d2d.drawOptions()
        self.opts.setBackgroundColour(smile_rgb(*SMILE_BG))
        self.opts.drawMolsSameScale = False

        self.opts.scalingFactor = 50
        self.opts.fixedFontSize = 20
        self.opts.bondLineWidth = 2.

        # rxn options
        self.opts.setSymbolColour((1, 1, 1))
        self.opts.setAnnotationColour((1,1,1))

    def __is_valid_smiles(self, smiles: str):
        try:
            return Chem.MolFromSmiles(smiles) is not None
        except:
            return False

    def __is_valid_smarts(self, rxn: str):
        try:
            return Chem.ReactionFromSmarts(rxn, useSmiles=True) is not None
        except:
            return False

    async def __render(self, ctx, mlcl, img):
        embed = discord.Embed(
            title=f'Render Complete!',
            description=f"{mlcl}"
        )
        embed.set_image(url="attachment://molecule.png")

        sent_message = await ctx.send(embed=embed, file=discord.File(img, filename="molecule.png"))
        
        await sent_message.add_reaction("❌")

        self.d2d.ClearDrawing()

    def __draw(self, drawFunc, mol, server_id, **drawFuncArgs):
        bg_color = self.db_handler.render_options.get_bgcolor(server_id)
        render_options = self.db_handler.get_render_option(server_id)

        # convert rgb to ratio
        bg_color = tuple(c / 255 for c in bg_color)
        
        self.opts.setBackgroundColour(bg_color)
        self.opts.setHighlightColour((0, 0, 1.0, 0.1))
        for key, value in render_options.items():
            setattr(self.opts, key, bool(value))

        self.loadAtomPalette(server_id)

        drawFunc(mol, **drawFuncArgs)

        self.d2d.FinishDrawing()

        bio = io.BytesIO(self.d2d.GetDrawingText())
        bio.seek(0)
        return bio

    def loadAtomPalette(self, server_id):
        pallette = self.db_handler.element_colors.get_element_colors(server_id)
        if pallette:
            pallette = transform_rgb_to_smile(pallette)
            pallette = DISCORD_DARK | pallette
        else:
            pallette = DISCORD_DARK

        self.opts.setAtomPalette(pallette)

    def create_molecule_image(self, mol, server_id, **drawFuncArgs):
        # not sure why this is needed, but otherwise it'll error
        self.d2d = Draw.MolDraw2DCairo(-1, -1)
        self.opts = self.d2d.drawOptions()

        try:
            Chem.Kekulize(mol, clearAromaticFlags=True)
        except:
            print("Kekulization failed, skipping.")
    
        return self.__draw(self.d2d.DrawMolecule, mol, server_id, **drawFuncArgs)

    def create_rxn_image(self, rxn, server_id, **drawFuncArgs):
        # not sure why this is needed, but otherwise it'll error
        self.d2d = Draw.MolDraw2DCairo(-1, -1)
        self.opts = self.d2d.drawOptions()


        return self.__draw(self.d2d.DrawReaction, rxn, server_id, **drawFuncArgs)

    async def render_molecule(self, ctx, molecule, server_id, **drawFuncArgs):
        molecule = molecule.strip()

        if not self.__is_valid_smiles(molecule):
            # check if molecule is identified by name
            try:
                molecule = cirpy.resolve(molecule, 'smiles')
            except:
                await ctx.send(
                    f"{molecule} is invalid, please try with a different compound ID or check for typos/erros!")

        mol = Chem.MolFromSmiles(molecule)
        loop = asyncio.get_running_loop()
        img = await loop.run_in_executor(
            None,
            functools.partial(
                self.create_molecule_image,
                mol,
                server_id, **drawFuncArgs
            )
        )

        await self.__render(ctx, molecule, img)

    async def render_reaction(self, ctx, reaction, server_id, **drawFuncArgs):
        reaction = reaction.strip()
        if not self.__is_valid_smarts(reaction):
            await ctx.send(
                f"{reaction} is an invalid reaction, please check for typos/erros!"
            )

        rxn = Reactions.ReactionFromSmarts(f'{reaction}', useSmiles=True)
        loop = asyncio.get_running_loop()
        img = await loop.run_in_executor(
            None,
            functools.partial(
                self.create_rxn_image,
                rxn,
                server_id, **drawFuncArgs
            )
        )

        await self.__render(ctx, reaction, img)