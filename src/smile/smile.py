import asyncio
import functools
import io
from logging import Logger, getLogger

import cirpy
import discord
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdChemReactions as Reactions

from constants import SMILE_BG, NAME
from db.db import DatabaseHandler
from util import complement_color, smile_rgb, transform_rgb_to_smile

from .pallette import DISCORD_DARK


class Smile(object):
    def __init__(self, database_handler: DatabaseHandler):
        self.db_handler:DatabaseHandler = database_handler
        self.logger:Logger = getLogger(NAME)
        self.logger.debug("Smile initialized.")

        self.d2d = Draw.MolDraw2DCairo(-1, -1)
        self.opts = self.d2d.drawOptions()

    def __is_valid_smiles(self, smiles: str) -> bool:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False
            Chem.SanitizeMol(mol)
            return True
        except Exception as e:
            self.logger.error(f"SMILES validation failed: {e}")
            return False

    def __is_valid_smarts(self, rxn: str) -> bool:
        try:
            rxn = Chem.ReactionFromSmarts(rxn, useSmiles=True)
            return rxn is not None
        except Exception as e:
            self.logger.error(f"SMILES validation failed: {e}")
            return False

    async def __render(self, ctx, mlcl, img) -> None:
        embed = discord.Embed(
            title=f'Render Complete!',
            description=f"{mlcl}"
        )
        embed.set_image(url="attachment://molecule.png")


        sent_message = await ctx.send(embed=embed, file=discord.File(img, filename="molecule.png"))
        
        await sent_message.add_reaction("❌")

        self.d2d.ClearDrawing()
    
    def __addAtomNumbers(self, mol) -> None:
        for atom in mol.GetAtoms():
            i = atom.GetIdx()
            self.opts.atomLabels[i] = str(i)

    def __loadRenderOptions(self, mols, server_id) -> None:
        if not isinstance(mols, list):
            mols = [mols]

        bg_color = self.db_handler.render_options.get_bgcolor(server_id)
        render_options = self.db_handler.get_render_option(server_id)

        # complementary color for legend & reaction plus and arrows
        complement_bg_color = complement_color(bg_color) + (1.0,)
        self.opts.setLegendColour(complement_bg_color)
        self.opts.setSymbolColour(complement_bg_color)

        # convert rgb to ratio
        bg_color = tuple(c / 255 for c in bg_color)
        self.opts.setBackgroundColour(bg_color)
        self.opts.setHighlightColour((0, 0, 1.0, 0.1))

        self.opts.setBackgroundColour(smile_rgb(*SMILE_BG))
        self.opts.drawMolsSameScale = False

        self.opts.scalingFactor = 50
        self.opts.fixedFontSize = 20
        self.opts.bondLineWidth = 2.

        # rxn options
        self.opts.setSymbolColour((1, 1, 1))
        self.opts.setAnnotationColour((1, 1, 1))

        if (render_options.get("includeAtomNumbers")):
            for mol in mols:
                self.__addAtomNumbers(mol)
            del render_options["includeAtomNumbers"]

        for key, value in render_options.items():
            setattr(self.opts, key, bool(value))

        self.loadAtomPalette(server_id)

    def __draw(self, drawFunc, mol, server_id, **drawFuncArgs) -> io.BytesIO:
        self.__loadRenderOptions(mol, server_id)

        drawFunc(mol, **drawFuncArgs)
        self.d2d.FinishDrawing()

        bio = io.BytesIO(self.d2d.GetDrawingText())
        bio.seek(0)
        return bio
    
    def __processLegend(self, legends, num_mols) -> list:
        render_legend = []
        if legends:
            legends = [legend.strip() for legend in legends.split(",")]
            for i in range(num_mols):
                if i < len(legends):
                    render_legend.append(legends[i])
                else:
                    render_legend.append("")
        else:
            render_legend = [""] * num_mols

        return render_legend
    
    @functools.lru_cache(maxsize=1000)
    def _resolve_name_to_smiles(self, name: str) -> str:
        return cirpy.resolve(name, 'smiles')

    def loadAtomPalette(self, server_id) -> None:
        pallette = self.db_handler.element_colors.get_element_colors(server_id)
        if pallette:
            pallette = transform_rgb_to_smile(pallette)
            pallette = DISCORD_DARK | pallette
        else:
            pallette = DISCORD_DARK

        self.opts.setAtomPalette(pallette)

    def create_molecule_image(self, mols, server_id, legends, **drawFuncArgs) -> io.BytesIO:
        if not isinstance(mols, list):
            mols = [mols]

        for mol in mols:
            try:
                Chem.Kekulize(mol, clearAromaticFlags=True)
            except:
                print("Kekulization failed, skipping.")

        self.d2d = Draw.MolDraw2DCairo(-1, -1)
        self.opts = self.d2d.drawOptions()

        mols_per_row = (len(mols) + 1) // 2

        highlight_atoms = drawFuncArgs.pop("highlightAtoms", None)
        self.__loadRenderOptions(mols, server_id)

        img_data = Draw.MolsToGridImage(
            mols,
            **drawFuncArgs,
            # **self.opts.__dict__,
            subImgSize=(300, 300),
            molsPerRow=mols_per_row,
            legends=legends,
            highlightAtomLists=[highlight_atoms] * len(mols) if highlight_atoms else None,
            returnPNG=True,
            drawOptions=self.opts,
        )

        bio = io.BytesIO(img_data)
        bio.seek(0)
        return bio

    def create_rxn_image(self, rxn, server_id, **drawFuncArgs) -> io.BytesIO:
        # not sure why this is needed, but otherwise it'll error
        self.d2d = Draw.MolDraw2DCairo(-1, -1)
        self.opts = self.d2d.drawOptions()


        return self.__draw(self.d2d.DrawReaction, rxn, server_id, **drawFuncArgs)

    async def render_molecule(self, ctx, molecule, server_id, legends, **drawFuncArgs) -> None:
        self.logger.info(f"smile.render_molecule(ctx, {molecule}, {server_id}, {legends})")

        molecules = [m.strip() for m in molecule.split(",")]
        mol_objects = []

        if len(molecules) > 4:
            await ctx.send("You can only render 4 molecules at a time.")
            return
        
        for mol in molecules:
            if not self.__is_valid_smiles(mol):
                try:
                    mol = self._resolve_name_to_smiles(mol)
                except:
                    await ctx.send(f"{mol} is invalid, please try another compound ID.")
                    return

            mol_obj = Chem.MolFromSmiles(mol)
            if mol_obj:
                mol_objects.append(mol_obj)
            else:
                await ctx.send(f"Could not parse {mol}, skipping.")

        if not mol_objects:
            await ctx.send("No valid molecules to render.")
            return

        loop = asyncio.get_running_loop()
        img = await loop.run_in_executor(
            None,
            functools.partial(
                self.create_molecule_image,
                mol_objects,
                server_id,
                legends=self.__processLegend(legends, len(mol_objects)),
                **drawFuncArgs
            )
        )

        await self.__render(ctx, ", ".join(molecules), img)


    async def render_reaction(self, ctx, reaction, server_id) -> None:
        self.logger.info(f"smile.render_reaction(ctx, {reaction}, {server_id})")

        reaction = reaction.strip()
        if not self.__is_valid_smarts(reaction):
            await ctx.send(
                f"{reaction} is an invalid reaction, please check for typos/erros!"
            )

        try:
            rxn = Reactions.ReactionFromSmarts(f'{reaction}', useSmiles=True)
        except Exception as e:
            self.logger.error(f"Reaction error: {e}")
            await ctx.send("Reaction error: {e}")
            return
        
        loop = asyncio.get_running_loop()
        img = await loop.run_in_executor(
            None,
            functools.partial(
                self.create_rxn_image,
                rxn,
                server_id
            )
        )

        await self.__render(ctx, reaction, img)