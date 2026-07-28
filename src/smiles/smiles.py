import asyncio
import functools
import io
import socket
from logging import Logger, getLogger
import cirpy
import discord
from rdkit import Chem
from rdkit.Chem import Draw, rdMolDescriptors, rdDepictor, rdDistGeom
from urllib.error import URLError, HTTPError

from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdChemReactions as Reactions

from constants import SMILE_BG, NAME, DISCORD_DARK_JSON
from db.db import DatabaseHandler
from util import complement_color, smile_rgb, transform_rgb_to_smile


class Smile(object):
    def __init__(self, database_handler: DatabaseHandler):
        self.db_handler:DatabaseHandler = database_handler
        self.logger:Logger = getLogger(NAME)
        self.logger.debug("Smile initialized.")
    
        self.d2d = rdMolDraw2D.MolDraw2DCairo(-1, -1) # not really being used
        self.opts = self.d2d.drawOptions()
        rdDepictor.LoadDefaultRingSystemTemplates()

    def __get_custom_keys(self, render_options: dict) -> set[str]:
        """Keys in the DB that aren't direct MolDrawOptions attributes."""
        test_opts = rdMolDraw2D.MolDraw2DCairo(-1, -1).drawOptions()
        custom = set()
        for key in render_options:
            if not hasattr(test_opts, key):
                custom.add(key)
        return custom

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

    # v2 functionality; in testing disabled for public
    # Function validates SMARTS notation before rendering
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

        # Check permission before attempting to react
        channel = sent_message.channel
        perms = channel.permissions_for(channel.guild.me)
        if perms.add_reactions:
            await sent_message.add_reaction("❌")
        else:
            self.logger.warning(f"Missing add_reactions permission in channel {channel.id}, skipping reaction.")

        self.d2d.ClearDrawing()
    
    def __addAtomNumbers(self, mol) -> None:
        for atom in mol.GetAtoms():
            i = atom.GetIdx()
            self.opts.atomLabels[i] = str(i)

    def __loadRenderOptions(self, mols, server_id) -> None:  # remove custom_keys param
        if not isinstance(mols, list):
            mols = [mols]

        bg_color = self.db_handler.render_options.get_bgcolor(server_id)
        render_options = self.db_handler.get_render_option(server_id)
        custom_keys = self.__get_custom_keys(render_options)

        complement_bg_color = tuple(c / 255 for c in complement_color(bg_color)) + (1.0,)
        bg_color = tuple(c / 255 for c in bg_color)

        self.opts.drawMolsSameScale = False
        self.d2d.SetFlexiMode(True)
        self.opts.scaleBondWidth = True
        self.opts.scaleHighlightBondWidth = True
        self.opts.legendFraction = 0.15
        self.opts.legendFontSize = 40

        self.opts.setSymbolColour(complement_bg_color)
        self.opts.setAnnotationColour(complement_bg_color)
        self.opts.setLegendColour(complement_bg_color)
        self.opts.setBackgroundColour(bg_color)
        self.opts.setHighlightColour((0, 0, 1.0, 0.1))

        if render_options.get("includeAtomNumbers"):
            for mol in mols:
                self.__addAtomNumbers(mol)

        for key, value in render_options.items():
            if key in custom_keys:
                continue
            setattr(self.opts, key, bool(value))

        self.loadAtomPalette(server_id)

    def __processLegend(self, legends, num_mols) -> tuple:
        render_legend = []
        if legends:
            legends = [legend.strip() for legend in legends.split(";")]
            for i in range(num_mols):
                if i < len(legends):
                    render_legend.append(legends[i])
                else:
                    render_legend.append("")
        else:
            render_legend = [""] * num_mols

        return tuple(render_legend)
    
    @functools.lru_cache(maxsize=1000)
    def _resolve_name_to_smiles_cached(self, name: str) -> str | None:
        return cirpy.resolve(name, "smiles")

    async def _resolve_name_to_smiles(self, ctx, name: str) -> str | None:
        try:
            return await asyncio.get_running_loop().run_in_executor(
                None, self._resolve_name_to_smiles_cached, name
            )
        except HTTPError as e:
            self.logger.warning(f"CIR HTTP error for {name!r}: {e.code} {e.reason}")
            await ctx.send(f"CIR returned an error ({e.code}): {e.reason}")
        except (URLError, socket.timeout) as e:
            self.logger.warning(f"CIR unreachable for {name!r}: {e}")
            await ctx.send(
                "The chemical name resolver (CIR) appears to be down or unreachable right now. Try again later.")
        except Exception as e:
            # Catch-all: something unexpected (bad response format, cirpy internals, etc.)
            self.logger.exception(f"Unexpected error resolving {name!r}")
            await ctx.send(f"Unexpected error resolving '{name}': {e}")
        return None

    def loadAtomPalette(self, server_id) -> None:
        palette = self.db_handler.element_colors.get_element_colors(server_id)
        if palette:
            palette = transform_rgb_to_smile(palette)
            palette = DISCORD_DARK_JSON | palette
        else:
            palette = DISCORD_DARK_JSON

        self.opts.setAtomPalette(palette)

    def create_molecule_image(self, mols, server_id, legends, **drawFuncArgs) -> io.BytesIO:
        if not isinstance(mols, list):
            mols = [mols]

        # Reinitialize drawer first
        self.d2d = rdMolDraw2D.MolDraw2DCairo(-1, -1)
        self.opts = self.d2d.drawOptions()

        for mol in mols:
            rdDepictor.Compute2DCoords(mol, useRingTemplates=True)
            Chem.SanitizeMol(mol)
            Chem.Kekulize(mol, clearAromaticFlags=True)
            rdDepictor.NormalizeDepiction(mol)
            rdDepictor.StraightenDepiction(mol)

        mols_per_row = (len(mols) + 1) // 2

        highlight_atoms = drawFuncArgs.pop("highlightAtoms", None)
        if highlight_atoms:
            highlight_atoms = tuple(int(a) for a in highlight_atoms)

        # Load DB settings onto self.opts BEFORE passing it to MolsToGridImage
        self.__loadRenderOptions(mols, server_id)

        img_data = Draw.MolsToGridImage(
            mols,
            molsPerRow=mols_per_row,
            subImgSize=(960, 540),
            legends=list(legends),
            highlightAtomLists=tuple([tuple(highlight_atoms)] * len(mols)) if highlight_atoms else None,
            returnPNG=True,
            drawOptions=self.opts,
        )

        if not isinstance(img_data, bytes):
            buf = io.BytesIO()
            img_data.save(buf, format="PNG")
            img_data = buf.getvalue()

        bio = io.BytesIO(img_data)
        bio.seek(0)
        return bio

    # v2 functionality; in testing disabled for public
    # Creates image of chemical reaction
    def create_rxn_image(self, rxn, server_id) -> io.BytesIO:
        self.d2d = rdMolDraw2D.MolDraw2DCairo(-1, -1)
        self.opts = self.d2d.drawOptions()

        self.__loadRenderOptions(rxn, server_id)
        # not sure why this is needed, but otherwise it'll error

        Reactions.SanitizeRxn(rxn)
        coords = Reactions.Compute2DCoordsForReaction(rxn)

        self.d2d.DrawReaction(rxn, confIds=coords)
        self.d2d.FinishDrawing()
        bio = io.BytesIO(self.d2d.GetDrawingText())

        bio.seek(0)
        return bio

    async def render_molecule(self, ctx, molecule, server_id, legends="", **drawFuncArgs) -> None:
        self.logger.info(f"smiles.render_molecule(ctx, {molecule}, {server_id}, {legends})")

        molecules = [m.strip() for m in molecule.split(";")]
        mol_objects = []

        if len(molecules) > 4:
            await ctx.send("You can only render 4 molecules at a time.")
            return

        for mol in molecules:
            if not self.__is_valid_smiles(mol):
                try:
                    mol = await self._resolve_name_to_smiles(ctx, mol)
                    if not mol:
                        raise ValueError("Could not resolve name to SMILES.")
                except Exception as e:
                    self.logger.error(f"Failed to resolve: {molecule}, error: {e}")
                    await ctx.send(f"Failed to resolve: {molecule}, error: {e}")
                    return

            mol_obj = Chem.MolFromSmiles(mol)

            if mol_obj:
                mol_objects.append(mol_obj)
            else:
                await ctx.send(f"Could not parse {mol}, skipping.")

        if not mol_objects:
            await ctx.send("No valid molecules to render.")
            return

        try:
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
        except Exception as e:
            self.logger.error(f"Rendering error: {e}")
            await ctx.send(f"Rendering error: {e}")
            return

    # v2 functionality; in testing disabled for public
    # Renders chemical reaction image

    async def render_reaction(self, ctx, reaction, server_id) -> None:
        self.logger.info(f"smiles.render_reaction(ctx, {reaction}, {server_id})")

        reaction = reaction.strip()
        if not self.__is_valid_smarts(reaction):
            await ctx.send(
                f"{reaction} is an invalid reaction, please check for typos/erros!"
            )
            return

        try:
            rxn = Reactions.ReactionFromSmarts(f'{reaction}', useSmiles=True)

            loop = asyncio.get_running_loop()
            img = await loop.run_in_executor(
                None,
                functools.partial(
                    self.create_rxn_image,
                    rxn,
                    server_id
                )
            )

        except Exception as e:
            self.logger.error(f"Reaction error: {e}")
            await ctx.send(f"Reaction error: {e}")
            return

        await self.__render(ctx, reaction, img)