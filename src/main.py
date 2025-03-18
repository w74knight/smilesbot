import discord
from discord.ext import commands

from constants import TOKEN, OWNER_ID, AUTO_DETECT_PATTERN
from cogs import __all__ as command_modules

from db.db import DatabaseHandler
from util import get_prefix
from smile.smile import Smile

class Bot(commands.Bot):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.owner_id = OWNER_ID
        self.help_command = None

        self.db_handler = DatabaseHandler()
        self.smile = Smile()

    async def on_ready(self):
        await self.change_presence(
            status=discord.Status.online, 
            activity=discord.Activity(
                    type=discord.ActivityType.watching,
                    name="/smiles"
                )
            )
        await self.tree.sync()

        print(f"Logged in as {self.user}")
        
    async def setup_hook(self):
        for module_name in command_modules:
            await self.load_extension(f'cogs.{module_name}')

    async def on_message(self, message):
        if message.author.bot:
            return
        
        settings = self.db_handler.get_server_setting(str(message.guild.id))
        if settings.get('auto_smile', False) and (match := AUTO_DETECT_PATTERN.search(message.content)):
            await self.smile.render_molecule(message.channel, match.group(1))
        
        await self.process_commands(message)

    async def close(self):
        self.db_handler.close()
        await super().close()

intents = discord.Intents.default()
intents.message_content = True
bot = Bot(command_prefix=get_prefix, intents=intents)

bot.run(TOKEN)
