import discord
import logging
import logging.handlers
from discord.ext import commands

from cogs import __all__ as command_modules
from constants import AUTO_DETECT_PATTERN, OWNER_ID, TOKEN, NAME
from db.db import DatabaseHandler
from smile.smile import Smile
from util import get_prefix

from DiscordLogger import DiscordLoggerHandler

class SmileBot(commands.Bot):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.owner_id = OWNER_ID
        self.help_command = None
        self.logger = logging.getLogger(NAME)

        self.db_handler:DatabaseHandler = None
        self.smile:Smile = None

    async def on_ready(self):
        self.logger.info("SmileBot initialized.")

        self.db_handler = DatabaseHandler()
        self.smile = Smile(self.db_handler)

        await self.change_presence(
            status=discord.Status.online, 
            activity=discord.Activity(
                    type=discord.ActivityType.watching,
                    name="/render"
                )
            )
        await self.tree.sync()

        self.logger.info(f"Logged in as {self.user} (ID: {self.user.id})")

    async def setup_hook(self):
        for module_name in command_modules:
            await self.load_extension(f'cogs.{module_name}')
        self.logger.info("All cogs loaded.")

    async def on_message(self, message):
        if message.author.bot:
            return
        
        settings = self.db_handler.get_server_setting(str(message.guild.id))
        if bool(settings.get('auto_smile', 0)) and (match := AUTO_DETECT_PATTERN.search(message.content)):
            self.logger.info(f"Auto-detecting SMILES/SMARTS in message: {message.content}")
            if ">>" in message.content:
                await self.smile.render_reaction(message.channel, match.group(1), str(message.guild.id))
            else:
                await self.smile.render_molecule(message.channel, match.group(1), str(message.guild.id))

        await self.process_commands(message)

    async def on_raw_reaction_add(self, payload):
        if payload.user_id == self.user.id:
            return

        if payload.emoji.name == "❌":
            channel = self.get_channel(payload.channel_id)
            message = await channel.fetch_message(payload.message_id)
            if message.author.id == self.user.id:
                return await message.delete()
            
    async def on_command_error(self, ctx, error):
        if isinstance(error, commands.CommandNotFound):
            return

        self.logger.error(f"Error in command {ctx.command}: {error}")

    async def close(self):
        self.db_handler.close()
        await super().close()
        self.logger.info("SmileBot closed.")


if __name__ == "__main__":
    logger = logging.getLogger(NAME)
    logger.setLevel(logging.INFO)
    logging.getLogger('discord.http').setLevel(logging.INFO)

    dt_fmt = '%Y-%m-%d %H:%M:%S'
    formatter = logging.Formatter('[{asctime}] [{levelname:<8}] {module}: {message}', dt_fmt, style='{')

    # Handlers
    file_handler = logging.handlers.RotatingFileHandler(
        filename='smilebot.log',
        encoding='utf-8',
        maxBytes=32 * 1024 * 1024,  # 32 MiB
        backupCount=5,  # Rotate through 5 files
    )
    discord_handler = DiscordLoggerHandler(level=logging.INFO)
    console_handler = logging.StreamHandler()

    # Set formatter for handlers
    file_handler.setFormatter(formatter)
    discord_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(discord_handler)
    logger.addHandler(console_handler)

    discord_logger = logging.getLogger('discord')
    discord_logger.setLevel(logging.INFO)
    discord_logger.addHandler(file_handler)
    discord_logger.addHandler(discord_handler)
    discord_logger.addHandler(console_handler)

    intents = discord.Intents.default()
    intents.message_content = True
    smilebot = SmileBot(command_prefix=get_prefix, intents=intents)

    smilebot.run(TOKEN, log_handler=None)
