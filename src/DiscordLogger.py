import logging
import discord
import aiohttp
from constants import WEBHOOK_URL
import asyncio

class DiscordLoggerHandler(logging.Handler):
    def __init__(self, level=logging.NOTSET):
        super().__init__(level)
        self.webhook_url = WEBHOOK_URL

    async def send_log(self, message):
        async with aiohttp.ClientSession() as session:
            webhook = discord.Webhook.from_url(self.webhook_url, session=session)
            await webhook.send(message)

    def emit(self, record):
        """Format the log record and send it asynchronously."""
        try:
            # Get the log message in a readable format
            log_message = self.format(record)

            # Use the running event loop if available, or create a new one
            loop = asyncio.get_event_loop()

            if loop.is_running():
                # Use create_task to run the async function without blocking the event loop
                loop.create_task(self.send_log(log_message))
            else:
                # If no event loop is running, start one
                asyncio.run(self.send_log(log_message))
        except Exception as e:
            self.handleError(record)