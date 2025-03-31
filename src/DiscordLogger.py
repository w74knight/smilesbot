import logging
import aiohttp
import discord
from constants import WEBHOOK_URL
from discord.ext import tasks
import threading

class DiscordLoggerHandler(logging.Handler):
    def __init__(self, bot, level=logging.NOTSET):
        super().__init__(level)
        self.webhook_url = WEBHOOK_URL
        self.log_queue = []
        self.bot = bot
        self.max_queue_size = 1000

    async def send_log(self, message, session):
        webhook = discord.Webhook.from_url(self.webhook_url, session=session)
        try:
            await webhook.send(message)
        except Exception as e:
            print(f"Failed to send log: {e}")

    def emit(self, record):
        try:
            log_message = self.format(record)
            if len(self.log_queue) >= self.max_queue_size:
                self.log_queue.pop(0)
            self.log_queue.append(log_message)

            if not self.send_logs.is_running():
                self.send_logs.start()
        except Exception as e:
            print(f"Log emission failed: {e}")
            self.handleError(record)

    @tasks.loop(seconds=10.0)
    async def send_logs(self):
        try:
            if not self.log_queue:
                return

            async with aiohttp.ClientSession() as session:
                while self.log_queue:  # Process until queue is empty
                    batch = self._get_batch()
                    if batch:  # In case all messages were too long
                        await self.send_log(batch, session)
        except Exception as e:
            print(f"Error in send_logs task: {e}")

    @send_logs.before_loop
    async def before_send_logs(self):
        await self.bot.wait_until_ready()

    def _get_batch(self):
        batch = ""
        while self.log_queue:
            next_log = self.log_queue[0]
            if len(batch) + len(next_log) + 1 > 2000:
                break
            batch += self.log_queue.pop(0) + "\n"
        return batch.strip()

    def close(self):
        self.send_logs.cancel()
        super().close()
