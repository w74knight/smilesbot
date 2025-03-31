import sqlite3
from logging import Logger, getLogger

from constants import NAME

from .element_colors import ElementColors
from .render_options import RenderOptions
from .server_settings import ServerSettings


class DatabaseHandler:
    def __init__(self, db_name='database.db'):
        self.logger:Logger = getLogger(NAME)
        self.connection:sqlite3.Connection = sqlite3.connect(db_name, check_same_thread=False)
        self.connection.row_factory = sqlite3.Row
        self.cursor = self.connection.cursor()
        
        self.render_options:RenderOptions = RenderOptions(self.connection)
        self.element_colors:ElementColors = ElementColors(self.connection)
        self.server_settings:ServerSettings = ServerSettings(self.connection)

        self.logger.info("DatabaseHandler initialized.")

    def clear(self, server_id) -> None:
        """
        Clear all settings for a given server ID.
        This includes render options, element colors, and server settings.
        """

        self.render_options.clear(server_id)
        self.element_colors.clear(server_id)
        self.server_settings.clear(server_id)

    def get_render_option(self, server_id) -> dict:
        """
        Get render options for a given server ID.
        """

        return self.render_options.get_render_option(server_id)
    
    def get_element_colors(self, server_id) -> dict:
        """
        Get element colors for a given server ID.
        """

        return self.element_colors.get_element_colors(server_id)
    
    def get_server_setting(self, server_id) -> dict:
        """
        Get server settings for a given server ID.
        """

        return self.server_settings.get_server_setting(server_id)
    
    def close(self) -> None:
        '''
        Close the database connection.
        '''
        
        self.logger.info("Closing database connection.")
        self.connection.close()