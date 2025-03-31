import sqlite3

from .element_colors import ElementColors
from .render_options import RenderOptions
from .server_settings import ServerSettings


class DatabaseHandler:
    def __init__(self, db_name='database.db'):
        self.connection = sqlite3.connect(db_name, check_same_thread=False)
        self.connection.row_factory = sqlite3.Row
        self.cursor = self.connection.cursor()
        
        self.render_options = RenderOptions(self.connection)
        self.element_colors = ElementColors(self.connection)
        self.server_settings = ServerSettings(self.connection)

    def clear(self, server_id):
        self.render_options.clear(server_id)
        self.element_colors.clear(server_id)
        self.server_settings.clear(server_id)

    def get_render_option(self, server_id):
        return self.render_options.get_render_option(server_id)
    
    def get_element_colors(self, server_id):
        return self.element_colors.get_element_colors(server_id)
    
    def get_server_setting(self, server_id):
        return self.server_settings.get_server_setting(server_id)

    def close(self):
        self.connection.close()