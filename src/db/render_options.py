import sqlite3
from logging import getLogger

from constants import NAME, SMILE_BG

DEFAULT_RENDER_OPTIONS = {
    "includeAtomNumbers": False,
    "addStereoAnnotations": False,
    "explicitMethyl": False,
    "atomLabelDeuteriumTritium": False
}

class RenderOptions:
    def __init__(self, connection):
        self.logger = getLogger(NAME)
        self.connection = connection
        self.cursor = connection.cursor()
        self.init()
        self.logger.info("RenderOptions initialized.")

    def init(self):
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS render_options (
                server_id TEXT PRIMARY KEY,
                background_color TEXT,
                includeAtomNumbers BOOLEAN,
                addStereoAnnotations BOOLEAN,
                explicitMethyl BOOLEAN,
                atomLabelDeuteriumTritium BOOLEAN,
                dummiesAreAttachments BOOLEAN
            )
        ''')
        self.connection.commit()

    # Get render options
    def get_render_option(self, server_id):
        self.cursor.execute('''
            SELECT * FROM render_options WHERE server_id=?
        ''', (server_id,))
        result = self.cursor.fetchone()
        return dict(result) if result else DEFAULT_RENDER_OPTIONS

    # Setter/Getter for background color
    def set_bgcolor(self, server_id, color):
        self.cursor.execute('''
            INSERT INTO render_options (server_id, background_color)
            VALUES (?, ?)
            ON CONFLICT(server_id) DO UPDATE SET
            background_color=excluded.background_color
        ''', (server_id, color))
        self.connection.commit()

    def get_bgcolor(self, server_id):
        self.cursor.execute('''
            SELECT background_color FROM render_options WHERE server_id=?
        ''', (server_id,))
        result = self.cursor.fetchone()
        
        if result and result["background_color"]:
            bg_color = result["background_color"].split(",")
            return tuple(int(c) for c in bg_color)

        return SMILE_BG

    def set(self, server_id, key, value):
        query = f'''
            INSERT INTO render_options (server_id, {key})
            VALUES (?, ?)
            ON CONFLICT(server_id) DO UPDATE SET
            {key}=excluded.{key}
        '''
        self.cursor.execute(query, (server_id, value))
        self.connection.commit()

    def get(self, server_id, key):
        self.cursor.execute(f'''SELECT {key} FROM render_options WHERE server_id=?''', (server_id,))
        result = self.cursor.fetchone()
        return result[key] if result else False
    
    def clear(self, server_id):
        self.cursor.execute('''
        DELETE FROM render_options WHERE server_id = ?;
        ''', (server_id,))
        self.connection.commit()